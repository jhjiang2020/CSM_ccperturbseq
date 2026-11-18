#!/usr/bin/env python3
"""Versioned four-library assignment benchmark; never writes production GEX."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
LIBS = ['C1', 'C2', 'T1', 'T2']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')


def run(command, **kwargs):
    print('RUN', ' '.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True, **kwargs)


def prepare(core):
    import yaml
    import pandas as pd
    import anndata as ad
    if (core/'prepared/_SUCCESS').exists():
        raise FileExistsError('Benchmark already prepared')
    core.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((ROOT/'config.yaml').read_text())
    assert cfg['gex_alignment']['backend'] == 'star'
    assert cfg['gex_alignment']['gene_assignment'] == 'unique'
    reference = pd.read_csv(cfg['paths']['grna_reference']).fillna('')
    reference['original_grna_name'] = reference.grna_name
    reference['grna_name'] = [
        'g_' + hashlib.sha256(f'{r.grna_seq}:{r.position}:{r.target_gene}'.encode()).hexdigest()[:16]
        for r in reference.itertuples()]
    refpath = core/'reference.csv'
    reference.to_csv(refpath, index=False)
    manifest = json.loads((ROOT/'output/pilot_4_sublib_audit/work/cohort_manifest.json').read_text())
    cells = []
    for lib in LIBS:
        m = manifest['libraries'][lib]
        cells += [dict(cell_id=f'{m["global_index"]}_{bc}', base_barcode=bc,
                       sublibrary=lib, sample='Control' if lib.startswith('C') else 'Treatment')
                  for bc in m['gex_good_base_barcodes']]
    cells = pd.DataFrame(cells).sort_values('cell_id')
    a = ad.read_h5ad(ROOT/'output/combined_adata.h5ad', backed='r')
    assert cells.cell_id.is_unique and set(cells.cell_id) <= set(a.obs_names)
    cells.to_csv(core/'cells.csv.gz', index=False)
    a.file.close()
    cfg['paths']['grna_reference'] = str(refpath)
    cfg['paths']['log_dir'] = str(core/'jobs')
    cfg.setdefault('guide_design', {})['enabled'] = True
    cfg['barcode_extraction']['step1a_process_template'] = str(ROOT/'output/gex_{library}/process')
    cfg['guide_assignment']['compute_saturation_curve'] = False
    cfg['grna_structure'].update(relaxed_structure='auto', relaxed_max_5prime_loss=3,
        relaxed_min_retained=17, relaxed_max_tso_shift=3, relaxed_min_quality=30,
        relaxed_require_ibar=False)
    cfg.pop('grna_recovery', None)
    d = cfg['gex_decontamination']
    d['enabled'] = True
    d['existing_alignment']['enabled'] = False
    d['step1c']['enabled'] = True
    d['combine_with_grna_arm'].update(enabled=True, combine_before_umi_collapse=True,
                                     retain_source_provenance=True)
    for arm in ['recovered', 'standard']:
        c = copy.deepcopy(cfg)
        dest = core/arm
        dest.mkdir(exist_ok=True)
        c['paths']['output_dir'] = str(dest)
        if arm == 'standard':
            c['grna_structure']['relaxed_structure'] = False
            c['gex_decontamination']['step1c']['enabled'] = False
            c['gex_decontamination']['combine_with_grna_arm']['enabled'] = False
            (dest/'grna_barcodes').symlink_to(core/'recovered/grna_barcodes', target_is_directory=True)
        (core/f'{arm}.yaml').write_text(yaml.safe_dump(c, sort_keys=False))
    inputs = [ROOT/'output/combined_adata.h5ad', ROOT/'config.yaml',
              Path(cfg['paths']['grna_reference'])]
    for lib in LIBS:
        inputs += [ROOT/f'output/gex_{lib}/process/{name}'
                   for name in ['barcode_data.csv','run_proc_def.json']]
    save(core/'prepared/manifest.json', dict(cells=len(cells), libraries=LIBS,
         guide_units=reference.grna_name.nunique(), pairs=reference.pair_id.nunique(),
         expression_changed=False, expression_backend='Parse pre + STAR + unique',
         inputs=[dict(path=str(p),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,
                      sha256=digest(p) if p.suffix!='.h5ad' else None) for p in inputs],
         sources={str(p.relative_to(ROOT)):digest(p) for p in (ROOT/'scripts').glob('*.py')},
         reference_name_policy='designed spacer + position + target; aliases retained'))
    (core/'prepared/_SUCCESS').touch()


def capture(core, library):
    """Full Parse preprocessing to scratch; only guide evidence enters benchmark."""
    import yaml
    import pandas as pd
    assert (core/'prepared/_SUCCESS').exists()
    cfgpath = core/'recovered.yaml'
    cfg = yaml.safe_load(cfgpath.read_text())
    lib = next(x for x in cfg['libraries'] if x['sublibrary_name'] == library)
    fqroot = Path(cfg['paths']['fastq_dir'])
    first = sorted(fqroot.rglob(f'*{lib["gex_prefix"]}*_R1_*fastq.gz'))
    second = [Path(str(p).replace('_R1_', '_R2_')) for p in first]
    assert first and all(p.exists() for p in second)
    dest = core/'recovered/gex_grna_filter'/library
    dest.mkdir(parents=True, exist_ok=True)
    assert not (dest/'_SUCCESS').exists()
    identity = dest/'step1a_input_identity.json'
    run([sys.executable,ROOT/'scripts/write_fastq_input_identity.py','--library',library,
         '--output',identity,*sum((['--r1',p] for p in first),[]),
         *sum((['--r2',p] for p in second),[])])
    scratch = Path('/broad/hptmp/jj')/core.name/f'gex_pre_{library}'
    assert not scratch.exists(), f'Refusing existing preprocessing scratch: {scratch}'
    env = dict(os.environ, CSM_GEX_DECONTAM_CONFIG=str(cfgpath),
        CSM_GEX_DECONTAM_LIBRARY=library, CSM_GEX_DECONTAM_OUTPUT_DIR=str(dest),
        CSM_GEX_DECONTAM_INPUT_MANIFEST=str(identity),
        CSM_GEX_DECONTAM_HOOK=str(ROOT/'scripts/gex_grna_preprocess_hook.py'),
        CSM_GEX_DECONTAM_RUN_ID=f'assign.{os.environ.get("JOB_ID","local")}.{library}')
    run([sys.executable,ROOT/'scripts/splitpipe_gex_wrapper.py',
         '--wrapper_splitpipe_dir',cfg['paths']['splitpipe_dir'],
         '--wrapper_patch',ROOT/'patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch',
         '--mode','pre','--one_step','--kit',cfg['splitpipe']['kit'],
         '--chemistry',cfg['splitpipe']['chemistry'],'--genome_dir',cfg['paths']['genome_dir'],
         '--fq1',*first,'--fq2',*second,'--output_dir',scratch,
         '--sample',lib['sample_name'],cfg['splitpipe']['sample_wells'],
         '--nthreads',os.environ.get('NSLOTS','1'),'--kit_score_skip'], env=env)
    # Coordinate identity must agree with the unchanged expression arm.
    from barcode_contract import _normalize_table
    old = _normalize_table(pd.read_csv(ROOT/f'output/gex_{library}/process/barcode_data.csv',comment='#'))
    new = _normalize_table(pd.read_csv(scratch/'process/barcode_data.csv',comment='#'))
    joined = old.merge(new, on=['bc_round','sequence'],suffixes=('_old','_new'),validate='one_to_one')
    assert len(joined)==len(old)
    for col in ['bci','wind','well']:
        assert (joined[col+'_old']==joined[col+'_new']).all(), col
    import shutil
    for name in ['barcode_data.csv','run_proc_def.json']:
        shutil.copy2(scratch/'process'/name, dest/name)
    assert (dest/'rejected_pairs.complete.json').exists()
    save(dest/'benchmark_capture.json',dict(scope='all raw GEX reads; Parse pre only',
         scratch=str(scratch),canonical_barcode_rows=len(old),runtime_barcode_rows=len(new),
         expression_changed=False))
    (dest/'_SUCCESS').touch()


def validate_standard(core):
    """Verify the completed standard product before a recovered-only retry."""
    import yaml
    import relaxed_structure
    cfg = yaml.safe_load((core/'standard.yaml').read_text())
    directory = core/'standard/grna_assignment'
    assert (directory/'grna_dedicated.complete.json').is_file()
    relaxed_structure.validate_assignment(cfg, directory, 'dedicated')
    print('Completed standard counts and their input/config identities verified')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','capture','validate-standard'])
    p.add_argument('--core',required=True,type=Path)
    p.add_argument('--library',choices=LIBS)
    a=p.parse_args()
    if a.stage=='prepare': prepare(a.core.resolve())
    elif a.stage=='validate-standard': validate_standard(a.core.resolve())
    else: capture(a.core.resolve(),a.library)


if __name__=='__main__': main()
