#!/usr/bin/env python3
"""Guide-level adapters for all eleven crispat families and CSM comparators."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import traceback

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

FAMILIES = ['umi','max','ratio','gauss','poisson_gauss','2beta','3beta',
            'poisson','negative_binomial','binomial','quantiles']
LIBS = ['C1','C2','T1','T2']


def write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')


def require(path):
    if not path.is_file(): raise RuntimeError(f'Incomplete dependency: {path}')


def canonical_id(spacer, position, target):
    return 'g_'+hashlib.sha256(f'{spacer}:{position}:{target}'.encode()).hexdigest()[:16]


def counts_matrix(frame, cells, guides):
    ci = pd.Index(cells).get_indexer(frame.bc_wells)
    gi = pd.Index(guides).get_indexer(frame.grna_name)
    if (gi<0).any(): raise ValueError('Unrecognized guide identity')
    keep = ci>=0
    values = frame.umi_count.to_numpy()
    if (values<0).any() or not np.equal(values,np.floor(values)).all():
        raise ValueError('Guide counts must be nonnegative integers')
    return sparse.coo_matrix((values[keep],(ci[keep],gi[keep])),
        shape=(len(cells),len(guides))).tocsr().astype(np.int64)


def identifiable_counts(frame, guides):
    """Separate explicit unpositioned design units from named guide counts.

    A spacer compatible with both arms can retain a unit_* identity when its
    position is unresolved. Keep that evidence in an audit ledger instead of
    distributing its UMI to two canonical guides. Unexpected names still fail.
    """
    known = frame.grna_name.isin(guides)
    unknown = frame.loc[~known]
    if len(unknown):
        required = {'guide_unit_id', 'position'}
        if not required <= set(unknown.columns):
            raise ValueError('Unrecognized guide identity without unit provenance')
        allowed = (unknown.grna_name.eq(unknown.guide_unit_id) &
                   unknown.grna_name.str.startswith('unit_') &
                   unknown.position.astype(int).eq(0))
        if not allowed.all():
            raise ValueError('Unrecognized positioned guide identity')
    return frame.loc[known].copy(), unknown.copy()


def save_calls(core, name, calls, meta):
    dest = core/'assignments'/name
    dest.mkdir(parents=True,exist_ok=False)
    sparse.save_npz(dest/'calls.npz',sparse.csr_matrix(calls,dtype=np.uint8))
    from scipy.io import mmwrite
    mmwrite(dest/'calls.mtx',sparse.csr_matrix(calls,dtype=np.int32).T)
    write_json(dest/'method.json',dict(method=name,**meta))
    (dest/'_SUCCESS').touch()


def pair_calls(frame, cells, guides, ref, impute=False):
    result = np.zeros((len(cells),len(guides)),bool)
    ci={v:i for i,v in enumerate(cells)}; gi={v:i for i,v in enumerate(guides)}
    mapping = {(r.pair_id,int(r.position)):r.grna_name for r in ref.itertuples()}
    for r in frame.itertuples():
        if r.bc_wells not in ci: continue
        for pos in [1,2]:
            if impute or getattr(r,f'pos{pos}_umi',0)>0:
                g=mapping.get((r.pair_id,pos))
                if g is not None: result[ci[r.bc_wells],gi[g]]=True
    return result


def export(core, include_historical=True):
    require(core/'_COUNTS_SUCCESS')
    cells=pd.read_csv(core/'cells.csv.gz',keep_default_na=False)
    ref=pd.read_csv(core/'reference.csv',keep_default_na=False)
    guides=ref[['grna_name','grna_seq','position','target_gene']].drop_duplicates().sort_values('grna_name')
    aliases=ref.groupby('grna_name').original_grna_name.agg(lambda x:';'.join(sorted(set(x))))
    guides['aliases']=guides.grna_name.map(aliases)
    dest=core/'matrices'; dest.mkdir(exist_ok=True)
    guides.to_csv(dest/'guides.csv',index=False)
    ids=cells.cell_id.tolist(); names=guides.grna_name.tolist()
    paths={'standard':core/'standard/grna_assignment/grna_umi_counts.csv.gz',
           'dedicated':core/'recovered/grna_assignment/grna_umi_counts.csv.gz',
           'combined':core/'recovered/grna_assignment/grna_umi_counts_by_source.csv.gz'}
    summaries=[]
    for arm,path in paths.items():
        frame=pd.read_csv(path,keep_default_na=False)
        frame,unpositioned=identifiable_counts(frame,names)
        unpositioned.to_csv(dest/f'unpositioned_units_{arm}.csv.gz',index=False)
        x=counts_matrix(frame,ids,names)
        sparse.save_npz(dest/f'{arm}.npz',x)
        cohort_unpositioned=unpositioned[unpositioned.bc_wells.isin(ids)]
        summaries.append(dict(arm=arm,umis=int(x.sum()),cells_with_umis=int((x.sum(axis=1)>0).sum()),
            unpositioned_umis=int(cohort_unpositioned.umi_count.sum()),
            unpositioned_cells=int(cohort_unpositioned.bc_wells.nunique())))
        if arm=='combined':
            frame[frame.bc_wells.isin(ids)].to_csv(dest/'pair_units.csv.gz',index=False)
        apath=path.parent/('grna_cell_assignment_combined.csv' if arm=='combined' else 'grna_cell_assignment.csv')
        pairs=pd.read_csv(apath,keep_default_na=False)
        if arm=='combined': pairs[pairs.bc_wells.isin(ids) & pairs.pair_id.ne('')].to_csv(dest/'pair5.csv.gz',index=False)
        save_calls(core,f'{arm}_pair5_observed',pair_calls(pairs,ids,names,ref),
                   dict(family='legacy_pair_sum',recovery=arm,impute=False))
        if arm=='combined':
            save_calls(core,'combined_pair5_imputed',pair_calls(pairs,ids,names,ref,True),
                       dict(family='legacy_pair_sum',recovery=arm,impute=True))
    # Original historical pair calls are an explicit comparator, not the new baseline.
    original=Path(__file__).resolve().parents[1]/'output/grna_assignment/grna_cell_assignment.csv'
    if include_historical and original.is_file():
        save_calls(core,'historical_pair5_observed',
                   pair_calls(pd.read_csv(original,keep_default_na=False),ids,names,ref),
                   dict(family='historical',recovery='historical',impute=False))
    x=sparse.load_npz(dest/'combined.npz')
    a=ad.AnnData(x,obs=cells.set_index('cell_id'),var=guides.set_index('grna_name'))
    a.obs['batch']=a.obs.sublibrary.map({v:i+1 for i,v in enumerate(LIBS)}).astype(int)
    a.write_h5ad(dest/'combined.h5ad')
    # Rescue only an actually observed guide in the same resolved construct.
    units=pd.read_csv(dest/'pair_units.csv.gz',keep_default_na=False)
    if 'structure_status' in units: units=units[units.structure_status.eq('concordant')]
    if 'ibar_status' in units: units=units[units.ibar_status.eq('present_resolved')]
    u=units[units.pair_id.ne('')].groupby(['bc_wells','pair_id','position','grna_name'],as_index=False).umi_count.sum()
    strong=u[u.umi_count>=5][['bc_wells','pair_id','position']].copy()
    strong['position']=3-strong.position.astype(int)
    ci={v:i for i,v in enumerate(ids)}; gi={v:i for i,v in enumerate(names)}
    for weak in [1,2]:
        selected=u[u.umi_count>=weak].merge(strong,on=['bc_wells','pair_id','position'])
        calls=(x.toarray()>=5)
        for r in selected.itertuples(): calls[ci[r.bc_wells],gi[r.grna_name]]=True
        save_calls(core,f'ibar_observed_rescue{weak}',calls,
                   dict(family='CSM_pair_corroboration',strong_umi=5,weak_umi=weak,impute=False))
    pd.DataFrame(summaries).to_csv(dest/'recovery_summary.csv',index=False)
    write_json(dest/'manifest.json',dict(cells=len(cells),guides=len(guides),
        zero_rows_retained=True,guide_identity='designed spacer + position + target',
        counts_before_assignment_threshold=True))
    (dest/'_SUCCESS').touch()


def read_crispat_assignments(dest, cells, guides):
    result=np.zeros((len(cells),len(guides)),bool)
    ci={c:i for i,c in enumerate(cells)}; gi={g:i for i,g in enumerate(guides)}
    files=sorted(dest.rglob('assignments*.csv'))
    if not files: raise RuntimeError(f'crispat produced no assignment files: {dest}')
    for file in files:
        try: frame=pd.read_csv(file)
        except pd.errors.EmptyDataError: continue
        if not len(frame): continue
        if not {'cell','gRNA'}<=set(frame): raise ValueError(f'Bad output schema: {file}')
        for cell,g in zip(frame.cell.astype(str),frame.gRNA.astype(str)):
            if cell not in ci or g not in gi: raise ValueError('Unknown crispat ID')
            result[ci[cell],gi[g]]=True
    return result


def methods(core, task):
    import crispat
    import torch
    torch.set_num_threads(1)
    require(core/'matrices/_SUCCESS')
    family=FAMILIES[task-1]
    a=ad.read_h5ad(core/'matrices/combined.h5ad')
    cells=a.obs_names.tolist(); guides=a.var_names.tolist()
    raw=core/'fits'/family; raw.mkdir(parents=True,exist_ok=False)
    variants=[('native',{})]
    if family=='umi': variants=[(f't{t}',{'thresholds':[t]}) for t in [1,2,3,5,8,10,15,20]]
    if family=='ratio': variants=[(f't{t}',{'thresholds':[t]}) for t in [0.1,0.2,0.3,0.5,0.7,0.9]]
    if family=='quantiles': variants=[(f'q{q}',{'thresholds':[q]}) for q in [0.5,0.7,0.9]]
    if family=='gauss': variants=[('all',{'nonzero':False}),('nonzero',{'nonzero':True})]
    modes=['global','position'] if family in ['max','ratio','2beta','3beta'] else ['global']
    failures=[]; completed=[]
    for mode in modes:
        for suffix,options in variants:
            name=f'crispat_{family}_{mode}_{suffix}'
            t0=time.monotonic(); combined=np.zeros(a.shape,bool)
            try:
                # Count models support an explicit batch term; other models fit per library.
                libraries=[None] if family in ['poisson','negative_binomial','binomial'] else LIBS
                for lib in libraries:
                    for pos in ([1,2] if mode=='position' else [None]):
                        subset=a[(a.obs.sublibrary==lib).to_numpy() if lib else slice(None),
                                 (a.var.position.astype(int)==pos).to_numpy() if pos else slice(None)].copy()
                        subset=subset[:,np.asarray(subset.X.sum(axis=0)).ravel()>0].copy()
                        if not subset.n_vars or not subset.n_obs: continue
                        if lib: subset.obs['batch']=1
                        folder=raw/name/f'{lib or "all"}_{pos or "all"}'; folder.mkdir(parents=True)
                        inp=folder/'input.h5ad'; subset.write_h5ad(inp)
                        kwargs=dict(options)
                        if family not in ['umi']: kwargs['UMI_threshold']=1
                        if family in ['poisson','negative_binomial','binomial']:
                            kwargs.update(parallelize=False)
                        if family in ['2beta','3beta','poisson','negative_binomial']:
                            (folder/'fitted_model_plots').mkdir(exist_ok=True)
                        if family in ['poisson','negative_binomial']:
                            (folder/'loss_plots').mkdir(exist_ok=True)
                        if family in ['gauss','poisson_gauss']: kwargs.update(n_jobs=1,make_plots=False)
                        function=getattr(crispat,'ga_'+family)
                        function(str(inp),output_dir=str(folder)+'/',**kwargs)
                        combined |= read_crispat_assignments(folder,cells,guides)
                if (combined & (a.X.toarray()==0)).any():
                    raise ValueError('Unobserved guides assigned by observed-only method')
                save_calls(core,name,combined,dict(family=family,mode=mode,options=options,
                    impute=False,seconds=time.monotonic()-t0,crispat_commit='a9372b9cf9b82b7409cc3be4d192056ae57a2673',
                    fit_scope='UMI only; transductive; no expression outcomes',
                    default_eligible=family not in ['max','ratio','2beta','3beta'],
                    selection_note='Dominant-guide families are diagnostic comparators; per-position variants still limit construct multiplicity.'))
                completed.append(name)
            except Exception:
                failure=dict(method=name,error=traceback.format_exc(),seconds=time.monotonic()-t0)
                failures.append(failure); print(failure['error'],flush=True)
    write_json(raw/'status.json',dict(family=family,completed=completed,failures=failures))
    # Execution finished, including explicitly reported failed fits; no fabricated calls.
    (raw/'_FINISHED').touch()
    if not completed: raise RuntimeError(f'All {family} variants failed')


def agreement(a,b):
    intersection=int((a & b).sum()); union=int((a | b).sum())
    return dict(intersection=intersection,union=union,
                positive_jaccard=intersection/union if union else float('nan'),
                a_only=int((a & ~b).sum()),b_only=int((b & ~a).sum()),
                exact_cell_set=float((a==b).all(axis=1).mean()),
                n_calls_a=int(a.sum()),n_calls_b=int(b.sum()))


def concordance(core):
    for family in FAMILIES: require(core/'fits'/family/'_FINISHED')
    guides=pd.read_csv(core/'matrices/guides.csv'); cells=pd.read_csv(core/'cells.csv.gz')
    dirs=sorted(p.parent for p in (core/'assignments').glob('*/_SUCCESS'))
    matrices={p.name:sparse.load_npz(p/'calls.npz').toarray().astype(bool) for p in dirs}
    rows=[]; byguide=[]; summary=[]; pair_summary=[]; target_rows=[]
    ref=pd.read_csv(core/'reference.csv',keep_default_na=False)
    units=pd.read_csv(core/'matrices/pair_units.csv.gz',keep_default_na=False)
    units=units[units.pair_id.ne('')]
    pairs=units.groupby(['bc_wells','pair_id','position']).umi_count.sum().unstack('position',fill_value=0).reset_index()
    ci=pd.Index(cells.cell_id).get_indexer(pairs.bc_wells)
    lookup={(r.pair_id,int(r.position)):r.grna_name for r in ref.itertuples()}
    pi=[pd.Index(guides.grna_name).get_indexer([lookup[(p,pos)] for p in pairs.pair_id]) for pos in [1,2]]
    target_names=sorted(set(guides.target_gene)-{'ctrl'})
    target_calls={name:np.column_stack([x[:,guides.target_gene.eq(t).to_numpy()].any(axis=1) for t in target_names])
                  for name,x in matrices.items()}
    for name,x in matrices.items():
        for lib in ['all',*LIBS]:
            keep=np.ones(len(cells),bool) if lib=='all' else cells.sublibrary.eq(lib).to_numpy()
            n=x[keep].sum(axis=1)
            summary.append(dict(method=name,library=lib,cells=int(keep.sum()),unassigned=int((n==0).sum()),
                one_guide=int((n==1).sum()),two_guides=int((n==2).sum()),more_than_two=int((n>2).sum()),
                positive_calls=int(n.sum())))
        called=np.column_stack([x[ci,pi[pos-1]] for pos in [1,2]])
        observed=called & np.column_stack([pairs.get(pos,pd.Series(0,index=pairs.index)).to_numpy()>0 for pos in [1,2]])
        for lib in ['all',*LIBS]:
            keep=np.ones(len(pairs),bool) if lib=='all' else cells.iloc[ci].sublibrary.eq(lib).to_numpy()
            count=observed[keep].sum(axis=1)
            pair_summary.append(dict(method=name,library=lib,resolved_pair_candidates=int(keep.sum()),
                neither_observed_called=int((count==0).sum()),one_observed_called=int((count==1).sum()),
                both_observed_called=int((count==2).sum()),
                called_without_same_pair_observation=int((called[keep] & ~observed[keep]).sum())))
    names=list(matrices)
    for i,name in enumerate(names):
        for other in names[i:]:
            for lib in ['all',*LIBS]:
                keep=np.ones(len(cells),bool) if lib=='all' else cells.sublibrary.eq(lib).to_numpy()
                rows.append(dict(method_a=name,method_b=other,library=lib,
                                 **agreement(matrices[name][keep],matrices[other][keep])))
            for j,g in enumerate(guides.grna_name):
                a=matrices[name][:,j]; b=matrices[other][:,j]; union=int((a|b).sum())
                byguide.append(dict(method_a=name,method_b=other,guide=g,union=union,
                    jaccard=float((a&b).sum()/union) if union else np.nan))
            target_rows.append(dict(method_a=name,method_b=other,**agreement(target_calls[name],target_calls[other])))
    dest=core/'concordance'; dest.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(dest/'method_concordance.csv',index=False)
    pd.DataFrame(byguide).to_csv(dest/'guide_concordance.csv.gz',index=False)
    pd.DataFrame(summary).to_csv(dest/'assignment_yield.csv',index=False)
    pd.DataFrame(pair_summary).to_csv(dest/'pair_completeness.csv',index=False)
    pd.DataFrame(target_rows).to_csv(dest/'target_concordance.csv',index=False)
    failed=[f for family in FAMILIES for f in json.loads((core/'fits'/family/'status.json').read_text())['failures']]
    write_json(dest/'manifest.json',dict(methods=names,failed_variants=failed,accuracy_ground_truth=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    frame=pd.DataFrame(rows); table=frame[frame.library=='all'].pivot(index='method_a',columns='method_b',values='positive_jaccard')
    table=table.combine_first(table.T).reindex(index=names,columns=names)
    fig,ax=plt.subplots(figsize=(15,13)); im=ax.imshow(table.to_numpy(),vmin=0,vmax=1,cmap='viridis')
    ax.set_xticks(range(len(names)),names,rotation=90,fontsize=6); ax.set_yticks(range(len(names)),names,fontsize=6)
    fig.colorbar(im,ax=ax,label='Positive cell–guide Jaccard'); fig.tight_layout()
    fig.savefig(dest/'concordance.png',dpi=180); plt.close(fig)
    (dest/'_SUCCESS').touch()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('stage',choices=['export','methods','concordance'])
    p.add_argument('--core',type=Path,required=True); p.add_argument('--task',type=int,default=0)
    args=p.parse_args()
    if args.stage=='methods': methods(args.core,args.task)
    else: globals()[args.stage](args.core)
