"""Preserve prior guide evidence and warn before rebuilding integrated AnnData."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4


def run_id():
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid4().hex[:8]


def archive_step1c(namespace, library):
    root = Path(namespace)
    names = [root/'manifests'/f'{library}.step1c.complete.json']
    names += [root/'grna_matches'/f'{library}{suffix}' for suffix in (
        '_grna_matches.tsv.gz', '_grna_matches_unfiltered.tsv.gz', '_matching_stats.yaml',
        '_recovery.json', '_design.yaml', '_relaxed.complete.json')]
    existing = [path for path in names if path.exists()]
    if not existing:
        return None
    archive = root/'restarts'/f'{run_id()}_{library}'
    archive.mkdir(parents=True)
    print(f'WARNING: rebuilding Step 1c for {library}; previous guide evidence is preserved in {archive}', flush=True)
    moves = []
    for path in existing:
        dest = archive/path.relative_to(root)
        dest.parent.mkdir(parents=True, exist_ok=True)
        path.rename(dest)
        moves.append({'original': str(path), 'archived': str(dest)})
    (archive/'relocation.json').write_text(json.dumps(moves, indent=2)+'\n')
    return archive


def write_anndata(adata, destination):
    """Write successfully before moving the old file; retain all old annotations."""
    import h5py
    destination = Path(destination)
    identifier = run_id()
    archive = destination.parent/'grna_restarts'/identifier
    summary = {}
    if destination.exists():
        with h5py.File(destination, 'r') as handle:
            for slot in ('obs', 'var'):
                columns = list(handle[slot].attrs.get('column-order', [])) if slot in handle else []
                summary[slot] = [str(x.decode() if isinstance(x, bytes) else x) for x in columns]
            for slot in ('obsm', 'varm', 'layers', 'uns', 'obsp', 'varp'):
                summary[slot] = list(handle[slot]) if slot in handle else []
        print('WARNING: Step 5 rebuilds AnnData from GEX and guide products. Existing downstream '
              'annotations, embeddings, layers and analyses are not merged into the rebuilt file.', flush=True)
        print(f'The complete previous AnnData will be preserved at {archive/destination.name}', flush=True)
        print('Previous annotation inventory: ' + json.dumps(summary), flush=True)
    temporary = destination.with_name(f'.{destination.stem}.{identifier}.h5ad')
    adata.write_h5ad(temporary)
    if destination.exists():
        archive.mkdir(parents=True)
        destination.rename(archive/destination.name)
        (archive/'annotations.json').write_text(json.dumps(summary, indent=2)+'\n')
    os.replace(temporary, destination)
