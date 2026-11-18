"""Translate the Step 4 evidence and calls into a versioned AnnData interface."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


SCHEMA = 1


def digest(path):
    with Path(path).open('rb') as stream:
        result = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024**2), b''):
            result.update(chunk)
        return result.hexdigest()


def guide_table(names, evidence, reference):
    """One record per count-matrix column; ambiguous names stay visible but cannot be rescued."""
    reference = reference.fillna('')
    ref_groups = {str(name): group for name, group in reference.groupby('grna_name')}

    def values(table, column):
        return {str(value) for value in table[column] if str(value)} if column in table else set()

    def unique_by_guide(column):
        if column not in evidence:
            return {}
        compact = evidence[['grna_name',column]].drop_duplicates().fillna('')
        return {name: values(group,column) for name,group in compact.groupby('grna_name')}

    observed_targets = unique_by_guide('target_gene')
    observed_positions = unique_by_guide('position')

    rows = []
    for name in names:
        ref = ref_groups.get(name, pd.DataFrame())
        sequences = values(ref, 'grna_seq')
        targets = values(ref, 'target_gene') | observed_targets.get(name,set())
        positions = {int(float(value.removeprefix('pos'))) for value in
                     (values(ref, 'position') or observed_positions.get(name,set()))}
        positions.discard(0)
        barcodes = sorted(values(ref, 'ibar_seq'))
        status = 'not_in_reference' if ref.empty else 'ambiguous' if (
            len(sequences) != 1 or len(targets) != 1 or len(positions) > 1) else 'resolved'
        barcode_design = ('unknown' if 'ibar_seq' not in ref else
                          'none' if not barcodes else 'mixed' if ref.ibar_seq.eq('').any() else 'present')
        rows.append(dict(guide_id=name, target_gene=next(iter(targets)) if len(targets) == 1 else '',
                         spacer=next(iter(sequences)) if len(sequences) == 1 else '',
                         position=next(iter(positions)) if len(positions) == 1 else 0,
                         identity_status=status, ibar_design=barcode_design))
    result = pd.DataFrame(rows, columns=['guide_id','target_gene','spacer','position','identity_status','ibar_design'])
    result.index = pd.Index(names, dtype=str, name='guide_key')
    return result


def assignment_matrix(rows, counts, cells, guides, dual):
    """Use accepted Step 4 rows, never a second threshold on aggregate counts."""
    if rows.empty:
        return sparse.csr_matrix(counts.shape, dtype=np.uint8)
    rows = rows.loc[rows.bc_wells.isin(cells)].fillna('')
    if dual:
        rows = rows.loc[rows.pair_id.ne('')]
        if 'detection_status' in rows:
            rows = rows.loc[rows.detection_status.isin(['complete','inferred'])]
        hits = pd.concat([rows.loc[pd.to_numeric(rows[f'pos{arm}_umi']).gt(0),
                                   ['bc_wells',f'pos{arm}_grna']].rename(columns={f'pos{arm}_grna':'grna_name'})
                          for arm in (1, 2)], ignore_index=True)
    else:
        hits = rows.loc[pd.to_numeric(rows.umi_count).gt(0), ['bc_wells','grna_name']]
    hits = hits.drop_duplicates()
    rr = pd.Index(cells).get_indexer(hits.bc_wells)
    cc = pd.Index(guides).get_indexer(hits.grna_name)
    if len(hits) and ((cc < 0).any() or (np.asarray(counts[rr, cc]).ravel() <= 0).any()):
        raise ValueError('Step 4 call is missing from the cached UMI matrix')
    return sparse.csr_matrix((np.ones(len(hits), dtype=np.uint8), (rr, cc)), shape=counts.shape)


def describe_assignments(adata, calls, features, tier='primary'):
    calls = calls.tocsr()
    targets = features.target_gene.to_numpy()
    ids = features.guide_id.to_numpy()
    indices = (calls.indices[start:stop] for start,stop in zip(calls.indptr[:-1],calls.indptr[1:]))
    values = [(json.dumps(sorted(ids[cols])), json.dumps(sorted(set(targets[cols]) - {''}))) for cols in indices]
    adata.obs[f'grna_{tier}_n_guides'] = np.diff(calls.indptr)
    adata.obs[f'grna_{tier}_guides'] = [entry[0] for entry in values]
    adata.obs[f'grna_{tier}_targets'] = [entry[1] for entry in values]


def evidence_totals(evidence, cells, names):
    retained = evidence.loc[evidence.bc_wells.isin(cells)]
    totals = retained.groupby('bc_wells').umi_count.sum().reindex(cells,fill_value=0).to_numpy(dtype=np.int64)
    unknown = retained.loc[~retained.grna_name.isin(names)]
    unresolved = unknown.groupby('bc_wells').umi_count.sum().reindex(cells,fill_value=0).to_numpy(dtype=np.int64)
    return totals, unresolved


def attach_contract(adata, counts, calls, evidence, assignment, reference, cache, options, source, source_keys,
                    cache_relative, cache_totals=None):
    names = list(adata.uns['gRNA_names'])
    features = guide_table(names, evidence, reference)
    if 'pair_id' in reference:
        constructs = reference.loc[reference.pair_id.ne('')].groupby('pair_id').target_gene.agg(
            lambda values: sorted(set(values)-{''}))
        adata.uns['grna_construct_targets'] = {str(key): value for key,value in constructs.items()}
    total, omitted = (evidence_totals(evidence, adata.obs_names, names) if cache_totals is None else cache_totals)
    if not np.array_equal(np.asarray(counts.sum(axis=1)).ravel(), total - omitted):
        raise ValueError('Step 4 guide matrix and UMI ledger totals do not match')
    adata.uns['grna_features'] = features
    adata.obs['grna_observed_umis'] = total
    adata.obs['grna_observed_n_guides'] = counts.getnnz(axis=1)
    adata.obs['grna_unindexed_umis'] = omitted
    ambiguous = features.identity_status.ne('resolved').to_numpy()
    adata.obs['grna_unresolved_identity_umis'] = np.asarray(counts[:,ambiguous].sum(axis=1)).ravel() + omitted
    describe_assignments(adata, calls, features)
    if options.get('dual_guide_counting_mode',False):
        if 'pair_id' in assignment:
            compatible = assignment.loc[assignment.pair_id.fillna('').ne('')].copy()
        else:
            compatible = assignment.iloc[:0].copy()
        linked = compatible.groupby('bc_wells').pair_id.nunique() if 'pair_id' in compatible else pd.Series(dtype=int)
        adata.obs['grna_primary_n_resolved_constructs'] = linked.reindex(adata.obs_names,fill_value=0).astype(int)
        partners = []
        for arm in (1, 2):
            col = f'inferred_pos{arm}_grna'
            if col in compatible:
                kept = compatible.loc[compatible[col].fillna('').ne(''), ['bc_wells','pair_id',col]].copy()
                kept.columns = ['bc_wells','construct_id','guide_id']
                kept['position'] = arm
                partners.append(kept)
        inferred = pd.concat(partners,ignore_index=True).drop_duplicates() if partners else pd.DataFrame()
        if not inferred.empty:
            inferred = inferred.loc[inferred.bc_wells.isin(adata.obs_names)].reset_index(drop=True)
            inferred.index = inferred.index.astype(str)
            adata.uns['grna_inferred_partners'] = inferred.rename(columns={'bc_wells':'cell_id'})
    adata.uns['grna_schema'] = dict(
        version=SCHEMA, source=source, count_key='gRNA_counts', primary_key='gRNA_assignment',
        feature_key='grna_features', source_suffixes=source_keys,
        assignment_mode='pair_sum' if options.get('dual_guide_counting_mode',False) else 'single',
        primary_umi_threshold=int(options['min_umi_threshold']),
        multiplet_strategy=options.get('multiplet_strategy','all'),
        cache_path=cache_relative, cache_sha256=digest(cache),
    )


def summarize_cache(path, cells, names):
    """Stream the ledger instead of keeping reads for all non-GEX barcodes in memory."""
    keys = {'bc_wells','grna_name','target_gene','position','umi_count'}
    indexed = set(names)
    cell_index = pd.Index(cells)
    total = np.zeros(len(cells),dtype=np.int64)
    omitted = total.copy()
    metadata = None
    for frame in pd.read_csv(path,usecols=lambda col:col in keys,keep_default_na=False,
                             dtype={'bc_wells':str,'grna_name':str},chunksize=300_000):
        rows = cell_index.get_indexer(frame.bc_wells)
        use = rows>=0
        valid = frame.loc[use]
        unknown = ~valid.grna_name.isin(indexed).to_numpy()
        np.add.at(total,rows[use],valid.umi_count.to_numpy(dtype=np.int64))
        np.add.at(omitted,rows[use][unknown],valid.umi_count.to_numpy(dtype=np.int64)[unknown])
        observed = valid.drop(columns=['bc_wells','umi_count']).drop_duplicates()
        metadata = observed if metadata is None else pd.concat([metadata,observed]).drop_duplicates()
    return metadata if metadata is not None else pd.DataFrame(columns=['grna_name']), (total, omitted)
