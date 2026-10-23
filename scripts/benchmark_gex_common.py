"""Matrix loading and agreement helpers for the v2 Step 1a benchmark."""

import numpy as np
import pandas as pd
from scipy import io as sio, sparse, stats


def correlation(first, second, method="pearson"):
    first, second = np.asarray(first), np.asarray(second)
    finite = np.isfinite(first) & np.isfinite(second)
    first, second = first[finite], second[finite]
    if len(first) < 3 or np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(stats.spearmanr(first, second).statistic if method == "spearman"
                 else np.corrcoef(first, second)[0, 1])


def sparse_cell_correlations(first, second):
    """Per-cell Pearson across all genes, including shared implicit zeros."""
    first, second = first.copy().astype(float).tocsr(), second.copy().astype(float).tocsr()
    first.data = np.log1p(first.data)
    second.data = np.log1p(second.data)
    n = first.shape[1]
    sum1, sum2 = np.asarray(first.sum(axis=1)).ravel(), np.asarray(second.sum(axis=1)).ravel()
    cross = np.asarray(first.multiply(second).sum(axis=1)).ravel() - sum1 * sum2 / n
    var1 = np.asarray(first.multiply(first).sum(axis=1)).ravel() - sum1**2 / n
    var2 = np.asarray(second.multiply(second).sum(axis=1)).ravel() - sum2**2 / n
    denominator = np.sqrt(np.maximum(var1 * var2, 0))
    return np.divide(cross, denominator, out=np.full(len(cross), np.nan), where=denominator > 0)


def load_product(path, process=False):
    if process:
        matrix_name, obs_name = "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv"
    else:
        matrix_name, obs_name = "count_matrix.mtx", "cell_metadata.csv"
    matrix = sparse.csr_matrix(sio.mmread(path / matrix_name))
    obs = pd.read_csv(path / obs_name, index_col="bc_wells")
    genes = pd.read_csv(path / "all_genes.csv", dtype=str)
    obs.index = obs.index.astype(str)
    if matrix.shape != (len(obs), len(genes)) or not obs.index.is_unique or genes["gene_id"].duplicated().any():
        raise ValueError(f"Malformed benchmark product: {path}")
    return matrix, obs, genes


def markdown_table(frame, columns):
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        lines.append("| " + " | ".join(f"{value:.3g}" if isinstance(value, (float, np.floating)) else str(value)
                                       for value in row) + " |")
    return "\n".join(lines)
