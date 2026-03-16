#!/usr/bin/env python3
"""
QC Q2: Recombination & Background Noise Estimation

Part A — Ambient Noise Estimation & Guide UMI Distributions:
  - estimate_background_noise: computes per-cell off-target UMI stats.
  - compute_guide_umi_distributions: estimates mean/std of log10(UMI) per guide.

Part B — Recombination Model (MLE) & Decoupling QC:
  - Generative model estimating dropout (d), silencing (s), and recombination (r)
    rates from observed cell category counts.
  - analyze_grna_ibar_pairing: checks gRNA-iBAR decoupling events.

Part C — UMI Threshold Sweep:
  - Re-runs pair assignment across UMI threshold values, estimates recombination
    parameters for each threshold, and plots sensitivity.

Usage:
  python scripts/qc_q2_recombination.py --config config.yaml
"""

import argparse
import os
import sys
import warnings
from collections import defaultdict
from math import factorial, exp, log
import importlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.optimize import minimize
from tqdm import tqdm

_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from scripts.qc_config_loader import load_qc_config

# Import assignment function for Part C threshold sweep
assign_mod = importlib.import_module("scripts.04_assign_grna_to_cells")
assign_dual_guide_pairs = assign_mod.assign_dual_guide_pairs


# ============================================================================
# Part A: Ambient Noise Estimation & Guide UMI Distributions
# ============================================================================

def estimate_background_noise(adata, assignment_df, umi_df, metric="total_umi", cell_category="perfect_singlet"):
    """Estimate background noise floor (off-target UMI counts).
    
    Args:
        adata: AnnData object
        assignment_df: Cell assignment DataFrame
        umi_df: UMI counts DataFrame
        metric: "total_umi", "n_guides", or "mean_umi".
        cell_category: "perfect_singlet" or "unassigned".
        
    Returns:
        pd.Series of per-cell noise values.
    """
    print(f"\n  --- Part A.1: Background Noise Estimation ({cell_category}, {metric}) ---")
    
    valid_bcs = set(adata.obs.index)
    ass_df = assignment_df[assignment_df['bc_wells'].isin(valid_bcs)]
    
    if cell_category == "perfect_singlet":
        target_cells = ass_df[(ass_df['detection_status'] == 'complete') & 
                              (ass_df['combo_type'] == 'singlet')]
        if target_cells.empty:
            print("    WARNING: No perfect singlets found for noise estimation.")
            return pd.Series(dtype=float)
            
        print(f"    Found {len(target_cells):,} perfect singlets for noise ")
        
        target_bcs = set(target_cells['bc_wells'])
        sub_umi = umi_df[umi_df['bc_wells'].isin(target_bcs)]
        
        merged = sub_umi.merge(
            target_cells[['bc_wells', 'pos1_grna', 'pos2_grna']], 
            on='bc_wells', 
            how='left'
        )
        is_off_target = (merged['grna_name'] != merged['pos1_grna']) & (merged['grna_name'] != merged['pos2_grna'])
        off_target_reads = merged[is_off_target]
        
    elif cell_category == "unassigned":
        assigned_bcs = set(ass_df['bc_wells'])
        unassigned_bcs = valid_bcs - assigned_bcs
        
        if not unassigned_bcs:
            print("    WARNING: No unassigned valid GEX cells found for noise estimation.")
            return pd.Series(dtype=float)
            
        print(f"    Found {len(unassigned_bcs):,} unassigned cells for noise estimation.")
        
        target_bcs = unassigned_bcs
        off_target_reads = umi_df[umi_df['bc_wells'].isin(target_bcs)]
        
    else:
        raise ValueError(f"Unknown cell_category '{cell_category}'")
    
    if metric == "total_umi":
        per_cell_noise = off_target_reads.groupby('bc_wells')['umi_count'].sum()
    elif metric == "n_guides":
        per_cell_noise = off_target_reads.groupby('bc_wells')['grna_name'].nunique()
    elif metric == "mean_umi":
        per_cell_noise = off_target_reads.groupby('bc_wells')['umi_count'].mean()
        per_cell_noise = per_cell_noise.fillna(0)
    else:
        print(f"    WARNING: Unknown metric '{metric}', defaulting to 'total_umi'")
        per_cell_noise = off_target_reads.groupby('bc_wells')['umi_count'].sum()
        
    per_cell_noise = per_cell_noise.reindex(list(target_bcs), fill_value=0)
    
    mean_val = per_cell_noise.mean()
    median_val = per_cell_noise.median()
    p95_val = np.percentile(per_cell_noise, 95)
    
    print(f"    Global Off-Target Noise Distribution (per cell):")
    print(f"      Mean:            {mean_val:.2f}")
    print(f"      Median:          {median_val:.2f}")
    print(f"      95th Percentile: {p95_val:.2f}")
    
    return per_cell_noise


def compute_guide_umi_distributions(assignment_df, umi_df, valid_bcs=None):
    """Compute log10(UMI) distribution stats per guide for GMM refinement.

    Analyzes perfect singlet cells to find the typical UMI count 
    distribution for each successfully assigned guide unit.

    Args:
        assignment_df: Cell assignment DataFrame
        umi_df: UMI counts DataFrame
        valid_bcs: Optional set of valid GEX cell barcodes

    Returns:
        DataFrame: [guide_name, mean_log10_umi, std_log10_umi, n_cells]
    """
    print(f"\n  --- Part A.2: Guide UMI Distributions ---")
    
    if valid_bcs is not None:
        ass_df = assignment_df[assignment_df['bc_wells'].isin(valid_bcs)]
    else:
        ass_df = assignment_df
        
    perfect_singlets = ass_df[(ass_df['detection_status'] == 'complete') & 
                              (ass_df['combo_type'] == 'singlet')]
                              
    if perfect_singlets.empty:
        print("    WARNING: No perfect singlets found for UMI distributions.")
        return pd.DataFrame()
        
    print(f"    Using {len(perfect_singlets):,} perfect singlets.")
    
    # Extract pos1 and pos2 UMI counts
    pos1_counts = perfect_singlets[['pos1_grna', 'pos1_umi']].rename(columns={'pos1_grna': 'guide_name', 'pos1_umi': 'umi'})
    pos2_counts = perfect_singlets[['pos2_grna', 'pos2_umi']].rename(columns={'pos2_grna': 'guide_name', 'pos2_umi': 'umi'})
    
    all_counts = pd.concat([pos1_counts, pos2_counts], ignore_index=True)
    all_counts = all_counts[all_counts['umi'] > 0].copy()
    all_counts['log10_umi'] = np.log10(all_counts['umi'])
    
    stats = all_counts.groupby('guide_name').agg(
        mean_log10_umi=('log10_umi', 'mean'),
        std_log10_umi=('log10_umi', 'std'),
        n_cells=('log10_umi', 'count')
    ).reset_index()
    
    # Handle std for guides with 1 cell
    stats['std_log10_umi'] = stats['std_log10_umi'].fillna(0.0)
    
    print(f"    Computed UMI distributions for {len(stats):,} unique guides.")
    return stats


# ============================================================================
# Part B: Recombination Model (MLE) & Decoupling QC
# ============================================================================

def per_construct_probs(d, s, r_swap, r12, r_del, c=0.0):
    q = 1.0 - s           # prob active
    p = 1.0 - d           # prob guide detected
    p_int = 1.0 - 2*r_swap - r12 - 2*r_del   # prob intact
    coupled = 1.0 - c

    P_C = q * p_int * p * p
    P_I = q * (p_int * 2*p*d + 2*r_swap * coupled * p*d + r12 * coupled * 2*p*d + 2*r_del * p)
    P_M = q * (2*r_swap + r12) * coupled * p * p
    P_X = 1.0 - P_C - P_I - P_M

    return P_C, P_I, P_M, P_X

def per_construct_probs_simple(d, s, r, c=0.0):
    q = 1.0 - s
    p = 1.0 - d
    p_int = 1.0 - r
    coupled = 1.0 - c

    P_C = q * p_int * p * p
    P_I = q * (p_int * 2*p*d + r * coupled * 2*p*d)
    P_M = q * r * coupled * p * p
    P_X = 1.0 - P_C - P_I - P_M

    return P_C, P_I, P_M, P_X

def _truncated_poisson(k, lam):
    if k < 1:
        return 0.0
    return (lam**k * exp(-lam) / factorial(k)) / (1.0 - exp(-lam))

def cell_probs_single(P_C, P_I, P_M, P_X, moi, k_max=4):
    probs = {'perfect_singlet': 0.0, 'inferred_singlet': 0.0,
             'multiplet': 0.0, 'unassigned': 0.0}

    for k in range(1, k_max + 1):
        pk = _truncated_poisson(k, moi)
        if pk < 1e-15:
            continue

        p_ps = k * P_C * P_X**(k-1)
        p_is = k * P_I * P_X**(k-1)
        p_un = P_X**k

        p_ms = k * P_M * P_X**(k-1)
        P_visible = P_C + P_I + P_M
        p_0_visible = P_X**k
        p_1_visible = k * P_visible * P_X**(k-1)
        p_2plus_visible = 1.0 - p_0_visible - p_1_visible

        p_mp = p_ms + p_2plus_visible

        probs['perfect_singlet'] += pk * p_ps
        probs['inferred_singlet'] += pk * p_is
        probs['multiplet'] += pk * p_mp
        probs['unassigned'] += pk * p_un

    return probs

def doublet_convolution(p_single):
    ps = p_single['perfect_singlet']
    pi = p_single['inferred_singlet']
    pm = p_single['multiplet']
    pu = p_single['unassigned']

    return {
        'perfect_singlet': 2 * ps * pu,
        'inferred_singlet': 2 * pi * pu,
        'multiplet': 1.0 - 2*ps*pu - 2*pi*pu - pu*pu,
        'unassigned': pu * pu,
    }

def observable_probs(theta, moi, delta, extended=False, decoupling_rate=0.0, granular=False):
    if granular:
        return observable_probs_v2(theta, moi, delta)
    if extended:
        d, s, r_swap, r12, r_del = theta
        P_C, P_I, P_M, P_X = per_construct_probs(d, s, r_swap, r12, r_del, c=decoupling_rate)
    else:
        d, s, r = theta
        P_C, P_I, P_M, P_X = per_construct_probs_simple(d, s, r, c=decoupling_rate)

    p_single = cell_probs_single(P_C, P_I, P_M, P_X, moi)
    p_doublet = doublet_convolution(p_single)

    return {
        cat: (1.0 - delta) * p_single[cat] + delta * p_doublet[cat]
        for cat in p_single
    }


# --- v2 granular model (7 categories, position-specific deletion) ---

# Category keys used throughout v2
_V2_CATS = [
    'perfect_singlet',
    'inferred_singlet_pos1_missing',
    'inferred_singlet_pos2_missing',
    'multiplet_all_complete',
    'multiplet_all_inferred',
    'multiplet_mixed',
    'unassigned',
]

def per_construct_probs_v2(d, s, r_swap, r12, r_del1, r_del2):
    """Per-construct signal probabilities with position-specific deletion.

    Returns (P_C, P_I1, P_I2, P_M, P_X) — 5 construct-level signal types:
        P_C:  complete pair (1 pair_id, both positions detected)
        P_I1: inferred, pos1 missing (only pos2 detected)
        P_I2: inferred, pos2 missing (only pos1 detected)
        P_M:  multiplet signal (2 pair_ids from recombination, both half-detected)
        P_X:  invisible (silenced or all guides dropped)
    """
    q = 1.0 - s
    p = 1.0 - d
    p_int = 1.0 - 2*r_swap - r12 - r_del1 - r_del2

    P_C = q * p_int * p * p
    # P_I1: pos1 missing. Sources: intact+pos1_drop, any swap+pos1_drop, pos1_del
    P_I1 = q * p * ((p_int + 2*r_swap + r12) * d + r_del1)
    # P_I2: pos2 missing. Sources: intact+pos2_drop, any swap+pos2_drop, pos2_del
    P_I2 = q * p * ((p_int + 2*r_swap + r12) * d + r_del2)
    # P_M: both detected but from different pair_ids (swap/cross-swap)
    P_M = q * (2*r_swap + r12) * p * p
    P_X = 1.0 - P_C - P_I1 - P_I2 - P_M

    return P_C, P_I1, P_I2, P_M, P_X


def cell_probs_single_v2(P_C, P_I1, P_I2, P_M, P_X, moi, k_max=4):
    """Cell-level probabilities for 7 observable categories (no doublet).

    For k constructs, each independently produces one of {C, I1, I2, M, X}.
    Cell classification:
        - 1 C visible, rest X  → perfect_singlet
        - 1 I1 visible, rest X → inferred_singlet_pos1_missing
        - 1 I2 visible, rest X → inferred_singlet_pos2_missing
        - ≥2 C visible, no I/M → multiplet_all_complete (MOI)
        - ≥1 M or ≥2 I visible, no C → multiplet_all_inferred
        - ≥1 C and ≥1 I/M visible → multiplet_mixed
        - all X → unassigned
    """
    probs = {cat: 0.0 for cat in _V2_CATS}

    for k in range(1, k_max + 1):
        pk = _truncated_poisson(k, moi)
        if pk < 1e-15:
            continue

        Xk = P_X ** k  # all invisible
        CX_k = (P_C + P_X) ** k  # all constructs either C or X
        notC_k = (1.0 - P_C) ** k  # no construct produces C

        p_ps = k * P_C * P_X ** (k - 1)
        p_i1 = k * P_I1 * P_X ** (k - 1)
        p_i2 = k * P_I2 * P_X ** (k - 1)
        p_mc = CX_k - p_ps - Xk  # ≥2 C visible, nothing else
        p_mi = notC_k - Xk - k * (P_I1 + P_I2) * P_X ** (k - 1)  # no C, ≥2 inferred pair_ids
        p_mm = 1.0 - CX_k - notC_k + Xk  # ≥1 C and ≥1 non-C visible
        p_un = Xk

        probs['perfect_singlet'] += pk * p_ps
        probs['inferred_singlet_pos1_missing'] += pk * p_i1
        probs['inferred_singlet_pos2_missing'] += pk * p_i2
        probs['multiplet_all_complete'] += pk * p_mc
        probs['multiplet_all_inferred'] += pk * p_mi
        probs['multiplet_mixed'] += pk * p_mm
        probs['unassigned'] += pk * p_un

    return probs


def doublet_convolution_v2(p_single):
    """Doublet merge for 7 categories.

    Two independent cells merge. Classification rules:
        - singlet + unassigned → same singlet type
        - multiplet + unassigned → same multiplet type
        - 2 complete singlets → multiplet_all_complete
        - 2 inferred singlets (any) → multiplet_all_inferred
        - complete singlet + inferred singlet → multiplet_mixed
        - complete + all_inferred multiplet → multiplet_mixed
        - inferred singlet + all_complete multiplet → multiplet_mixed
        etc.
    """
    ps = p_single['perfect_singlet']
    pi1 = p_single['inferred_singlet_pos1_missing']
    pi2 = p_single['inferred_singlet_pos2_missing']
    mc = p_single['multiplet_all_complete']
    mi = p_single['multiplet_all_inferred']
    mm = p_single['multiplet_mixed']
    pu = p_single['unassigned']

    pi = pi1 + pi2  # total inferred singlet

    d_ps = 2 * ps * pu
    d_i1 = 2 * pi1 * pu
    d_i2 = 2 * pi2 * pu
    # All-complete: ps+ps, ps+mc, mc+mc, mc+pu
    d_mc = ps*ps + 2*ps*mc + mc*mc + 2*mc*pu
    # All-inferred: pi+pi combos, pi+mi, mi+mi, mi+pu
    d_mi = pi*pi + 2*pi*mi + mi*mi + 2*mi*pu
    d_un = pu * pu
    # Mixed: everything else
    d_mm = 1.0 - d_ps - d_i1 - d_i2 - d_mc - d_mi - d_un

    return {
        'perfect_singlet': d_ps,
        'inferred_singlet_pos1_missing': d_i1,
        'inferred_singlet_pos2_missing': d_i2,
        'multiplet_all_complete': d_mc,
        'multiplet_all_inferred': d_mi,
        'multiplet_mixed': d_mm,
        'unassigned': d_un,
    }


def observable_probs_v2(theta, moi, delta):
    """Observable probabilities for the 7-category granular model.

    Args:
        theta: (d, s, r_swap, r12, r_del1, r_del2)
        moi: MOI (lambda)
        delta: doublet rate
    """
    d, s, r_swap, r12, r_del1, r_del2 = theta
    P_C, P_I1, P_I2, P_M, P_X = per_construct_probs_v2(d, s, r_swap, r12, r_del1, r_del2)

    p_single = cell_probs_single_v2(P_C, P_I1, P_I2, P_M, P_X, moi)
    p_doublet = doublet_convolution_v2(p_single)

    return {
        cat: (1.0 - delta) * p_single[cat] + delta * p_doublet[cat]
        for cat in _V2_CATS
    }

def neg_log_likelihood(theta, counts, moi, delta, extended=False, decoupling_rate=0.0):
    eps = 1e-10
    if extended:
        d, s, r_swap, r12, r_del = theta
        if (d < eps or d > 1 - eps or s < eps or s > 1 - eps or
                r_swap < 0 or r12 < 0 or r_del < 0 or
                2*r_swap + r12 + 2*r_del > 1 - eps):
            return 1e20
    else:
        d, s, r = theta
        if d < eps or d > 1 - eps or s < eps or s > 1 - eps or r < 0 or r > 1 - eps:
            return 1e20

    try:
        probs = observable_probs(theta, moi, delta, extended, decoupling_rate=decoupling_rate)
    except Exception:
        return 1e20

    nll = 0.0
    for cat, n in counts.items():
        p = max(probs.get(cat, eps), eps)
        nll -= n * log(p)

    return nll

def estimate_stepwise(counts, moi, delta):
    n1 = counts['perfect_singlet']
    n2 = counts['inferred_singlet']
    n3 = counts['multiplet']
    n4 = counts['unassigned']
    N = n1 + n2 + n3 + n4

    R = n1 / max(n1 + n2, 1)
    d = max(0.01, min(0.5, (1.0 - R) / (1.0 + R)))

    p_unassigned = n4 / max(N, 1)
    d2 = d * d
    s = max(0.01, min(0.5, (p_unassigned - d2) / max(1.0 - d2, 0.01)))

    q = 1.0 - s
    p = 1.0 - d
    P_C_no_r = q * p * p
    P_I_no_r = q * 2 * p * d
    P_X_no_r = 1.0 - P_C_no_r - P_I_no_r

    p_single_no_r = cell_probs_single(P_C_no_r, P_I_no_r, 0.0, P_X_no_r, moi)
    p_doublet_no_r = doublet_convolution(p_single_no_r)
    p_multi_no_r = (1.0 - delta) * p_single_no_r['multiplet'] + delta * p_doublet_no_r['multiplet']

    n_multi_expected = N * p_multi_no_r
    n_surplus = max(0, n3 - n_multi_expected)

    p_k1 = _truncated_poisson(1, moi)
    denom = max(N * p_k1 * q * p * p, 1)
    r = max(0.001, min(0.3, n_surplus / denom))

    return d, s, r

def fit_mle(counts, moi, delta, extended=False, decoupling_rate=0.0, verbose=True):
    d0, s0, r0 = estimate_stepwise(counts, moi, delta)

    if extended:
        x0 = [d0, s0, r0/5, r0/5, r0/5]
        bounds = [(0.001, 0.99), (0.001, 0.99),
                  (0.0, 0.3), (0.0, 0.3), (0.0, 0.3)]
        param_names = ['d', 's', 'r_swap', 'r12', 'r_del']
    else:
        x0 = [d0, s0, r0]
        bounds = [(0.001, 0.99), (0.001, 0.99), (0.0, 0.5)]
        param_names = ['d', 's', 'r']

    result = minimize(
        neg_log_likelihood,
        x0,
        args=(counts, moi, delta, extended, decoupling_rate),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': 1000, 'ftol': 1e-12}
    )

    theta_hat = result.x
    probs_hat = observable_probs(theta_hat, moi, delta, extended, decoupling_rate=decoupling_rate)

    if verbose:
        N = sum(counts.values())
        print(f"\n    MLE Results (extended={extended}):")
        print(f"    {'Parameter':<12} {'Estimate':>10}")
        print(f"    {'─'*24}")
        for name, val in zip(param_names, theta_hat):
            print(f"    {name:<12} {val:>10.4f}")
        print(f"\n    NLL: {result.nll if hasattr(result, 'nll') else result.fun:.2f}")
        print(f"    Converged: {result.success}")

        print(f"\n    {'Category':<22} {'Observed':>10} {'Expected':>10} {'Residual':>10}")
        print(f"    {'─'*54}")
        for cat in counts:
            obs = counts[cat]
            exp_n = probs_hat.get(cat, 0) * N
            resid = obs - exp_n
            print(f"    {cat:<22} {obs:>10,} {exp_n:>10,.0f} {resid:>+10,.0f}")

    return {
        'params': dict(zip(param_names, theta_hat)),
        'param_names': param_names,
        'theta': theta_hat,
        'nll': result.fun,
        'success': result.success,
        'probs': probs_hat,
        'counts': counts,
    }


def neg_log_likelihood_v2(theta, counts, moi, delta):
    """Negative log-likelihood for the 7-category granular model."""
    eps = 1e-10
    d, s, r_swap, r12, r_del1, r_del2 = theta
    if (d < eps or d > 1 - eps or s < eps or s > 1 - eps or
            r_swap < 0 or r12 < 0 or r_del1 < 0 or r_del2 < 0 or
            2*r_swap + r12 + r_del1 + r_del2 > 1 - eps):
        return 1e20

    try:
        probs = observable_probs_v2(theta, moi, delta)
    except Exception:
        return 1e20

    nll = 0.0
    for cat, n in counts.items():
        p = max(probs.get(cat, eps), eps)
        nll -= n * log(p)
    return nll


def estimate_stepwise_v2(counts, moi, delta):
    """Step-wise initializer for the 7-category granular model.

    Returns (d, s, r_swap, r12, r_del1, r_del2).
    """
    n1 = counts['perfect_singlet']
    n2 = counts['inferred_singlet_pos1_missing']
    n3 = counts['inferred_singlet_pos2_missing']
    n_multi = (counts['multiplet_all_complete']
              + counts['multiplet_all_inferred']
              + counts['multiplet_mixed'])
    n_un = counts['unassigned']
    N = n1 + n2 + n3 + n_multi + n_un

    # Step 1: Dropout d from perfect vs total singlet ratio
    n_inf = n2 + n3
    R = n1 / max(n1 + n_inf, 1)
    d = max(0.01, min(0.5, (1.0 - R) / (1.0 + R)))

    # Step 2: Silencing s from unassigned fraction
    p_unassigned = n_un / max(N, 1)
    d2 = d * d
    s = max(0.01, min(0.5, (p_unassigned - d2) / max(1.0 - d2, 0.01)))

    # Step 3: Position-specific deletion from pos1_missing vs pos2_missing asymmetry
    # Under intact construct: P_I1/P_I2 ≈ (pint*d + r_del1)/(pint*d + r_del2)
    # Use fractional surplus to estimate r_del1, r_del2
    frac2 = n2 / max(n2 + n3, 1)  # fraction pos1_missing
    frac3 = n3 / max(n2 + n3, 1)  # fraction pos2_missing
    # Baseline: equal → frac2 = frac3 = 0.5. Excess gives deletion estimate.
    p = 1.0 - d
    q = 1.0 - s
    base_inf_rate = q * p * d  # baseline inferred rate per pos (no recomb/del)
    # Rough: r_del1 ≈ (n2/N - base_inf_rate) / (q*p) etc.
    r_del1_est = max(0.0, (n2 / max(N, 1) - base_inf_rate) / max(q * p, 0.01))
    r_del2_est = max(0.0, (n3 / max(N, 1) - base_inf_rate) / max(q * p, 0.01))
    r_del1 = max(0.001, min(0.15, r_del1_est))
    r_del2 = max(0.001, min(0.15, r_del2_est))

    # Step 4: Recombination from multiplet surplus
    # At r=0, the only multiplets come from k≥2 (MOI)
    P_C_no_r = q * p * p
    P_I_no_r = q * 2 * p * d
    P_X_no_r = 1.0 - P_C_no_r - P_I_no_r
    p_single_no_r = cell_probs_single(P_C_no_r, P_I_no_r, 0.0, P_X_no_r, moi)
    p_doublet_no_r = doublet_convolution(p_single_no_r)
    p_multi_no_r = ((1.0 - delta) * p_single_no_r['multiplet']
                    + delta * p_doublet_no_r['multiplet'])
    n_multi_expected = N * p_multi_no_r
    n_surplus = max(0, n_multi - n_multi_expected)
    p_k1 = _truncated_poisson(1, moi)
    denom = max(N * p_k1 * q * p * p, 1)
    r_total = max(0.001, min(0.3, n_surplus / denom))

    # Split r_total into r_swap and r12 (rough: r_swap=r_total/3, r12=r_total/3)
    r_swap = max(0.001, r_total / 3.0)
    r12 = max(0.001, r_total / 3.0)

    return d, s, r_swap, r12, r_del1, r_del2


def fit_mle_v2(counts, moi, delta, verbose=True):
    """Fit the 7-category granular model via MLE.

    Parameters: d, s, r_swap, r12, r_del1, r_del2
    """
    d0, s0, r_swap0, r12_0, r_del1_0, r_del2_0 = estimate_stepwise_v2(counts, moi, delta)

    x0 = [d0, s0, r_swap0, r12_0, r_del1_0, r_del2_0]
    bounds = [
        (0.001, 0.99),  # d
        (0.001, 0.99),  # s
        (0.0, 0.3),     # r_swap
        (0.0, 0.3),     # r12
        (0.0, 0.3),     # r_del1
        (0.0, 0.3),     # r_del2
    ]
    param_names = ['d', 's', 'r_swap', 'r12', 'r_del1', 'r_del2']

    result = minimize(
        neg_log_likelihood_v2,
        x0,
        args=(counts, moi, delta),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': 2000, 'ftol': 1e-12}
    )

    theta_hat = result.x
    probs_hat = observable_probs_v2(theta_hat, moi, delta)

    # Compute BIC for model comparison
    n_params = len(param_names)
    N = sum(counts.values())
    bic = n_params * log(max(N, 1)) + 2 * result.fun

    if verbose:
        print(f"\n    MLE Results (granular, 7-category):")
        print(f"    {'Parameter':<12} {'Estimate':>10}")
        print(f"    {'─'*24}")
        for name, val in zip(param_names, theta_hat):
            print(f"    {name:<12} {val:>10.4f}")
        print(f"\n    NLL: {result.fun:.2f}  BIC: {bic:.2f}")
        print(f"    Converged: {result.success}")

        print(f"\n    {'Category':<32} {'Observed':>10} {'Expected':>10} {'Residual':>10}")
        print(f"    {'─'*64}")
        for cat in _V2_CATS:
            obs = counts.get(cat, 0)
            exp_n = probs_hat.get(cat, 0) * N
            resid = obs - exp_n
            print(f"    {cat:<32} {obs:>10,} {exp_n:>10,.0f} {resid:>+10,.0f}")

    return {
        'params': dict(zip(param_names, theta_hat)),
        'param_names': param_names,
        'theta': theta_hat,
        'nll': result.fun,
        'bic': bic,
        'success': result.success,
        'probs': probs_hat,
        'counts': counts,
        'granular': True,
    }


def count_model_categories(assignment_df, valid_bcs=None, granular=False):
    """Count cells in each observable category for the model.

    Args:
        assignment_df: Cell assignment DataFrame from Step 4.
        valid_bcs: Set of valid GEX cell barcodes (for unassigned count).
        granular: If True, return 7-category dict; if False, return legacy 4-category.

    Returns:
        dict of category counts.
    """
    if valid_bcs is not None:
        valid_bcs = set(valid_bcs)
        df_valid = assignment_df[assignment_df['bc_wells'].isin(valid_bcs)]
        n_unassigned = len(valid_bcs) - df_valid['bc_wells'].nunique()
    else:
        df_valid = assignment_df
        n_unassigned = 0

    cell_meta = df_valid.drop_duplicates('bc_wells')

    # --- Singlets ---
    singlets = df_valid[df_valid['combo_type'] == 'singlet']
    n_perfect = int((singlets['detection_status'] == 'complete').sum())

    if not granular:
        n_inferred = int((singlets['detection_status'] == 'inferred').sum())
        n_multiplets = int((cell_meta['combo_type'] == 'multiplet').sum())
        return {
            'perfect_singlet': n_perfect,
            'inferred_singlet': n_inferred,
            'multiplet': n_multiplets,
            'unassigned': n_unassigned,
        }

    # --- Granular: split inferred singlets by which position is missing ---
    inferred_singlets = singlets[singlets['detection_status'] == 'inferred']
    # pos1_missing: pos1_umi == 0 (pos2 detected, pos1 inferred from pair_id)
    # pos2_missing: pos2_umi == 0 (pos1 detected, pos2 inferred from pair_id)
    if 'pos1_umi' in inferred_singlets.columns and 'pos2_umi' in inferred_singlets.columns:
        n_inf_pos1 = int((inferred_singlets['pos1_umi'] == 0).sum())
        n_inf_pos2 = int((inferred_singlets['pos2_umi'] == 0).sum())
    else:
        # Fallback: check inferred_pos columns
        has_inf_p1 = 'inferred_pos1_grna' in inferred_singlets.columns
        has_inf_p2 = 'inferred_pos2_grna' in inferred_singlets.columns
        n_inf_pos1 = 0
        n_inf_pos2 = 0
        if has_inf_p1:
            n_inf_pos1 = int((inferred_singlets['inferred_pos1_grna'].fillna('') != '').sum())
        if has_inf_p2:
            n_inf_pos2 = int((inferred_singlets['inferred_pos2_grna'].fillna('') != '').sum())

    # --- Granular: split multiplets by detection_status composition ---
    multiplet_bcs = cell_meta[cell_meta['combo_type'] == 'multiplet']['bc_wells']
    if len(multiplet_bcs) > 0:
        multi_rows = df_valid[df_valid['bc_wells'].isin(multiplet_bcs)]
        # Per-cell: check if ALL pair-rows are complete, all inferred, or mixed
        cell_status = multi_rows.groupby('bc_wells')['detection_status'].agg(
            lambda x: (
                'all_complete' if (x == 'complete').all()
                else ('all_inferred' if (x == 'inferred').all()
                      else 'mixed')
            )
        )
        n_mc = int((cell_status == 'all_complete').sum())
        n_mi = int((cell_status == 'all_inferred').sum())
        n_mm = int((cell_status == 'mixed').sum())
    else:
        n_mc = n_mi = n_mm = 0

    return {
    'multiplet_mixed': n_mm,
    'unassigned': n_unassigned,
}

def analyze_grna_ibar_pairing(umi_counts_path, ref_df, out_dir, valid_bcs=None, umi_thresholds=None, save_decoupled=True, return_df=True):
    """Estimate gRNA-iBAR decoupling rate with subcategories across UMI thresholds.

    Decoupling = the observed (grna_name, ibar_seq) pair is NOT in the reference,
    meaning the guide unit identity has been disrupted.

    Args:
        umi_counts_path: str (path to CSV) or pandas.DataFrame containing UMI counts
        ref_df: reference DataFrame
        out_dir: output directory for the summary CSV
        valid_bcs: optional list/set of valid barcodes to restrict analysis to
        umi_thresholds: single int or list of ints to filter minimum UMIs. If None,
                        defaults to a list from config or [1, 2, 3, 4, 5, 10, 20].
        save_decoupled: whether to save a csv of all decoupled guides. Only applies
                        if exactly one UMI threshold is specified.
        return_df: whether to return the detached guides dataframe. Only applies
                   if exactly one UMI threshold is specified.

    Returns:
        dict with 'summary_df' (DataFrame) containing stats per threshold, 
        and optionally 'decoupled_df' (DataFrame if applicable).
    """
    print("\n  Part B.2: gRNA-iBAR Pairing QC (Decoupling)")

    if umi_thresholds is None:
        try:
            from scripts.qc_config_loader import load_qc_config
            cfg = load_qc_config("config.yaml") # Try default config first
            umi_thresholds = cfg.get('recombination', {}).get('umi_thresholds', [1, 2, 3, 4, 5, 10, 20])
        except Exception:
            umi_thresholds = [1, 2, 3, 4, 5, 10, 20]

    if isinstance(umi_thresholds, (int, float)):
        umi_thresholds = [int(umi_thresholds)]
        
    umi_thresholds = sorted(list(set(umi_thresholds)))
    is_single_thresh = len(umi_thresholds) == 1
    
    if not is_single_thresh:
        save_decoupled = False
        return_df = False

    # Build lookups from reference
    ref_pairs = set()           
    ref_ibar_to_grna = {}       
    ref_all_ibars = set()       
    for _, row in ref_df.iterrows():
        ref_pairs.add((row['grna_name'], row['ibar_seq']))
        ref_all_ibars.add(row['ibar_seq'])
        ref_ibar_to_grna.setdefault(row['ibar_seq'], set()).add(row['grna_name'])

    if isinstance(umi_counts_path, pd.DataFrame):
        header = umi_counts_path.columns.tolist()
        is_df = True
    else:
        header = pd.read_csv(umi_counts_path, nrows=0).columns.tolist()
        is_df = False

    if 'ibar_seq_corrected' in header:
        ibar_col = 'ibar_seq_corrected'
    elif 'ibar_seq' in header:
        ibar_col = 'ibar_seq'
    else:
        print("    WARNING: No ibar column in UMI counts. Skipping.")
        return {'summary_df': pd.DataFrame()}

    needed_cols = ['bc_wells', 'grna_name', ibar_col, 'umi_count']
    if 'pair_id' in header: needed_cols.append('pair_id')
    if 'sample_name' in header: needed_cols.append('sample_name')
    
    available = [c for c in needed_cols if c in header]

    n_correct = {t: 0 for t in umi_thresholds}
    n_wrong_target = {t: 0 for t in umi_thresholds}
    n_ibar_missing = {t: 0 for t in umi_thresholds}
    n_novel_ibar = {t: 0 for t in umi_thresholds}
    
    cell_has_decoupled = {t: set() for t in umi_thresholds}
    cell_all = {t: set() for t in umi_thresholds}
    
    decoupled_records = []

    if is_df:
        chunks = [umi_counts_path[available]]
    else:
        chunk_size = 1_000_000
        chunks = pd.read_csv(umi_counts_path, usecols=available, chunksize=chunk_size)

    for chunk in chunks:
        if valid_bcs is not None:
            chunk = chunk[chunk['bc_wells'].isin(valid_bcs)]

        g_arr = chunk['grna_name'].fillna('').values
        ib_arr = chunk[ibar_col].fillna('').values
        bc_arr = chunk['bc_wells'].values
        umi_arr = chunk['umi_count'].values
        
        if save_decoupled or return_df:
            pid_arr = chunk['pair_id'].values if 'pair_id' in chunk.columns else np.full(len(chunk), '')
            sam_arr = chunk['sample_name'].values if 'sample_name' in chunk.columns else np.full(len(chunk), '')

        for i in range(len(chunk)):
            g = g_arr[i]
            ib = ib_arr[i]
            bc = bc_arr[i]
            umi = umi_arr[i]
            
            if g == '':
                continue
            
            cat = None
            if ib == '': cat = 'ibar_missing'
            elif (g, ib) in ref_pairs: cat = 'correct'
            elif ib in ref_all_ibars: cat = 'wrong_target'
            else: cat = 'novel_ibar'
            
            for t in umi_thresholds:
                if umi >= t:
                    cell_all[t].add(bc)
                    if cat == 'correct':
                        n_correct[t] += 1
                    else:
                        if cat == 'ibar_missing': n_ibar_missing[t] += 1
                        elif cat == 'wrong_target': n_wrong_target[t] += 1
                        elif cat == 'novel_ibar': n_novel_ibar[t] += 1
                        cell_has_decoupled[t].add(bc)
                
            if is_single_thresh and (save_decoupled or return_df) and cat is not None and cat != 'correct':
                t = umi_thresholds[0]
                if umi >= t:
                    grna_from_ibar = ','.join(ref_ibar_to_grna.get(ib, [])) if ib != '' else ''
                    decoupled_records.append({
                        'bc_wells': bc,
                        'grna_name': g,
                        'ibar_seq_corrected': ib,
                        'pair_id': pid_arr[i],
                        'grna_name_from_ibar': grna_from_ibar,
                        'sample_name': sam_arr[i],
                        'decoupling_category': cat,
                        'umi_count': umi
                    })

    rows = []
    for t in umi_thresholds:
        nc = n_correct[t]
        nwt = n_wrong_target[t]
        nim = n_ibar_missing[t]
        nni = n_novel_ibar[t]
        
        ntot = nc + nwt + nim + nni
        ndec = nwt + nim + nni
        ncells = len(cell_all[t])
        ncells_dec = len(cell_has_decoupled[t])
        rate = ncells_dec / ncells if ncells > 0 else 0.0
        
        rows.append({
            'umi_threshold': t,
            'total': ntot,
            'correct': nc,
            'decoupled_total': ndec,
            'wrong_target': nwt,
            'ibar_missing': nim,
            'novel_ibar': nni,
            'pct_correct': (nc / ntot * 100) if ntot > 0 else 0,
            'pct_decoupled': (ndec / ntot * 100) if ntot > 0 else 0,
            'n_cells': ncells,
            'n_cells_decoupled': ncells_dec,
            'cell_decoupling_rate': rate
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(out_dir, 'q2_decoupling_qc_summary.csv'), index=False)
    
    # Print basic output for terminal tracking (showing threshold=1 or latest)
    def _pct(n, total): return f"{n/total*100:.2f}%" if total > 0 else "N/A"
    t_print = umi_thresholds[0]  # Just print for the lowest threshold
    print(f"    Total UMI entries (with gRNA match, >= {t_print} UMIs): {summary.loc[summary['umi_threshold']==t_print, 'total'].values[0]:,}")
    print(f"    Cell-level decoupling rate at >= {t_print} UMIs:         {summary.loc[summary['umi_threshold']==t_print, 'cell_decoupling_rate'].values[0]:.4f} "
          f"({summary.loc[summary['umi_threshold']==t_print, 'n_cells_decoupled'].values[0]:,}/{summary.loc[summary['umi_threshold']==t_print, 'n_cells'].values[0]:,} cells)")
          
    if not is_single_thresh:
        print(f"    Tested {len(umi_thresholds)} UMI thresholds. Use plot_decoupling_vs_umi_threshold() to visualize.")

    ret = {'summary_df': summary}
    # For backward compatibility with existing main script, always return a default scalar rate (lowest threshold)
    ret['cell_decoupling_rate'] = summary.loc[0, 'cell_decoupling_rate']
    
    if is_single_thresh:
        decoupled_df = pd.DataFrame(decoupled_records) if decoupled_records else pd.DataFrame()
        if return_df:
            ret['decoupled_df'] = decoupled_df
        if save_decoupled and not decoupled_df.empty:
            out_file = os.path.join(out_dir, f"q2_decoupled_guides_min_umi_{umi_thresholds[0]}.csv")
            decoupled_df.to_csv(out_file, index=False)
            print(f"    Saved {len(decoupled_df):,} decoupled guides to {out_file}")

    return ret

def plot_decoupling_vs_umi_threshold(summary_df, save_path=None):
    """Plot cell decoupling rate and total valid cells vs UMI threshold."""
    if 'umi_threshold' not in summary_df.columns:
        print("    WARNING: 'umi_threshold' column missing from summary.")
        return None
        
    fig, ax1 = plt.subplots(figsize=(8, 5))
    
    color = 'tab:red'
    ax1.set_xlabel('Minimum UMI Threshold', fontsize=12)
    ax1.set_ylabel('Cell Decoupling Rate', color=color, fontsize=12)
    ax1.plot(summary_df['umi_threshold'], summary_df['cell_decoupling_rate'], 
             color=color, marker='o', linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--', alpha=0.6)
    
    # Second y-axis for total cells remaining
    ax2 = ax1.twinx()  
    color2 = 'tab:blue'
    ax2.set_ylabel('Total Valid Cells (Log10)', color=color2, fontsize=12)  
    ax2.plot(summary_df['umi_threshold'], np.log10(np.maximum(summary_df['n_cells'], 1)), 
             color=color2, marker='s', linestyle=':', linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    plt.title('Decoupling Rate vs Minimum UMI Threshold', fontsize=14)
    fig.tight_layout()  
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"    Saved plot to {save_path}")
        
    return fig


def plot_observed_vs_expected(res, save_path):
    """Plot observed vs expected category fractions. Auto-detects 4 vs 7 categories."""
    counts = res['counts']
    probs = res['probs']
    sample = res.get('sample', 'All Cells')
    N = sum(counts.values())

    if res.get('granular', False):
        cats = _V2_CATS
        labels = ['Perfect\nsinglet', 'Inf. pos1\nmissing', 'Inf. pos2\nmissing',
                  'Multi.\nall compl.', 'Multi.\nall inf.', 'Multi.\nmixed', 'Unassigned']
        figsize = (10, 5)
    else:
        cats = ['perfect_singlet', 'inferred_singlet', 'multiplet', 'unassigned']
        labels = ['Perfect\nsinglet', 'Inferred\nsinglet', 'Multiplet', 'Unassigned']
        figsize = (7, 5)

    fig, ax = plt.subplots(figsize=figsize)
    obs_frac = [counts.get(c, 0) / N for c in cats]
    exp_frac = [probs.get(c, 0) for c in cats]

    x = np.arange(len(cats))
    w = 0.35

    bars1 = ax.bar(x - w/2, obs_frac, w, label='Observed',
                    color='#3498db', alpha=0.8, edgecolor='white')
    bars2 = ax.bar(x + w/2, exp_frac, w, label='Expected',
                    color='#e74c3c', alpha=0.8, edgecolor='white')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8 if res.get('granular') else 9)
    ax.set_ylabel('Fraction of cells', fontsize=11)
    ax.set_title(f'{sample} (N={N:,})', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2, axis='y')

    for bar, frac, cat in zip(bars1, obs_frac, cats):
        n = counts.get(cat, 0)
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{n:,}', ha='center', va='bottom', fontsize=6 if res.get('granular') else 7)

    params = res['params']
    param_str = ', '.join(f'{k}={v:.3f}' for k, v in params.items())
    bic_str = f"\nBIC={res['bic']:.1f}" if 'bic' in res else ''
    ax.text(0.02, 0.97, param_str + bic_str, transform=ax.transAxes,
            fontsize=7, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {save_path}")


def run_recombination_model(adata, out_dir, moi=0.3, doublet_rate=0.03, extended=False, decoupling_rate=0.0, verbose=True, split_by_sample=False):
    """Run the recombination MLE model on the dataset, optionally split by sample."""
    print(f"\n  --- Part B.1: Recombination MLE Model (MOI={moi}, extended={extended}) ---")
    valid_bcs = set(adata.obs.index)
    results = []
    
    # Dynamically load assignment_df from out_dir parent
    output_dir = os.path.dirname(out_dir)
    assign_path = os.path.join(output_dir, 'grna_assignment', 'grna_cell_assignment.csv')
    assignment_df = pd.read_csv(assign_path)
    
    if not split_by_sample:
        counts = count_model_categories(assignment_df, valid_bcs)
        try:
            res = fit_mle(counts, moi, doublet_rate, extended=extended, decoupling_rate=decoupling_rate, verbose=verbose)
            res['sample'] = 'All Cells'
            plot_observed_vs_expected(res, os.path.join(out_dir, 'q2_observed_vs_expected.png'))
            results.append(res)
        except Exception as e:
            print(f"    WARNING: MLE failed for 'All Cells': {e}")
    else:
        samples = adata.obs['sample'].unique() if 'sample' in adata.obs.columns else ['All Cells']
        for s in samples:
            if verbose:
                print(f"\n    Sample: {s}")
            sample_bcs = set(adata.obs[adata.obs['sample'] == s].index)
            counts = count_model_categories(assignment_df, sample_bcs)
            try:
                res = fit_mle(counts, moi, doublet_rate, extended=extended, decoupling_rate=decoupling_rate, verbose=verbose)
                res['sample'] = str(s)
                safe_s = str(s).replace('/', '_').replace(' ', '_')
                plot_observed_vs_expected(res, os.path.join(out_dir, f'q2_observed_vs_expected_{safe_s}.png'))
                results.append(res)
            except Exception as e:
                print(f"    WARNING: MLE failed for sample {s}: {e}")
                
    if results:
        params_list = []
        for r in results:
            row = {'sample': r['sample']}
            row.update(r['params'])
            row['nll'] = r['nll']
            row['success'] = r['success']
            params_list.append(row)
            
        pd.DataFrame(params_list).to_csv(os.path.join(out_dir, 'q2_recombination_params.csv'), index=False)
        
    return results


# ============================================================================
# Part C: UMI Threshold Sweep
# ============================================================================

def sweep_umi_threshold(deduped_df, ref_df, valid_bcs, out_dir, cfg_recomb, thresholds=None):
    """Re-run dual-guide assignment and MLE across multiple UMI thresholds.
    
    Args:
        deduped_df: UMI deduped read counts.
        ref_df: reference DataFrame.
        valid_bcs: set of all valid single-cell barcodes.
        out_dir: Output directory
        cfg_recomb: Recombination config section (contains moi, doublet_rate, etc.)
        thresholds: List of thresholds to test. If None, uses defaults.
        
    Returns:
        DataFrame of MLE results per threshold.
    """
    print("\n  --- Part C: UMI Threshold Sweep ---")
    
    if thresholds is None:
        thresholds = cfg_recomb.get('umi_thresholds', [1, 2, 3, 4, 5, 8, 10])
        
    moi = cfg_recomb.get('moi', 0.3)
    delta = cfg_recomb.get('doublet_rate', 0.03)
    extended = cfg_recomb.get('extended_model', False)
    decoupling_rate = cfg_recomb.get('decoupling_rate', 0.0)
    
    print(f"    Sweep thresholds: {thresholds}")
    print(f"    Fixed params: MOI={moi}, delta={delta}, extended={extended}, c={decoupling_rate}")
    
    # Check iBAR column
    ibar_col = 'ibar_seq'
    if 'ibar_seq_corrected' in deduped_df.columns:
        ibar_col = 'ibar_seq_corrected'
        
    sweep_results = []
    
    for t in tqdm(thresholds, desc="    Sweeping Thresholds"):
        # 1. Assign guides
        with open(os.devnull, 'w') as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                ass_df = assign_dual_guide_pairs(deduped_df, ref_df, min_umi_threshold=t, ibar_col=ibar_col)
            finally:
                sys.stdout = old_stdout
                
        # 2. Count categories (granular)
        counts = count_model_categories(ass_df, valid_bcs, granular=True)
        
        # 3. Fit granular MLE
        try:
            res = fit_mle_v2(counts, moi, delta, verbose=False)
            row = {'threshold': t}
            row.update(res['params'])
            row.update({f'n_{c}': counts.get(c, 0) for c in _V2_CATS})
            row['nll'] = res['nll']
            row['bic'] = res['bic']
            row['success'] = res['success']
            sweep_results.append(row)
        except Exception as e:
            print(f"    WARNING: MLE failed for threshold {t}: {e}")
            
    sweep_df = pd.DataFrame(sweep_results)
    sweep_df.to_csv(os.path.join(out_dir, 'q2_threshold_sweep.csv'), index=False)
    
    # Plot recombination components
    y_cols = ['r_swap', 'r12', 'r_del1', 'r_del2']
        
    fig, ax = plt.subplots(figsize=(8, 5))
    for col in y_cols:
        if col in sweep_df.columns:
            ax.plot(sweep_df['threshold'], sweep_df[col], marker='o', label=col)
            
    ax.set_xlabel('UMI Threshold', fontsize=12)
    ax.set_ylabel('Estimated Recombination Rate', fontsize=12)
    ax.set_title('Recombination Parameter Sensitivity to UMI Threshold', fontsize=13)
    ax.set_xticks(sweep_df['threshold'])
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plot_path = os.path.join(out_dir, 'q2_threshold_sweep.png')
    plt.tight_layout()
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Plot saved to {plot_path}")
    
    return sweep_df


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="QC Q2: Recombination & Background Noise Estimation"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_qc_config(args.config)
    output_dir = cfg['paths']['output_dir']

    print(f"\n{'='*60}")
    print(f"  QC Q2: Recombination & Background Noise Estimation")
    print(f"{'='*60}")

    out_dir = os.path.join(output_dir, 'qc_analysis')
    os.makedirs(out_dir, exist_ok=True)
    
    # Inputs
    adata_path = os.path.join(output_dir, 'combined_adata.h5ad')
    assign_dir = os.path.join(output_dir, 'grna_assignment')
    assign_path = os.path.join(assign_dir, 'grna_cell_assignment.csv')
    umi_counts_path = os.path.join(assign_dir, 'grna_umi_counts.csv')
    ref_path = cfg['paths']['grna_reference']
    
    # 1. Load Data
    print(f"\n  Loading datastreams...")
    adata = sc.read_h5ad(adata_path)
    valid_bcs = set(adata.obs.index)
    print(f"    AnnData: {len(valid_bcs):,} valid cells")
    
    assignment_df = pd.read_csv(assign_path)
    umi_df = pd.read_csv(umi_counts_path)
    ref_df = pd.read_csv(ref_path)

    # --- Part A ---
    qc_cfg = cfg.get('qc', {})
    metric = qc_cfg.get('q1c_noise_metric', 'total_umi')
    category = qc_cfg.get('q1c_noise_category', 'perfect_singlet')
    noise = estimate_background_noise(adata, assignment_df, umi_df, metric=metric, cell_category=category)
    noise.to_csv(os.path.join(out_dir, 'q2_background_noise.csv'))
    
    dist_df = compute_guide_umi_distributions(assignment_df, umi_df, valid_bcs)
    if not dist_df.empty:
        dist_df.to_csv(os.path.join(out_dir, 'q2_guide_umi_distributions.csv'), index=False)
        
    # --- Part B ---
    decouple_result = analyze_grna_ibar_pairing(umi_counts_path, ref_df, out_dir, valid_bcs=valid_bcs)
    cell_decouple_rate = decouple_result['cell_decoupling_rate']
    plot_decoupling_vs_umi_threshold(decouple_result['summary_df'], os.path.join(out_dir, 'q2_decoupling_vs_umi_threshold.png'))
    
    recomb_cfg = cfg.get('recombination', {})
    moi = recomb_cfg.get('moi', 0.3)
    doublet_rate = recomb_cfg.get('doublet_rate', 0.03)
    extended = recomb_cfg.get('extended_model', False)
    decoupling_rate = recomb_cfg.get('decoupling_rate', 0.0)
    
    # --- Part B.1: Legacy 4-category model ---
    counts_legacy = count_model_categories(assignment_df, valid_bcs, granular=False)
    print("\n  --- Part B.1a: Legacy Recombination MLE (4-category) ---")
    res_legacy = fit_mle(counts_legacy, moi, doublet_rate, extended=extended, decoupling_rate=decoupling_rate)
    plot_observed_vs_expected(res_legacy, os.path.join(out_dir, 'q2_observed_vs_expected_legacy.png'))
    
    # --- Part B.2: Granular 7-category model ---
    counts_gran = count_model_categories(assignment_df, valid_bcs, granular=True)
    print("\n  --- Part B.1b: Granular Recombination MLE (7-category) ---")
    res_gran = fit_mle_v2(counts_gran, moi, doublet_rate)
    res_gran['sample'] = 'All Cells'
    plot_observed_vs_expected(res_gran, os.path.join(out_dir, 'q2_observed_vs_expected_granular.png'))
    
    # Print decoupling rate alongside r_swap for interpretation
    r_swap_hat = res_gran['params'].get('r_swap', 0)
    print(f"\n  --- Swap vs Decoupling Comparison ---")
    print(f"    MLE r_swap (swap + decoupling):       {r_swap_hat:.4f}")
    print(f"    Cell-level decoupling rate (from B.2): {cell_decouple_rate:.4f}")
    print(f"    Estimated true recombination swap:     {max(0, r_swap_hat - cell_decouple_rate):.4f}")
    
    # Save both results
    params_row = {'model': 'granular'}
    params_row.update(res_gran['params'])
    params_row['nll'] = res_gran['nll']
    params_row['bic'] = res_gran['bic']
    params_row['success'] = res_gran['success']
    params_row['cell_decoupling_rate'] = cell_decouple_rate
    legacy_row = {'model': 'legacy'}
    legacy_row.update(res_legacy['params'])
    legacy_row['nll'] = res_legacy['nll']
    legacy_row['success'] = res_legacy['success']
    pd.DataFrame([params_row, legacy_row]).to_csv(
        os.path.join(out_dir, 'q2_recombination_params.csv'), index=False)
    
    # --- Part C ---
    if recomb_cfg.get('run_threshold_sweep', True):
        sweep_umi_threshold(umi_df, ref_df, valid_bcs, out_dir, recomb_cfg)

    print(f"\n{'='*60}")
    print(f"  QC Q2 Complete")
    print(f"  All results saved to: {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
