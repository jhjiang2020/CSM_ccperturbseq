# CSM ccPerturbSeq — QC Analysis Overview

This document summarizes the Quality Control (QC) pipeline designed for the ccPerturbSeq dataset. The analysis is orchestrated via the `scripts/qc_analysis.ipynb` notebook and is divided into four major analytical modules (Q1–Q4).

## Environment setup

Q1–Q4, the GEX carry-over QC, and the analysis notebooks share the tracked
`csm_analysis` environment:

```bash
conda env create -f environment_qc.yml
# For an existing environment:
conda env update -n csm_analysis -f environment_qc.yml --prune

conda activate csm_analysis
python3 -c 'import edlib, matplotlib, pandas, pysam, scanpy, scipy, yaml'
```

Set `qc_conda_env: "csm_analysis"` in QC configurations and set
`environment.conda_prefix` in the referenced main config to the output of
`conda info --base`. UGER wrappers activate this environment inside compute
jobs; it does not need to remain activated on the login node.

---

## Q1: Bulk Guide-Seq Library Analysis (`qc_q1_guideseq.py`)

**Objective:**
Analyze pre-lentiviral bulk plasmid library sequencing to resolve guide unit identities, detect construct defects (dropout, swaps, deletions), and quantify per-guide and per-combo frequencies for downstream fitness analysis (Q4).

> **Note:** Q1 runs **standalone**, independent of the single-cell pipeline. It uses its own config (`qc_config_guideseq.yaml`) and does not require `combined_adata.h5ad`.

**Methodology:**
1.  **FastQ Processing & Guide Resolution:** Merges paired-end reads (R1 and R2rc) using the `tRNA_leader` as an anchor.
2.  **Multi-Step Disambiguation Logic:** Resolves guide unit identities through a rigorous sequential matching process: first matching the gRNA sequence, then leveraging internal stem architectures, and finally checking the iBAR to overcome overlapping guide ambiguities.
3.  **Intactness & Swap Classification:** Accurately classifies construct outcomes by requiring successful resolution of both construct positions to a valid pair ID for an "intact" baseline, and actively tracks structural swaps based on resolved identities across `pos1` and `pos2`.
4.  **Frequency Export:** `get_guide_frequencies()` and `get_combo_frequencies()` produce per-guide and per-combo frequency tables (ignoring recombination/deletion reads), feeding Q4 fitness analysis.

---

## Q2: Recombination & Background Noise Estimation (`qc_q2_recombination.py`)

**Objective:**
Estimate construct recombination events and ambient guide noise from single-cell data, providing foundational quality metrics for cell categorization and recombination-aware analysis.

**Methodology:**

### Part A — Ambient Noise Estimation & Guide UMI Distributions
1.  **Background Noise:** `estimate_background_noise()` computes per-cell off-target UMI statistics (total_umi, n_guides, or mean_umi) from perfect singlet or unassigned cell populations.
2.  **Guide UMI Distributions:** `compute_guide_umi_distributions()` estimates mean and std of log10(UMI) per guide unit from perfect singlets, producing inputs for optional GMM-based guide assignment refinement.

### Part B — Recombination Model (MLE) & Decoupling QC
1.  **Multinomial MLE Model:** Fits theoretically predicted observable category proportions against empirically observed cell counts, adjusting for MOI and split-pool doublet rates.
    *   **Legacy Model (4-category):** Fits `perfect_singlet`, `inferred_singlet`, `multiplet`, `unassigned` to estimate global dropout ($d$), silencing ($s$), and aggregate recombination ($r$).
    *   **Granular Model (7-category):** Expands observables to `perfect_singlet`, `inferred_singlet_pos1_missing`, `inferred_singlet_pos2_missing`, `multiplet_all_complete` (MOI), `multiplet_all_inferred` (swaps), `multiplet_mixed`, and `unassigned`. Accurately resolves recombination into swap ($r_\text{swap}$), 1→2 transfer ($r_{12}$), and position-specific deletion ($r_{\text{del}1}$, $r_{\text{del}2}$).
2.  **Decoupling QC:** `analyze_grna_ibar_pairing()` checks gRNA-iBAR pairing correctness as a post-MLE quality check, detecting decoupling events that could inflate recombination estimates.

### Part C — UMI Threshold Sweep
`sweep_umi_threshold()` re-runs dual-guide pair assignment at multiple UMI thresholds, fits the granular recombination model at each, and plots parameter sensitivity ($r_\text{swap}$, $r_{12}$, $r_{\text{del}1}$, $r_{\text{del}2}$) to the threshold choice.

---

## Q3: KD Efficacy + Inferred Source Classification (`qc_q3_kd_efficacy.py`)

**Objective:**
Determine how effectively each guide suppresses its designated target gene, classify inferred singlets by source, and validate that dropout-inferred cells show comparable KD efficacy to perfect singlets.

**Methodology:**

### Q3a — Gold-Standard Perfect Singlet Efficacy
1.  **Strict Control Matching:** Isolates `gRNA_ctrl` cells (one targeting guide + one ctrl) and `gRNA_nonctrl` cells (two active targeting guides).
2.  **Baseline Definition:** All comparisons use `ctrl_ctrl` singlet cell populations (cells expressing two non-targeting controls).
3.  **Metrics:** log2FC, Wilcoxon p-value, Cohen's d, pseudobulk CPM, % non-zero expression.

### Q3b — Validation Across Detection Modes
Evaluates knockdown consistency in `perfect_all`, `detected_in_inferred`, and `inferred_in_inferred` cell populations against the shared `ctrl_ctrl` baseline.

### Q3c — Inferred Source Classification + Source 1 DE Validation
1.  `classify_inferred_sources()` uses Q2's noise threshold to classify inferred singlets into:
    *   **Source 1 (Dropout):** No residual off-target signal at the missing position after noise masking → true biological dropout.
    *   **Source 2 (Recombinant/Noise):** Residual signal remains → likely recombination or ambient noise artifact.
2.  **Source 1 DE Validation:** Repeats KD efficacy tests exclusively on Source 1 (dropout) inferred singlets to confirm they show comparable knockdown to perfect singlets.

### Utilities
*   `compare_expression()` — General-purpose expression comparison between any two cell groups (heatmap, violin, or both).
*   `compute_custom_combo_de()` — Multi-combo × multi-gene DE analysis at gene-combo and guide-pair levels.

---

## Q4: Fitness Analysis (`qc_q4_fitness.py`)

**Objective:**
Compare guide-combo frequencies between pre-lentiviral library (from Q1) and single-cell data to estimate perturbation-driven survival/fitness effects.

**Methodology:**
1.  **Cell-Level Fitness (Part A):** `analyze_fitness()` computes `log2FC = cell_frequency / library_frequency` for each target combo, identifying enriched and depleted perturbations post-selection.
2.  **Novel Combo Detection (Part B):** `identify_novel_combos()` finds guide combinations observed in scRNA-seq perfect singlets but absent from the pre-lentiviral library, indicating possible lentiviral recombination artifacts.
3.  **Fitness Estimation (Part C):** `estimate_fitness()` computes per-combo and per-target-gene fitness effects (`log2FC(scRNA-seq freq / pre-lenti freq)`) with aggregate statistics and ranked visualizations.
