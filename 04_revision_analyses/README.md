# 04_revision_analyses — Peer Review Sensitivity Analyses and Robustness Checks

This directory contains every reproducibility and validation script developed in response to peer review critique for:
**"A stress-physiology gene co-occurrence network correlates with antibiotic-resistance genes in wastewater metagenomes"**
(Journal of Environmental Chemical Engineering, Ref: JECE-D-26-19867).

These scripts address specific questions raised regarding temporal replication, background network null models, predictive cross-validation boundaries, normalisation frameworks, mobile genetic elements, taxonomic contributions, and subtype distributions.

---

## Script Inventory & Reviewer Cross-Reference

| Script | Reviewer Comment | Manuscript / Supp Item | Key Finding |
| :--- | :--- | :--- | :--- |
| **`01_partial_corr_and_permutation.py`** | Reviewer 1 Comment 1 / Reviewer 3 Comment 2 | Section 3.9, Table S11 | Partial correlation controlling for total functional-gene abundance remains significant (rho = 0.70, p = 0.0018); Stringent response ranks at the 99.4th percentile against 1,000 size-matched background KO null sets (p_emp = 0.006). |
| **`02_normalization_reconciliation.py`** | Reviewer 1 Comment 2 / Reviewer 4 Major 3 | Section 2.8, Table S12 | Normalising module scores by total sequencing depth yields equivalent correlation (rho = 0.81, p = 0.00005), confirming that TPM vs 16S copies/cell normalisation mismatch does not alter rank order. |
| **`03_mge_analysis.py`** | Reviewer 1 Comment 4 | Section 3.8, Table S13 | Class 1 integron integrase (intI1) correlates with total ARG burden (rho = 0.73, p = 0.0006) and Stringent response (rho = 0.54, p = 0.021); insertion sequences and plasmids correlate with resistance burden. |
| **`04_repeated_measures_lmm.py`** | Reviewer 3 Comment 1 / Reviewer 4 Major 6 | Section 2.6, Section 3.7 | Site-level aggregation across the 9 independent physical sites confirms significant correlation (rho = 0.850, p = 0.0037); repeated-measures Linear Mixed-Effects Model (LMM) with random site intercepts confirms robust coupling (p < 0.0001). |
| **`05_network_null_models.py`** | Reviewer 3 Comment 3 / Reviewer 4 Major 4 | Section 3.1, Figure S3 | Evaluates full 903-KO co-occurrence network against 1,000 Erdős–Rényi random graphs and 1,000 degree-preserving Maslov–Sneppen rewired null graphs; empirical network exhibits significantly higher clustering coefficient (C = 0.584 vs 0.054 and 0.218, p < 0.001). |
| **`06_loso_elasticnet.py`** | Reviewer 3 Comment 4 / Reviewer 4 Major 8 | Section 3.6, Table S14 | Strict Leave-One-Site-Out (LOSO) cross-validation of 15-module Elastic Net yields R^2 = 0.471 (95% CI: -0.073 to 0.667) and rho = 0.684 (95% CI: 0.282 to 0.885), outperforming the univariate Stringent baseline (R^2 = 0.392, rho = 0.579). |
| **`07_clr_elasticnet.py`** | Reviewer 4 Major Issue 1(c) | Section 3.6, Section 3.9 | Under Centered-Log-Ratio (CLR) transformation where library-size scale is removed, 15-module Elastic Net retains out-of-sample predictive performance under LOSO (R^2 = 0.621, rho = 0.880). |
| **`08_taxonomic_stratification.py`** | Reviewer 4 Major Issue 5 | Section 3.7, Table S15 | Lineage-specific recalculation of Stringent response scores across 10 individual phyla and classes demonstrates that no single bacterial taxon explains the association; coupling reflects community-wide coordination. |
| **`09_leave_top_subtype_out_sensitivity.py`** | Reviewer 4 Major Issue 7 | Section 3.7, Table S16 | Leave-top-subtype-out tests across all multi-subtype ARG types confirm that associations are not driven by single dominant genes (all leave-out rho = 0.57 to 0.75, all p <= 0.0142). Bacitracin is explicitly documented as a single-subtype class (bacA). |

---

## Execution Instructions

All scripts require Python >= 3.8 and dependencies listed in `requirements.txt`:
```bash
python 01_partial_corr_and_permutation.py
python 02_normalization_reconciliation.py
python 03_mge_analysis.py
python 04_repeated_measures_lmm.py
python 05_network_null_models.py
python 06_loso_elasticnet.py
python 07_clr_elasticnet.py
python 08_taxonomic_stratification.py
python 09_leave_top_subtype_out_sensitivity.py
```
