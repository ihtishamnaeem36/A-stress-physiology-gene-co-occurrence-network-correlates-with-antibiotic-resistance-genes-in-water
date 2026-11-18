"""
02_comment2_normalization_reconciliation.py
Reproducibility script for Reviewer 1 Comment 2 / Reviewer 4 Major Issue 3.
Evaluates the impact of normalizing Stringent response module score by total KO sequencing depth
versus raw sum-of-logs against SARG Copies Per Cell (CPC).
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# 1. Load sample table
sample_table_path = os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "tables", "full_sample_table_with_day.csv")
samples_df = pd.read_csv(sample_table_path, index_col=0)
samples = list(samples_df.index)

# 2. Load unstratified KEGG profile
kegg_path = os.path.join(PROJECT_ROOT, "DATA used", "KEGG_abundance_mixed.xls")
kegg_df = pd.read_csv(kegg_path, sep="\t")
kegg_df = kegg_df.rename(columns={"SSLW3": "SCW2"})
if "MHW3" in kegg_df.columns:
    kegg_df = kegg_df.drop(columns=["MHW3"])

unstrat = kegg_df[~kegg_df["KEGG"].astype(str).str.contains("|", regex=False)]
total_ko_abundance = unstrat[samples].sum(axis=0) # raw sum of TPMs in sample
total_log_abundance = np.log10(unstrat[samples] + 1).sum(axis=0)

stringent_scores = samples_df["stringent_score"]
arg_burden_cpc = samples_df["total_ARG_burden"]

# Primary
r_prim, p_prim = spearmanr(stringent_scores, arg_burden_cpc)

# Depth-normalized (relative proportion in KO pool)
stringent_rel = stringent_scores / total_ko_abundance
r_rel, p_rel = spearmanr(stringent_rel, arg_burden_cpc)

# Log-depth normalized
stringent_rel_log = stringent_scores / total_log_abundance
r_rel_log, p_rel_log = spearmanr(stringent_rel_log, arg_burden_cpc)

# CLR transformed
clr_table_path = os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "tables", "module_scores_CLR_vs_ARG_burden.csv")
clr_df = pd.read_csv(clr_table_path)
clr_strin = clr_df.loc[clr_df["module"] == "Stringent response (ppGpp)"].iloc[0]

results = [
    {"Normalization Scheme": "Primary (Sum-of-logs TPM vs ARG CPC)", "rho": r_prim, "p_value": p_prim, "Details": "Original reported model (n=18)"},
    {"Normalization Scheme": "Depth-normalized (Stringent Score / Total KO TPM vs ARG CPC)", "rho": r_rel, "p_value": p_rel, "Details": "Strictly relative KO fraction in sample"},
    {"Normalization Scheme": "Log-depth normalized (Stringent Score / Total Log Abundance vs ARG CPC)", "rho": r_rel_log, "p_value": p_rel_log, "Details": "Normalized by total log-scale capacity"},
    {"Normalization Scheme": "Centered Log Ratio (CLR transformed across 903 KOs vs ARG CPC)", "rho": clr_strin["rho"], "p_value": clr_strin["p"], "Details": f"Compositional transformation (FDR q = {clr_strin['q']:.3f})"}
]

out_df = pd.DataFrame(results)
out_csv = os.path.join(PROJECT_ROOT, "Revesion", "tables", "comment2_normalization_reconciliation.csv")
out_df.to_csv(out_csv, index=False)
print(f"Saved {out_csv}")
print(out_df.to_string(index=False))
