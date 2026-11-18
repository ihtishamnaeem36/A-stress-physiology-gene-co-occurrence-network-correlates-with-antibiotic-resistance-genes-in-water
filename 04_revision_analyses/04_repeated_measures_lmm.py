"""
04_reviewer3_comment1_repeated_measures_lmm.py
Reproducibility script for Reviewer 3 Comment 1 and Reviewer 4 Major Issue 6.
Evaluates repeated-measures temporal structure (9 sites x 2 rounds = 18 samples):
1. Site-level aggregation (Day 1 + Day 30 average, n=9 strictly independent sites).
2. Separate round evaluations (Day 1 alone vs Day 30 alone, n=9 each).
3. Linear Mixed-Effects Model (LMM with random site intercepts via REML).
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import statsmodels.api as sm
import statsmodels.formula.api as smf

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# 1. Load sample table
sample_table_path = os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "tables", "full_sample_table_with_day.csv")
df = pd.read_csv(sample_table_path, index_col=0)

# 2. Site-level average (n=9)
site_avg = df.groupby("site")[["stringent_score", "total_ARG_burden"]].mean()
rho_site, p_rho_site = spearmanr(site_avg["stringent_score"], site_avg["total_ARG_burden"])
r_site, p_r_site = pearsonr(site_avg["stringent_score"], site_avg["total_ARG_burden"])

# 3. Round 1 (Day 1) and Round 2 (Day 30) alone (n=9 each)
d1 = df[df["day"] == "Day1"]
rho_d1, p_d1 = spearmanr(d1["stringent_score"], d1["total_ARG_burden"])

d30 = df[df["day"] == "Day30"]
rho_d30, p_d30 = spearmanr(d30["stringent_score"], d30["total_ARG_burden"])

# 4. Linear Mixed-Effects Model: total_ARG_burden ~ stringent_score + (1 | site)
lmm = smf.mixedlm("total_ARG_burden ~ stringent_score", df, groups=df["site"])
lmm_fit = lmm.fit()

# Standardized LMM
df_std = df.copy()
df_std["total_ARG_burden"] = (df["total_ARG_burden"] - df["total_ARG_burden"].mean()) / df["total_ARG_burden"].std()
df_std["stringent_score"] = (df["stringent_score"] - df["stringent_score"].mean()) / df["stringent_score"].std()
lmm_std = smf.mixedlm("total_ARG_burden ~ stringent_score", df_std, groups=df_std["site"])
lmm_std_fit = lmm_std.fit()

# Output results table
records = [
    {
        "Analysis / Test": "Site-level average (n=9 independent sites)",
        "Sample Size": 9,
        "Statistic": f"Spearman rho = {rho_site:.4f}, Pearson r = {r_site:.4f}",
        "p_value": p_rho_site,
        "Details": "Day1 and Day30 averaged per physical site"
    },
    {
        "Analysis / Test": "Round 2 alone (Day 30 only)",
        "Sample Size": 9,
        "Statistic": f"Spearman rho = {rho_d30:.4f}",
        "p_value": p_d30,
        "Details": "Independent single-timepoint cross-section"
    },
    {
        "Analysis / Test": "Round 1 alone (Day 1 only)",
        "Sample Size": 9,
        "Statistic": f"Spearman rho = {rho_d1:.4f}",
        "p_value": p_d1,
        "Details": "Independent single-timepoint cross-section"
    },
    {
        "Analysis / Test": "Linear Mixed-Effects Model (LMM)",
        "Sample Size": 18,
        "Statistic": f"Fixed effect beta = {lmm_fit.params['stringent_score']:.3f} (z = {lmm_fit.tvalues['stringent_score']:.3f})",
        "p_value": lmm_fit.pvalues["stringent_score"],
        "Details": "Random site intercept (REML, 9 groups, 2 obs/group)"
    },
    {
        "Analysis / Test": "Standardized LMM Fixed Effect",
        "Sample Size": 18,
        "Statistic": f"Standardized beta = {lmm_std_fit.params['stringent_score']:.3f} (z = {lmm_std_fit.tvalues['stringent_score']:.3f})",
        "p_value": lmm_std_fit.pvalues["stringent_score"],
        "Details": "Standardized predictor and response"
    }
]

out_df = pd.DataFrame(records)
out_csv = os.path.join(PROJECT_ROOT, "Revesion", "tables", "reviewer3_comment1_repeated_measures.csv")
out_df.to_csv(out_csv, index=False)
print(f"Saved {out_csv}")
print(out_df.to_string(index=False))
