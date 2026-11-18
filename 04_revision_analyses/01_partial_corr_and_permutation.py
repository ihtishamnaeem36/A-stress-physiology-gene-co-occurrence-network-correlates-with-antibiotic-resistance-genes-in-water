"""
01_comment1_partial_corr_and_permutation.py
Reproducibility script for Reviewer 1 Comment 1 / Reviewer 3 Comment 2.
1. Computes total functional-gene abundance across all ~10,389 unstratified KOs.
2. Computes partial Spearman correlation between Stringent response and total ARG burden controlling for total abundance.
3. Fits standardized OLS regression of ARG burden ~ Stringent Score + Total Functional Abundance.
4. Generates 1,000 size-matched (83 KOs) random KO sets from the background metagenome to compute empirical null distribution.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t
import statsmodels.api as sm

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
total_log_abundance = np.log10(unstrat[samples].astype(float) + 1).sum(axis=0)

stringent_scores = samples_df["stringent_score"]
arg_burden = samples_df["total_ARG_burden"]

# 3. Correlations
r_strin_arg, p_strin_arg = spearmanr(stringent_scores, arg_burden)
r_strin_tot, p_strin_tot = spearmanr(stringent_scores, total_log_abundance)
r_arg_tot, p_arg_tot = spearmanr(arg_burden, total_log_abundance)

# Partial spearman correlation
rx = stringent_scores.rank()
ry = arg_burden.rank()
rz = total_log_abundance.rank()

r_xy = np.corrcoef(rx, ry)[0, 1]
r_xz = np.corrcoef(rx, rz)[0, 1]
r_yz = np.corrcoef(ry, rz)[0, 1]

r_xy_z = (r_xy - r_xz * r_yz) / np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
df = len(samples) - 3
t_stat = r_xy_z * np.sqrt(df / (1 - r_xy_z**2))
p_val_partial = 2 * (1 - t.cdf(abs(t_stat), df))

# Standardized OLS
X = pd.DataFrame({"Total_Functional_Abundance": total_log_abundance, "Stringent_Score": stringent_scores})
X_norm = (X - X.mean()) / X.std()
X_norm = sm.add_constant(X_norm)
y_norm = (arg_burden - arg_burden.mean()) / arg_burden.std()
ols_model = sm.OLS(y_norm, X_norm).fit()

# Permutation test (1,000 random 83-KO sets)
np.random.seed(42)
ko_mat = unstrat[samples].astype(float)
rhos = []
for _ in range(1000):
    rand_idx = np.random.choice(ko_mat.shape[0], size=83, replace=False)
    rand_score = np.log10(ko_mat.iloc[rand_idx] + 1).sum(axis=0)
    r, _ = spearmanr(rand_score, arg_burden)
    rhos.append(r)

rhos = np.array(rhos)
p_empirical = np.mean(rhos >= r_strin_arg)

# Save results table
results = [
    {"Test": "Primary Stringent vs ARG burden", "rho": r_strin_arg, "p_value": p_strin_arg, "Details": "Unadjusted Spearman"},
    {"Test": "Total Functional Abundance vs ARG burden", "rho": r_arg_tot, "p_value": p_arg_tot, "Details": "Full ~10,389 unstratified KOs"},
    {"Test": "Stringent vs Total Functional Abundance", "rho": r_strin_tot, "p_value": p_strin_tot, "Details": "Abundance scale coupling"},
    {"Test": "Partial Spearman (Stringent vs ARG | Total Abundance)", "rho": r_xy_z, "p_value": p_val_partial, "Details": "Controlling for full functional profile"},
    {"Test": "Multiple OLS: Stringent Score Beta", "rho": np.nan, "p_value": ols_model.pvalues["Stringent_Score"], "Details": f"Beta = {ols_model.params['Stringent_Score']:.3f}"},
    {"Test": "Multiple OLS: Total Abundance Beta", "rho": np.nan, "p_value": ols_model.pvalues["Total_Functional_Abundance"], "Details": f"Beta = {ols_model.params['Total_Functional_Abundance']:.3f}"},
    {"Test": "Random 83-KO Permutation Null (n=1000)", "rho": np.mean(rhos), "p_value": p_empirical, "Details": f"Mean rho = {np.mean(rhos):.3f}, 95th pct = {np.percentile(rhos, 95):.3f}"}
]
out_df = pd.DataFrame(results)
out_csv = os.path.join(PROJECT_ROOT, "Revesion", "tables", "comment1_robustness_results.csv")
out_df.to_csv(out_csv, index=False)
print(f"Saved {out_csv}")
