"""
09_leave_top_subtype_out_sensitivity.py
Reproducibility script for Reviewer 4 Major Issue 7 (Supplementary Table S16).
Evaluates whether type-level associations between Stringent response and the 7 SARG
antibiotic resistance types are driven by single dominant gene families.
Performs leave-top-subtype-out sensitivity tests for all multi-subtype resistance types.
"""
import os, sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

sarg_candidates = [
    os.path.join(PROJECT_ROOT, "DATA used", "data SARG (resistance gene) (best for wasterwter metagenomic/MAM_asem_anno22fHiA9hpGDXDXTLWiNAPm_result", "All.SARG.Subtype.txt"),
    os.path.join(os.path.expanduser("~"), "Documents", "1. metagenomics", "thesis", "project 2.2", "DATA used", "data SARG (resistance gene) (best for wasterwter metagenomic/MAM_asem_anno22fHiA9hpGDXDXTLWiNAPm_result", "All.SARG.Subtype.txt"),
    "All.SARG.Subtype.txt"
]

sarg_path = None
for c in sarg_candidates:
    if os.path.exists(c):
        sarg_path = c
        break

sample_table_candidates = [
    os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "tables", "full_sample_table_with_day.csv"),
    os.path.join(os.path.expanduser("~"), "Documents", "1. metagenomics", "thesis", "project 2.2", "10. Final Manuscript (Network-Resistome Integration)", "tables", "full_sample_table_with_day.csv"),
    "full_sample_table_with_day.csv"
]

sample_table_path = None
for c in sample_table_candidates:
    if os.path.exists(c):
        sample_table_path = c
        break

def run_leave_top_subtype_out():
    if not sarg_path or not os.path.exists(sarg_path):
        print(f"Warning: {sarg_path} not found. Please provide All.SARG.Subtype.txt.")
        return
        
    df_sub = pd.read_csv(sarg_path, sep='	')
    df_sub = df_sub.rename(columns={'SSLW3': 'SCW2'}).drop(columns=['MHW3'], errors='ignore')

    samples_df = pd.read_csv(sample_table_path, index_col=0)
    samples = list(samples_df.index)
    stringent = samples_df['stringent_score']

    df_sub['Type'] = df_sub['# SARG Subtype'].apply(lambda x: x.split('__')[0])
    df_sub['Subtype'] = df_sub['# SARG Subtype'].apply(lambda x: x.split('__')[1])

    types = ['multidrug', 'sulfonamide', 'aminoglycoside', 'chloramphenicol', 'beta_lactam', 'trimethoprim']
    records = []

    print("=== Leave-Top-Subtype-Out Sensitivity Analysis (Table S16) ===")
    for t in types:
        sub_df = df_sub[df_sub['Type'] == t].copy()
        sub_df['mean_abund'] = sub_df[samples].mean(axis=1)
        sub_df = sub_df.sort_values('mean_abund', ascending=False)
        top_subtype = sub_df.iloc[0]['Subtype']
        
        # Full type correlation
        full_abund = sub_df[samples].sum(axis=0)
        rho_full, p_full = spearmanr(stringent, full_abund)
        
        # Leave top subtype out
        sub_no_top = sub_df.iloc[1:][samples].sum(axis=0)
        rho_no_top, p_no_top = spearmanr(stringent, sub_no_top)
        
        records.append({
            "Resistance_Type": t,
            "Subtype_Count": len(sub_df),
            "Top_Subtype": top_subtype,
            "Full_Spearman_rho": round(rho_full, 3),
            "Full_p_value": round(p_full, 4),
            "Leave_Top_Out_rho": round(rho_no_top, 3),
            "Leave_Top_Out_p_value": round(p_no_top, 4)
        })
        print(f"{t:15} | Top: {top_subtype:12} | Full rho={rho_full:.3f} (p={p_full:.4f}) | Leave-out rho={rho_no_top:.3f} (p={p_no_top:.4f})")

    res_df = pd.DataFrame(records)
    out_csv = os.path.join(os.path.dirname(__file__), "leave_top_subtype_out_results.csv")
    res_df.to_csv(out_csv, index=False)

if __name__ == '__main__':
    run_leave_top_subtype_out()
