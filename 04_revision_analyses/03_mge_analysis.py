"""
03_comment4_mge_analysis.py
Reproducibility script for Reviewer 1 Comment 4.
Analyzes Mobile Genetic Elements (MGEs) including intI1, Insertion Sequences (IS),
and Plasmids from the curated MGE pipeline output, and correlates them with
total ARG burden and Stringent response module score across the 18 wastewater metagenomes.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# 1. Load sample table
sample_table_path = os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "tables", "full_sample_table_with_day.csv")
samples_df = pd.read_csv(sample_table_path, index_col=0)
samples = list(samples_df.index)

# 2. Load MGE data
mge_base = "/home/ihtishamnaeem36/Documents/1. metagenomics/6 final metagenomic project intergrated analysis across all cities/3 Curated databases via NR protein/MGE (Transfer Element)/MAM_asem_annojLrEjqQqjzzk3jzoQYWjNN_result"

def load_and_fix(fname):
    df = pd.read_csv(os.path.join(mge_base, fname), sep="\t")
    df = df.rename(columns={"SSLW3": "SCW2"})
    if "MHW3" in df.columns:
        df = df.drop(columns=["MHW3"])
    return df

l1_df = load_and_fix("All.MGE.Level1.txt")
l2_df = load_and_fix("All.MGE.Level2.txt")
unstrat_df = load_and_fix("All.MGE.abundance_unstratified.tsv")

# Extract MGE metrics
mge_total = l1_df[samples].sum(axis=0).astype(float)

inti1_rows = unstrat_df[unstrat_df.iloc[:, 0].astype(str).str.contains("intI1", case=False, regex=True)]
inti1_abundance = inti1_rows[samples].sum(axis=0).astype(float)

plasmid_rows = l2_df[l2_df.iloc[:, 0].astype(str).str.startswith(("Inc", "Col"))]
plasmid_abundance = plasmid_rows[samples].sum(axis=0).astype(float)

is_rows = l1_df[l1_df.iloc[:, 0].astype(str).str.startswith("IS")]
is_abundance = is_rows[samples].sum(axis=0).astype(float)

stringent_scores = samples_df["stringent_score"]
arg_burden = samples_df["total_ARG_burden"]

# Build results table
records = [
    {
        "MGE Category": "Total MGE Burden",
        "Description": "Sum across all Level 1 transfer elements (CPC)",
        "rho_vs_ARG": spearmanr(mge_total, arg_burden)[0],
        "p_vs_ARG": spearmanr(mge_total, arg_burden)[1],
        "rho_vs_Stringent": spearmanr(mge_total, stringent_scores)[0],
        "p_vs_Stringent": spearmanr(mge_total, stringent_scores)[1]
    },
    {
        "MGE Category": "Class 1 Integron Integrase (intI1)",
        "Description": "Clinical integron integrase gene copies (CPC)",
        "rho_vs_ARG": spearmanr(inti1_abundance, arg_burden)[0],
        "p_vs_ARG": spearmanr(inti1_abundance, arg_burden)[1],
        "rho_vs_Stringent": spearmanr(inti1_abundance, stringent_scores)[0],
        "p_vs_Stringent": spearmanr(inti1_abundance, stringent_scores)[1]
    },
    {
        "MGE Category": "Plasmids",
        "Description": "Inc and Col replicon families (Level 2, CPC)",
        "rho_vs_ARG": spearmanr(plasmid_abundance, arg_burden)[0],
        "p_vs_ARG": spearmanr(plasmid_abundance, arg_burden)[1],
        "rho_vs_Stringent": spearmanr(plasmid_abundance, stringent_scores)[0],
        "p_vs_Stringent": spearmanr(plasmid_abundance, stringent_scores)[1]
    },
    {
        "MGE Category": "Insertion Sequences (IS)",
        "Description": "IS families (IS26, IS10, IS200, IS621, etc., CPC)",
        "rho_vs_ARG": spearmanr(is_abundance, arg_burden)[0],
        "p_vs_ARG": spearmanr(is_abundance, arg_burden)[1],
        "rho_vs_Stringent": spearmanr(is_abundance, stringent_scores)[0],
        "p_vs_Stringent": spearmanr(is_abundance, stringent_scores)[1]
    }
]

res_df = pd.DataFrame(records)
out_csv = os.path.join(PROJECT_ROOT, "Revesion", "tables", "comment4_mge_analysis.csv")
res_df.to_csv(out_csv, index=False)
print(f"Saved {out_csv}")
print(res_df.to_string(index=False))
