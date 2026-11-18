"""
08_taxonomic_stratification.py
Reproducibility script for Reviewer 4 Major Issue 5 (Supplementary Table S15).
Evaluates taxonomic stratification across 10 lineages (5 phyla, 5 classes):
Recalculates Stringent response module scores within individual clades and tests
Spearman correlation against total ARG burden across the 18 wastewater metagenomes.
"""
import os, sys
import pandas as pd
import numpy as np
import scipy.stats as stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# Try local or mounted paths
extract_candidates = [
    os.path.join(PROJECT_ROOT, "10. Final Manuscript (Network-Resistome Integration)", "data", "stringent_stratified_extract.tsv"),
    os.path.join(os.path.expanduser("~"), "Documents", "1. metagenomics", "thesis", "project 2.2", "10. Final Manuscript (Network-Resistome Integration)", "data", "stringent_stratified_extract.tsv"),
    "stringent_stratified_extract.tsv"
]

extract_path = None
for c in extract_candidates:
    if os.path.exists(c):
        extract_path = c
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

SAMPLE_COLS_18 = ['MCW1','MCW2','MHW1','MHW2','MSLW1','MSLW2','PCW1','PCW2',
                   'PHW1','PHW2','PSLW1','PSLW2','SCW1','SCW2','SHW1','SHW2',
                   'SSLW1','SSLW2']

def run_taxonomic_stratification():
    if not extract_path or not os.path.exists(extract_path):
        print(f"Warning: {extract_path} not found. Please provide stringent_stratified_extract.tsv.")
        return
    
    df = pd.read_csv(extract_path, sep='	')
    df['KO'] = df['KEGG'].str.split('|').str[0]
    df['taxon'] = df['KEGG'].str.split('|', n=1).str[1]
    df = df.rename(columns={'SSLW3': 'SCW2'}).drop(columns=['MHW3'], errors='ignore')

    samples_df = pd.read_csv(sample_table_path, index_col=0)
    arg = samples_df.loc[SAMPLE_COLS_18, 'total_ARG_burden']

    strat = df[df['taxon'].notna()].copy()
    classified = strat[strat['taxon'] != 'unclassified'].copy()

    def get_rank(t, prefix):
        for p in t.split('.'):
            if p.startswith(prefix):
                return p.replace(prefix, '')
        return 'unresolved'

    classified['phylum'] = classified['taxon'].apply(lambda x: get_rank(x, 'p__'))
    classified['class'] = classified['taxon'].apply(lambda x: get_rank(x, 'c__'))

    clades = [
        ("Phylum", "Pseudomonadota"),
        ("Phylum", "Bacteroidota"),
        ("Phylum", "Bacillota"),
        ("Phylum", "Actinomycetota"),
        ("Phylum", "Thermodesulfobacteriota"),
        ("Class", "Betaproteobacteria"),
        ("Class", "Gammaproteobacteria"),
        ("Class", "Bacteroidia"),
        ("Class", "Flavobacteriia"),
        ("Class", "Actinomycetes"),
    ]

    records = []
    for level, name in clades:
        col = 'phylum' if level == 'Phylum' else 'class'
        sub_df = classified[classified[col] == name]
        sub_ko = sub_df.groupby('KO')[SAMPLE_COLS_18].sum()
        score = np.log10(sub_ko + 1).sum(axis=0)
        rho, pval = stats.spearmanr(score, arg)
        records.append({
            "Taxonomic_Level": level,
            "Taxon_Name": name,
            "Spearman_rho": round(rho, 3),
            "p_value": round(pval, 4)
        })

    res_df = pd.DataFrame(records)
    out_csv = os.path.join(os.path.dirname(__file__), "taxonomic_stratification_results.csv")
    res_df.to_csv(out_csv, index=False)
    print("=== Taxonomic Stratification Results (Table S15) ===")
    print(res_df.to_string(index=False))

if __name__ == '__main__':
    run_taxonomic_stratification()
