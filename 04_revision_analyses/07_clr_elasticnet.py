"""
07_reviewer4_major1_clr_elasticnet.py

Purpose
-------
Addresses Reviewer 4 Major Issue 1(c):
Evaluates out-of-sample performance of Elastic Net under Centered-Log-Ratio (CLR)
transformation to determine whether multivariate prediction depends solely on total
functional-gene abundance scale or reflects genuine multi-module coordination.

Methods
-------
1. Reconstructs CLR module scores from 903 target KOs across 18 wastewater metagenomes.
2. Fits ElasticNetCV across 15 CLR module scores to predict total ARG burden.
3. Evaluates under:
   - Nested Leave-One-Out Cross-Validation (LOOCV, 18 folds)
   - Leave-One-Site-Out Cross-Validation (LOSO-site, 9 folds)
4. Computes 2,000 bootstrap resamples for 95% confidence intervals on R2 and Spearman rho.

Outputs
-------
- Revesion/tables/reviewer4_major1_clr_elasticnet.csv
"""

import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut, KFold
from sklearn.preprocessing import StandardScaler

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Load raw KO abundance and module map
abund = pd.read_excel(os.path.join(BASE_DIR, '1.step  Target module selection and k number recovery', 'step1 figures table', 'Step1_Abundance_Matrix.xlsx'), sheet_name='Full_Abundance_Matrix')
abund = abund.rename(columns={'# Gene Family':'KO_ID'}).drop(columns=['Module'])
sample_cols = [c for c in abund.columns if c != 'KO_ID']
freeze = pd.read_csv(os.path.join(BASE_DIR, 'BG', 'Project22_903KO_Module_Map.csv'))[['KO_ID','Module_name']]
data = abund.merge(freeze, on='KO_ID', how='inner')
ko_ids = data['KO_ID'].tolist()
module_of = dict(zip(data['KO_ID'], data['Module_name']))
modules = sorted(set(module_of.values()))

# CLR transformation across the 903 target KOs
X_raw = data[sample_cols].astype(float).values
n_ko, n_samp = X_raw.shape
X_clr = np.zeros_like(X_raw)
for j in range(n_samp):
    col = X_raw[:, j]
    nz = col[col > 0]
    pseudo = 0.5 * nz.min() if len(nz) else 1e-6
    col2 = np.where(col > 0, col, pseudo)
    logcol = np.log(col2)
    X_clr[:, j] = logcol - logcol.mean()

ko_to_row = {k: i for i, k in enumerate(ko_ids)}
clr_scores = pd.DataFrame(index=sample_cols, columns=modules, dtype=float)
for m in modules:
    rows = [ko_to_row[k] for k in ko_ids if module_of[k] == m]
    clr_scores[m] = X_clr[rows, :].sum(axis=0)

target_path = os.path.join(BASE_DIR, '10. Final Manuscript (Network-Resistome Integration)', 'tables', 'ARG_total_burden_per_sample.csv')
y = pd.read_csv(target_path, index_col=0).iloc[:, 0]
common = clr_scores.index.intersection(y.index)
X_clr = clr_scores.loc[common]
y = y.loc[common]

sites = pd.Series([s[:-1] for s in common], index=common)
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]

def bootstrap_ci(y_true, y_pred, n_boot=2000, seed=2026):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    r2_list, rho_list = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, yp = y_true[idx], y_pred[idx]
        if np.std(yt) == 0 or np.std(yp) == 0:
            continue
        r2 = 1 - np.sum((yt - yp)**2) / np.sum((yt - yt.mean())**2)
        rho, _ = spearmanr(yt, yp)
        r2_list.append(r2)
        rho_list.append(rho)
    r2_ci = np.percentile(r2_list, [2.5, 97.5])
    rho_ci = np.percentile(rho_list, [2.5, 97.5])
    return r2_ci, rho_ci

def main():
    # 1. CLR LOOCV
    loo = LeaveOneOut()
    preds_loo = np.zeros(len(y))
    for tr_idx, te_idx in loo.split(X_clr):
        Xtr, Xte = X_clr.iloc[tr_idx].values, X_clr.iloc[te_idx].values
        scaler = StandardScaler().fit(Xtr)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=0)
        model = ElasticNetCV(l1_ratio=L1_RATIOS, cv=inner_cv, max_iter=20000, random_state=0)
        model.fit(scaler.transform(Xtr), y.iloc[tr_idx].values)
        preds_loo[te_idx] = model.predict(scaler.transform(Xte))

    # 2. CLR LOSO-site
    logo = LeaveOneGroupOut()
    preds_loso = np.zeros(len(y))
    for tr_idx, te_idx in logo.split(X_clr, y, groups=sites):
        Xtr, Xte = X_clr.iloc[tr_idx].values, X_clr.iloc[te_idx].values
        scaler = StandardScaler().fit(Xtr)
        inner_cv = KFold(n_splits=4, shuffle=True, random_state=0)
        model = ElasticNetCV(l1_ratio=L1_RATIOS, cv=inner_cv, max_iter=20000, random_state=0)
        model.fit(scaler.transform(Xtr), y.iloc[tr_idx].values)
        preds_loso[te_idx] = model.predict(scaler.transform(Xte))

    results = []
    for cv_name, preds in [('LOOCV (18 folds)', preds_loo), ('LOSO-site (9 folds)', preds_loso)]:
        r2 = 1 - np.sum((y.values - preds)**2) / np.sum((y.values - y.values.mean())**2)
        rmse = np.sqrt(np.mean((y.values - preds)**2))
        rho, p = spearmanr(y.values, preds)
        r2_ci, rho_ci = bootstrap_ci(y.values, preds)
        results.append({
            'Model': 'Elastic Net (15 modules, CLR transform)',
            'Cross-Validation Scheme': cv_name,
            'R2': round(float(r2), 3),
            'R2_95CI': f'[{r2_ci[0]:.3f}, {r2_ci[1]:.3f}]',
            'Spearman_rho': round(float(rho), 3),
            'p_value': f'{p:.4e}' if p < 0.001 else f'{p:.4f}',
            'rho_95CI': f'[{rho_ci[0]:.3f}, {rho_ci[1]:.3f}]',
            'RMSE': round(float(rmse), 3)
        })

    df = pd.DataFrame(results)
    out_csv = os.path.join(BASE_DIR, 'Revesion', 'tables', 'reviewer4_major1_clr_elasticnet.csv')
    df.to_csv(out_csv, index=False)
    print('Saved:', out_csv)
    print(df.to_string())

if __name__ == '__main__':
    main()

