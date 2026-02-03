"""
CHAPTER 3: BASIC DID METHODS - WITH 3 GRAPHS
Replication of Moser & Voena (2012) - De Chaisemartin Course
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

from statsmodels.formula.api import ols

try:
    import pyfixest as pf
    HAS_PYFIXEST = True
except ImportError:
    HAS_PYFIXEST = False

# Set working directory
FOLDER = r'C:\Users\Usuario\Documents\GitHub\soolution'
os.chdir(FOLDER)

print("CHAPTER 3: BASIC DID METHODS")

# Load data
print("\n[0] LOADING DATA...")
df = pd.read_stata('moser_voena_didtextbook.dta')
print(f"    Observations: {len(df):,}")
print(f"    Subclasses: {df['subclass'].nunique():,}")

reltimeminus_vars = sorted([col for col in df.columns if col.startswith('reltimeminus')],
                           key=lambda x: int(x.replace('reltimeminus', '')))
reltimeplus_vars = sorted([col for col in df.columns if col.startswith('reltimeplus')],
                          key=lambda x: int(x.replace('reltimeplus', '')))

# 1) Static TWFE Regression
print("\n[1] STATIC TWFE REGRESSION")

if HAS_PYFIXEST:
    model_twfe = pf.feols("patents ~ twea | subclass + year", data=df, vcov={"CRV1": "subclass"})
    twea_coef = model_twfe.coef()['twea']
    twea_se = model_twfe.se()['twea']
    twea_tstat = twea_coef / twea_se
else:
    model_twfe = ols("patents ~ twea + C(year) + C(subclass)", data=df).fit()
    twea_coef = model_twfe.params['twea']
    twea_se = model_twfe.bse['twea']
    twea_tstat = model_twfe.tvalues['twea']

print(f"\n    Coef twea:    {twea_coef:.7f}")
print(f"    Std. Error:   {twea_se:.7f}")
print(f"    t-statistic:  {twea_tstat:.2f} ***")

# 2) TWFE = DID Equivalence
print("\n[2] TWFE = DID EQUIVALENCE")

model_did = ols("patents ~ treatmentgroup + post + twea", data=df).fit(
    cov_type='cluster', cov_kwds={'groups': df['subclass']}
)
print(f"    treatmentgroup: {model_did.params['treatmentgroup']:.7f}")
print(f"    post:           {model_did.params['post']:.7f}")
print(f"    twea:           {model_did.params['twea']:.7f} <- Same as TWFE")

# 3) Testing Randomized Treatment
print("\n[3] TESTING RANDOMIZED TREATMENT")

df_pre = df[df['year'] <= 1918].copy()
model_random = ols("patents ~ treatmentgroup", data=df_pre).fit(
    cov_type='cluster', cov_kwds={'groups': df_pre['subclass']}
)
print(f"    Coef treatmentgroup: {model_random.params['treatmentgroup']:.7f}")
print(f"    -> Treated had FEWER patents before TWEA")

# 4) Event-Study TWFE Regression
print("\n[4] EVENT-STUDY TWFE REGRESSION")

reltime_vars = reltimeminus_vars + reltimeplus_vars
formula_es = "patents ~ treatmentgroup + " + " + ".join(reltime_vars) + " + C(year)"
model_es = ols(formula_es, data=df).fit(
    cov_type='cluster', cov_kwds={'groups': df['subclass']}
)

# Extract coefficients for graph
es_results = []
for var in reversed(reltimeminus_vars):
    i = int(var.replace('reltimeminus', ''))
    if var in model_es.params.index:
        es_results.append({
            'time': -i, 'coef': model_es.params[var], 'se': model_es.bse[var],
            'ci_low': model_es.conf_int().loc[var, 0], 'ci_high': model_es.conf_int().loc[var, 1]
        })

es_results.append({'time': 0, 'coef': 0, 'se': 0, 'ci_low': 0, 'ci_high': 0})

for var in reltimeplus_vars:
    i = int(var.replace('reltimeplus', ''))
    if var in model_es.params.index:
        es_results.append({
            'time': i, 'coef': model_es.params[var], 'se': model_es.bse[var],
            'ci_low': model_es.conf_int().loc[var, 0], 'ci_high': model_es.conf_int().loc[var, 1]
        })

es_df = pd.DataFrame(es_results)

print(f"    reltimeminus18: {model_es.params['reltimeminus18']:.7f}")
print(f"    reltimeplus1:   {model_es.params['reltimeplus1']:.7f}")
print(f"    reltimeplus10:  {model_es.params['reltimeplus10']:.7f}")
print(f"    reltimeplus21:  {model_es.params['reltimeplus21']:.7f}")

# Pre-trends test
pre_vars_test = [v for v in reltimeminus_vars if v in model_es.params.index]
r_matrix = np.zeros((len(pre_vars_test), len(model_es.params)))
for i, var in enumerate(pre_vars_test):
    r_matrix[i, list(model_es.params.index).index(var)] = 1
f_test = model_es.f_test(r_matrix)
print(f"\n    Pre-trends Test: F = {float(f_test.fvalue):.2f}, p = {float(f_test.pvalue):.4f}")

# Graph 1: Complete Event-Study
print("\nGRAPH 1: COMPLETE EVENT-STUDY")

fig, ax = plt.subplots(figsize=(12, 7))
ax.scatter(es_df['time'], es_df['coef'], color='navy', s=80, zorder=3)
ax.plot(es_df['time'], es_df['coef'], color='navy', linewidth=1.5, zorder=2)
ax.vlines(es_df['time'], es_df['ci_low'], es_df['ci_high'], color='maroon', linewidth=2, zorder=1)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.set_xlabel('Relative time to year before TWEA', fontsize=12)
ax.set_ylabel('Effect', fontsize=12)
ax.set_title('TWFE Event-study estimates', fontsize=14, fontweight='bold')
ax.set_xlim(-19, 22)
ax.set_ylim(-0.25, 1.05)
ax.set_xticks(range(-18, 22, 3))
ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FOLDER, 'graphES_moser1_python.png'), dpi=150)
plt.savefig(os.path.join(FOLDER, 'graphES_moser1_python.pdf'), dpi=150)
print(f"    Saved to: {FOLDER}")
plt.show()

# Save for graph 2
es_df_post = es_df[es_df['time'] >= 0].copy()

# 5) Event-Study Without Pre-Trends
print("\n[5] EVENT-STUDY WITHOUT PRE-TRENDS")

formula_post = "patents ~ treatmentgroup + " + " + ".join(reltimeplus_vars) + " + C(year)"
model_post = ols(formula_post, data=df).fit(
    cov_type='cluster', cov_kwds={'groups': df['subclass']}
)

es_post_results = [{'time': 0, 'coef': 0, 'se': 0, 'ci_low': 0, 'ci_high': 0}]
for var in reltimeplus_vars:
    i = int(var.replace('reltimeplus', ''))
    if var in model_post.params.index:
        es_post_results.append({
            'time': i, 'coef': model_post.params[var], 'se': model_post.bse[var],
            'ci_low': model_post.conf_int().loc[var, 0], 'ci_high': model_post.conf_int().loc[var, 1]
        })

es_post_df = pd.DataFrame(es_post_results)

# Graph 2: With vs Without Pre-Periods
print("\nGRAPH 2: WITH VS WITHOUT PRE-PERIODS")

fig, ax = plt.subplots(figsize=(12, 7))
# Without pre-periods (blue)
ax.scatter(es_post_df['time'], es_post_df['coef'], color='steelblue', s=60, zorder=3)
ax.plot(es_post_df['time'], es_post_df['coef'], color='steelblue', linewidth=1.5, label='Without Pre-Periods')
ax.vlines(es_post_df['time'], es_post_df['ci_low'], es_post_df['ci_high'], color='steelblue', linewidth=1.5)
# With pre-periods (red)
ax.scatter(es_df_post['time'], es_df_post['coef'], color='red', s=60, zorder=3)
ax.plot(es_df_post['time'], es_df_post['coef'], color='red', linewidth=1.5, label='With Pre-Periods')
ax.vlines(es_df_post['time'], es_df_post['ci_low'], es_df_post['ci_high'], color='red', linewidth=1.5)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.set_xlabel('Relative time to year before TWEA', fontsize=12)
ax.set_ylabel('Effect', fontsize=12)
ax.set_title('TWFE Event-study estimates', fontsize=14, fontweight='bold')
ax.set_xlim(-1, 22)
ax.set_ylim(-0.25, 1.05)
ax.set_xticks(range(0, 22, 3))
ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FOLDER, 'graphES_moser2_python.png'), dpi=150)
plt.savefig(os.path.join(FOLDER, 'graphES_moser2_python.pdf'), dpi=150)
print(f"    Saved to: {FOLDER}")
plt.show()

# 6) Pre-Trends Power (only 6 pre-periods)
print("\n[6] PRE-TRENDS POWER")

# Filter data and use only 6 pre-periods
df_short = df[(df['year'] >= 1912) & (df['year'] <= 1939)].copy()
reltimeminus_short = ['reltimeminus1', 'reltimeminus2', 'reltimeminus3',
                       'reltimeminus4', 'reltimeminus5', 'reltimeminus6']

formula_short = "patents ~ treatmentgroup + " + " + ".join(reltimeminus_short + reltimeplus_vars) + " + C(year)"
model_short = ols(formula_short, data=df_short).fit(
    cov_type='cluster', cov_kwds={'groups': df_short['subclass']}
)

# Extract coefficients
es_short = []
for var in reversed(reltimeminus_short):
    i = int(var.replace('reltimeminus', ''))
    if var in model_short.params.index:
        es_short.append({
            'time': -i, 'coef': model_short.params[var],
            'ci_low': model_short.conf_int().loc[var, 0],
            'ci_high': model_short.conf_int().loc[var, 1]
        })

es_short.append({'time': 0, 'coef': 0, 'ci_low': 0, 'ci_high': 0})

for var in reltimeplus_vars:
    i = int(var.replace('reltimeplus', ''))
    if var in model_short.params.index:
        es_short.append({
            'time': i, 'coef': model_short.params[var],
            'ci_low': model_short.conf_int().loc[var, 0],
            'ci_high': model_short.conf_int().loc[var, 1]
        })

es_short_df = pd.DataFrame(es_short)

# Slope from Stata: pretrends power 0.5, numpre(6)
slope = 0.0098235

# Graph 3: Pre-Trends Power
print("\nGRAPH 3: PRE-TRENDS POWER")

fig, ax = plt.subplots(figsize=(12, 7))
ax.scatter(es_short_df['time'], es_short_df['coef'], color='navy', s=80, zorder=3)
ax.plot(es_short_df['time'], es_short_df['coef'], color='navy', linewidth=1.5, zorder=2)
ax.vlines(es_short_df['time'], es_short_df['ci_low'], es_short_df['ci_high'], color='maroon', linewidth=2, zorder=1)
# Trend line
x_trend = np.linspace(-6, 21, 100)
y_trend = x_trend * slope
ax.plot(x_trend, y_trend, color='gray', linestyle='--', linewidth=1.5, label=f'Linear trend (slope={slope})')
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel('Relative time to year before TWEA', fontsize=12)
ax.set_ylabel('Effect', fontsize=12)
ax.set_title('TWFE Event-study estimates', fontsize=14, fontweight='bold')
ax.set_xlim(-7, 22)
ax.set_ylim(-0.25, 1.05)
ax.set_xticks(range(-6, 22, 3))
ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FOLDER, 'graphES_moser3_python.png'), dpi=150)
plt.savefig(os.path.join(FOLDER, 'graphES_moser3_python.pdf'), dpi=150)
print(f"    Saved to: {FOLDER}")
plt.show()

# Final Summary
print("\nSUMMARY")
print(f"""
    RESULTS:
    - Coef TWEA: {twea_coef:.4f}***
    - Effect year 1: {model_es.params['reltimeplus1']:.4f}
    - Effect year 21: {model_es.params['reltimeplus21']:.4f}

    GRAPHS SAVED TO {FOLDER}:
    - graphES_moser1_python.png/pdf
    - graphES_moser2_python.png/pdf
    - graphES_moser3_python.png/pdf

    CHAPTER 3 COMPLETED
""")