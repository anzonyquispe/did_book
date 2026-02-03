"""
=============================================================================
CHAPTER 4, SECTION 1: Estimators with Controls
Moser & Voena (2012) - De Chaisemartin DID Textbook
=============================================================================

Stata equivalent:
    reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)
    test reltimeminus1 ... reltimeminus18
    reg patents1900 treatmentgroup if year==1900

Requirements:
    pip install pyfixest pandas numpy matplotlib statsmodels scipy
"""

import pandas as pd
import numpy as np
import pyfixest as pf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy import stats as sp_stats
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# LOAD DATA
# =============================================================================
DATA_PATH = 'moser_voena_didtextbook.dta'

df = pd.read_stata(DATA_PATH)
df['subclass'] = df['subclass'].astype(int)
df['year'] = df['year'].astype(int)
df['patents1900'] = df['patents1900'].astype(int)
df['treatmentgroup'] = df['treatmentgroup'].astype(int)

reltimeminus_vars = [f'reltimeminus{i}' for i in range(1, 19)]
reltimeplus_vars = [f'reltimeplus{i}' for i in range(1, 22)]
all_reltime = reltimeminus_vars + reltimeplus_vars

print("=" * 70)
print("CHAPTER 4, SECTION 1: TWFE Event-Study with Controls")
print("=" * 70)


# =============================================================================
# 1) TWFE with controls: year x patents1900 interaction FEs
# Stata: reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)
# =============================================================================
print("\n*** 1) TWFE Event-Study controlling for year x patents1900 ***")

fml = "patents ~ " + " + ".join(all_reltime) + " | treatmentgroup + year^patents1900"
model = pf.feols(fml, data=df, vcov={"CRV1": "subclass"})

coefs = model.coef()
ci = model.confint()

print(f"\n{'Variable':<20} {'Coef':>10} {'95% CI lo':>10} {'95% CI hi':>10}")
print("-" * 55)
for v in all_reltime:
    print(f"{v:<20} {coefs[v]:>10.6f} {ci.loc[v, '2.5%']:>10.6f} {ci.loc[v, '97.5%']:>10.6f}")


# =============================================================================
# Joint F-test for pre-trends
# Stata: test reltimeminus1 ... reltimeminus18
# =============================================================================
print("\n*** Joint F-test: all pre-treatment coefficients = 0 ***")

V = model._vcov
pre_idx = [list(coefs.index).index(v) for v in reltimeminus_vars]
pre_coefs = np.array([coefs[v] for v in reltimeminus_vars])
V_pre = V[np.ix_(pre_idx, pre_idx)]

wald = pre_coefs @ np.linalg.inv(V_pre) @ pre_coefs
k = len(pre_coefs)
f_stat = wald / k
p_value = 1 - sp_stats.chi2.cdf(wald, k)
print(f"  Wald chi2({k})  = {wald:.4f}")
print(f"  F({k},.)       = {f_stat:.4f}")
print(f"  p-value        = {p_value:.4f}")


# =============================================================================
# Event-study plot: graphES_moser_controls.pdf
# =============================================================================
print("\n*** Generating event-study plot ***")

# Build arrays: pre-treatment (reversed so -18 comes first), zero, post-treatment
time_vals = list(range(-18, 0)) + [0] + list(range(1, 22))
coef_vals = ([coefs[f'reltimeminus{i}'] for i in range(18, 0, -1)]
             + [0.0]
             + [coefs[f'reltimeplus{i}'] for i in range(1, 22)])
ci_lo_vals = ([ci.loc[f'reltimeminus{i}', '2.5%'] for i in range(18, 0, -1)]
              + [0.0]
              + [ci.loc[f'reltimeplus{i}', '2.5%'] for i in range(1, 22)])
ci_hi_vals = ([ci.loc[f'reltimeminus{i}', '97.5%'] for i in range(18, 0, -1)]
              + [0.0]
              + [ci.loc[f'reltimeplus{i}', '97.5%'] for i in range(1, 22)])

t = np.array(time_vals)
b = np.array(coef_vals)
lo = np.array(ci_lo_vals)
hi = np.array(ci_hi_vals)

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(t, b, 'o-', color='navy', markersize=5, linewidth=1.2, zorder=3)
ax.vlines(t, lo, hi, color='maroon', linewidth=1, zorder=2)
ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax.set_title('TWFE Event-study estimates (with controls)')
ax.set_xlabel('Relative time to year before TWEA')
ax.set_ylabel('Effect')
ax.set_xticks(range(-18, 22, 3))
ax.set_ylim(-0.25, 1.0)
ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
plt.tight_layout()
fig.savefig('graphES_moser_controls.pdf', bbox_inches='tight')
fig.savefig('graphES_moser_controls.png', dpi=150, bbox_inches='tight')
print("  Saved: graphES_moser_controls.pdf / .png")
plt.close()


# =============================================================================
# Testing that treatment group and covariate are correlated
# Stata: reg patents1900 treatmentgroup if year==1900
# =============================================================================
print("\n*** Testing correlation: patents1900 ~ treatmentgroup (year==1900) ***")

df_1900 = df[df['year'] == 1900].copy()
X = sm.add_constant(df_1900['treatmentgroup'])
y = df_1900['patents1900']
ols_corr = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': df_1900['subclass']})
print(f"  treatmentgroup coef = {ols_corr.params['treatmentgroup']:.6f}")
print(f"  t-stat              = {ols_corr.tvalues['treatmentgroup']:.4f}")
print(f"  p-value             = {ols_corr.pvalues['treatmentgroup']:.4f}")
print(f"  => {'Correlated' if ols_corr.pvalues['treatmentgroup'] < 0.05 else 'Not correlated'} at 5%")


# =============================================================================
# Reference TWFE (post-only, no pre-period dummies) — needed for later sections
# Stata: reg patents i.yearpost treatmentgroup reltimeplus*, cluster(subclass)
# =============================================================================
print("\n*** Reference TWFE post-only (for comparison in Sections 2-5) ***")

fml_ref = "patents ~ " + " + ".join(reltimeplus_vars) + " + treatmentgroup | year"
model_ref = pf.feols(fml_ref, data=df, vcov={"CRV1": "subclass"})
ref_coefs = model_ref.coef()
ref_ci = model_ref.confint()

for i in [1, 7, 14, 21]:
    v = f'reltimeplus{i}'
    print(f"  {v:<15} = {ref_coefs[v]:>9.6f}  [{ref_ci.loc[v, '2.5%']:.4f}, {ref_ci.loc[v, '97.5%']:.4f}]")

# Save for other sections to load
ref_data = {
    'time': np.arange(1, 22),
    'coef': np.array([ref_coefs[f'reltimeplus{i}'] for i in range(1, 22)]),
    'ci_lo': np.array([ref_ci.loc[f'reltimeplus{i}', '2.5%'] for i in range(1, 22)]),
    'ci_hi': np.array([ref_ci.loc[f'reltimeplus{i}', '97.5%'] for i in range(1, 22)]),
}
np.savez('twfe_reference.npz', **ref_data)
print("  Saved: twfe_reference.npz")


print("\n" + "=" * 70)
print("SECTION 1 COMPLETE")
print("=" * 70)
