!pip install pandas numpy matplotlib scipy statsmodels pyfixest cvxpy
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats, optimize, linalg
import warnings
import os
import time

warnings.filterwarnings("ignore")

# Optional imports
try:
    import pyfixest as pf
    HAS_PYFIXEST = True
except ImportError:
    HAS_PYFIXEST = False
    print("⚠  pyfixest not installed. pip install pyfixest")

try:
    import statsmodels.api as sm
    from statsmodels.formula.api import ols as sm_ols
    HAS_SM = True
except ImportError:
    HAS_SM = False
    print("⚠  statsmodels not installed. pip install statsmodels")

try:
    import cvxpy as cp
    HAS_CVXPY = True
except ImportError:
    HAS_CVXPY = False
    print("⚠  cvxpy not installed (needed for SC/SDID). pip install cvxpy")
    
#CHARGE DATA
path = Path('/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta')
df = pd.read_stata(path)
# CONFIG
DATA_FILE = "moser_voena_didtextbook.dta"
GRAPH_DIR = "graphs"
np.random.seed(1)

os.makedirs(GRAPH_DIR, exist_ok=True)
# HELPERS
def get_reltime_vars(df):
    minus_v = sorted(
        [c for c in df.columns if c.startswith("reltimeminus")],
        key=lambda x: int(x.replace("reltimeminus", "")),
    )
    plus_v = sorted(
        [c for c in df.columns if c.startswith("reltimeplus")],
        key=lambda x: int(x.replace("reltimeplus", "")),
    )
    return minus_v, plus_v


def cluster_se(X, y, resid, clusters):
    """OLS cluster-robust SEs (Stata-style)."""
    n, k = X.shape
    unique_cl = np.unique(clusters)
    G = len(unique_cl)
    bread = np.linalg.inv(X.T @ X)
    meat = np.zeros((k, k))
    for cl in unique_cl:
        idx = clusters == cl
        ei = resid[idx]
        Xi = X[idx]
        score = Xi.T @ ei
        meat += np.outer(score, score)
    scale = G / (G - 1) * (n - 1) / (n - k)
    V = bread @ meat @ bread * scale
    return np.sqrt(np.diag(V))


def ols_with_cluster(df_reg, y_var, x_vars, fe_vars, cluster_var):
    """Manual OLS with FE dummies and cluster-robust SEs."""
    df_r = df_reg.dropna(subset=[y_var] + x_vars + fe_vars + [cluster_var]).copy()
    y = df_r[y_var].values.astype(float)
    
    # Build X: x_vars + FE dummies (dropping first)
    X_parts = [df_r[x_vars].values.astype(float)]
    for fe in fe_vars:
        dummies = pd.get_dummies(df_r[fe], drop_first=True, dtype=float)
        X_parts.append(dummies.values)
    X = np.column_stack(X_parts)
    X = np.column_stack([np.ones(len(y)), X])  # intercept
    
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ beta
    se = cluster_se(X, y, resid, df_r[cluster_var].values)
    
    # Results for x_vars only (skip intercept at pos 0)
    coefs = beta[1 : 1 + len(x_vars)]
    ses = se[1 : 1 + len(x_vars)]
    return coefs, ses


def stata_style_plot():
    """Apply Stata s2color-like styling."""
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.grid": False,
        "axes.linewidth": 0.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "font.family": "sans-serif",
        "font.size": 10,
    })
minus_vars, plus_vars = get_reltime_vars(df)
all_reltime = minus_vars + plus_vars

#CHAPTER 3
# 1) Static TWFE regression
if HAS_PYFIXEST:
    m1 = pf.feols("patents ~ twea | subclass + year", data=df, vcov={"CRV1": "subclass"})
    print(m1.summary())
else:
    coefs, ses = ols_with_cluster(df, "patents", ["twea"], ["subclass", "year"], "subclass")
    t = coefs[0] / ses[0]
    p = 2 * (1 - stats.t.cdf(abs(t), df["subclass"].nunique() - 1))
    print(f"    twea: coef = {coefs[0]:.6f}, se = {ses[0]:.6f}, t = {t:.2f}, p = {p:.4f}")

# reghdfe equivalence
print("\n    reghdfe patents twea, absorb(subclass year) cluster(subclass)")
if HAS_PYFIXEST:
    m1b = pf.feols("patents ~ twea | subclass + year", data=df, vcov={"CRV1": "subclass"})
    print(f"    twea: {m1b.coef()['twea']:.6f}  (same as above)")

# 2) Equivalence between static TWFE regression and DID
print("\n*** 2) Equivalence between static TWFE and DID")
print("    reg patents treatmentgroup post twea, cluster(subclass)")

if HAS_PYFIXEST:
    m2 = pf.feols("patents ~ treatmentgroup + post + twea", data=df, vcov={"CRV1": "subclass"})
    print(m2.summary())
else:
    coefs, ses = ols_with_cluster(df, "patents", ["treatmentgroup", "post", "twea"], [], "subclass")
    for i, v in enumerate(["treatmentgroup", "post", "twea"]):
        print(f"    {v}: coef = {coefs[i]:.6f}, se = {ses[i]:.6f}")

# 3) Testing randomized treatment
print("\n*** 3) Testing randomized treatment")
print("    reg patents treatmentgroup if year<=1918, cluster(subclass)")

df_pre = df[df["year"] <= 1918].copy()
if HAS_PYFIXEST:
    m3 = pf.feols("patents ~ treatmentgroup", data=df_pre, vcov={"CRV1": "subclass"})
    print(m3.summary())
else:
    coefs, ses = ols_with_cluster(df_pre, "patents", ["treatmentgroup"], [], "subclass")
    print(f"    treatmentgroup: coef = {coefs[0]:.6f}, se = {ses[0]:.6f}")
    
# 4) Event-study TWFE regression
print("\n*** 4) Event-study TWFE regression")
print("    reg patents i.year treatmentgroup reltimeminus* reltimeplus*, cluster(subclass)")
if HAS_PYFIXEST:
    fml_es = "patents ~ " + " + ".join(all_reltime) + " + treatmentgroup | year"
    m4 = pf.feols(fml_es, data=df, vcov={"CRV1": "subclass"})
    print(m4.summary())
else:
    coefs, ses = ols_with_cluster(
        df, "patents", all_reltime + ["treatmentgroup"], ["year"], "subclass"
    )
    for i, v in enumerate(all_reltime):
        print(f"    {v}: {coefs[i]:.6f} ({ses[i]:.6f})")

# F-test sobre coeficientes pre-tendencia
print("\n    F-test on all reltimeminus coefficients:")
if HAS_PYFIXEST:
    c = m4.coef()
    v = pd.DataFrame(m4._vcov, index=c.index, columns=c.index)
    pre_names = [v_ for v_ in minus_vars if v_ in c.index]
    beta_pre = c[pre_names].values
    V_pre = v.loc[pre_names, pre_names].values
    F_stat = beta_pre @ np.linalg.inv(V_pre) @ beta_pre / len(pre_names)
    p_val = 1 - stats.f.cdf(F_stat, len(pre_names), df["subclass"].nunique() - 1)
    print(f"    F({len(pre_names)}, {df['subclass'].nunique()-1}) = {F_stat:.4f}")
    print(f"    Prob > F = {p_val:.4f}")
#  Event-study plot (Graph 1: graphES_moser1) 
print("\n    Generating graphES_moser1...")
stata_style_plot()

if HAS_PYFIXEST:
    fml_plot = "patents ~ " + " + ".join(all_reltime) + " + treatmentgroup | year"
    m4p = pf.feols(fml_plot, data=df, vcov={"CRV1": "subclass"})
    coefs_d = m4p.coef()
    se_d = m4p.se()

    reltime = [0]
    est = [0.0]
    ci_lo = [0.0]
    ci_hi = [0.0]

    for v in minus_vars:
        k = -int(v.replace("reltimeminus", ""))
        if v in coefs_d.index:
            reltime.append(k)
            est.append(coefs_d[v])
            ci_lo.append(coefs_d[v] - 1.96 * se_d[v])
            ci_hi.append(coefs_d[v] + 1.96 * se_d[v])
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_d.index:
            reltime.append(k)
            est.append(coefs_d[v])
            ci_lo.append(coefs_d[v] - 1.96 * se_d[v])
            ci_hi.append(coefs_d[v] + 1.96 * se_d[v])

    order = np.argsort(reltime)
    reltime = np.array(reltime)[order]
    est = np.array(est)[order]
    ci_lo = np.array(ci_lo)[order]
    ci_hi = np.array(ci_hi)[order]

    # Store post-period results for later comparison (Chapter 4)
    mask_post = reltime >= 0
    res_post = {
        "reltime": reltime[mask_post],
        "est": est[mask_post],
        "ci_lo": ci_lo[mask_post],
        "ci_hi": ci_hi[mask_post],
    }
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(reltime, est, "o-", color="navy", markersize=5, linewidth=1)
    ax.vlines(reltime, ci_lo, ci_hi, color="maroon", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="-")
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(-18, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser1.pdf"))
    print("    ✓ Saved graphs/graphES_moser1.pdf")
    plt.show()

#5) Event-study without pre-trends estimates
print("\n*** 5) Event-study TWFE without pre-trends estimates")
print("    reg patents i.yearpost treatmentgroup reltimeplus*, cluster(subclass)")

if HAS_PYFIXEST:
    fml_nopretrend = "patents ~ " + " + ".join(plus_vars) + " + treatmentgroup | yearpost"
    m5 = pf.feols(fml_nopretrend, data=df, vcov={"CRV1": "subclass"})
    print(m5.summary())

    # Build results for "without pre-periods"
    coefs_5 = m5.coef()
    se_5 = m5.se()
    rt_nopre = [0]
    est_nopre = [0.0]
    ci_lo_nopre = [0.0]
    ci_hi_nopre = [0.0]
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_5.index:
            rt_nopre.append(k)
            est_nopre.append(coefs_5[v])
            ci_lo_nopre.append(coefs_5[v] - 1.96 * se_5[v])
            ci_hi_nopre.append(coefs_5[v] + 1.96 * se_5[v])

    order2 = np.argsort(rt_nopre)
    rt_nopre = np.array(rt_nopre)[order2]
    est_nopre = np.array(est_nopre)[order2]
    ci_lo_nopre = np.array(ci_lo_nopre)[order2]
    ci_hi_nopre = np.array(ci_hi_nopre)[order2]

    # Graph 2: graphES_moser2 
    print("\n    Generating graphES_moser2...")
    fig, ax = plt.subplots(figsize=(8, 5))

    # Without pre-periods (blue)
    ax.plot(rt_nopre, est_nopre, "o-", color="steelblue", markersize=4, linewidth=1,
            label="Without Pre-Periods")
    ax.vlines(rt_nopre, ci_lo_nopre, ci_hi_nopre, color="steelblue", linewidth=1)

    # With pre-periods (red) — from res_post
    ax.plot(res_post["reltime"], res_post["est"], "o-", color="red", markersize=4,
            linewidth=1, label="With Pre-Periods")
    ax.vlines(res_post["reltime"], res_post["ci_lo"], res_post["ci_hi"],
              color="red", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(0, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    ax.legend(loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser2.pdf"))
    
    print("    ✓ Saved graphs/graphES_moser2.pdf")
    plt.show()

    # Numerically equivalent to Borusyak et al / having all year FEs
    print("\n    Note: Numerically equivalent to Borusyak et al (did_imputation)")
    print("    Also equivalent to: reg patents i.year treatmentgroup reltimeplus*")

# 6) Linear pre-trends we could fail to detect (pretrends power)
print("\n*** 6) Linear pre-trends we could fail to detect")

if HAS_PYFIXEST:
    # Restrict to 6 pre-periods (year >= 1912 & year <= 1939)
    df_sub = df[(df["year"] >= 1912) & (df["year"] <= 1939)].copy()
    short_minus = [f"reltimeminus{i}" for i in range(1, 7)]
    short_minus_present = [v for v in short_minus if v in df_sub.columns]

    # reghdfe patents reltimeminus* reltimeplus*, absorb(treatmentgroup year) cluster(subclass)
    fml_pre = "patents ~ " + " + ".join(short_minus_present + plus_vars) + " | treatmentgroup + year"
    m6_full = pf.feols(fml_pre, data=df, vcov={"CRV1": "subclass"})

    # Pretrends power calculation (Roth 2022)
    coefs_pre = m6_full.coef()
    vcov_pre = pd.DataFrame(m6_full._vcov, index=coefs_pre.index, columns=coefs_pre.index)

    # Extract pre-trend coefficients (first 6)
    pre_names_6 = [v for v in short_minus_present if v in coefs_pre.index][:6]
    beta_pre6 = coefs_pre[pre_names_6].values
    V_pre6 = vcov_pre.loc[pre_names_6, pre_names_6].values

    # Pre-trends power: slope that gives 50% power for the joint F-test
    n_pre = len(pre_names_6)
    dfn = n_pre
    dfd = df["subclass"].nunique() - 1

    def power_at_slope(slope):
        """Power of F-test under linear trend with given slope."""
        mu = np.array([slope * (-int(v.replace("reltimeminus", ""))) for v in pre_names_6])
        ncp = mu @ np.linalg.inv(V_pre6) @ mu
        cv = stats.f.ppf(0.95, dfn, dfd)
        power = 1 - stats.ncf.cdf(cv, dfn, dfd, ncp)
        return power

    # Search for slope giving 50% power
    try:
        result = optimize.brentq(lambda s: power_at_slope(s) - 0.5, 0.0001, 1.0)
        slope_50 = result
        print(f"    Slope for 50% power: {slope_50:.6f}")
    except Exception:
        slope_50 = None
        print("    Could not find slope for 50% power (search failed)")
 # Event-study on restricted sample for the graph
    fml_sub = "patents ~ " + " + ".join(short_minus_present + plus_vars) + " + treatmentgroup | year"
    m6_sub = pf.feols(fml_sub, data=df_sub, vcov={"CRV1": "subclass"})
    coefs_6 = m6_sub.coef()
    se_6 = m6_sub.se()

    rt6 = [0]
    est6 = [0.0]
    ci_lo6 = [0.0]
    ci_hi6 = [0.0]

    for v in short_minus_present:
        k = -int(v.replace("reltimeminus", ""))
        if v in coefs_6.index:
            rt6.append(k)
            est6.append(coefs_6[v])
            ci_lo6.append(coefs_6[v] - 1.96 * se_6[v])
            ci_hi6.append(coefs_6[v] + 1.96 * se_6[v])
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_6.index:
            rt6.append(k)
            est6.append(coefs_6[v])
            ci_lo6.append(coefs_6[v] - 1.96 * se_6[v])
            ci_hi6.append(coefs_6[v] + 1.96 * se_6[v])

    order6 = np.argsort(rt6)
    rt6 = np.array(rt6)[order6]
    est6 = np.array(est6)[order6]
    ci_lo6 = np.array(ci_lo6)[order6]
    ci_hi6 = np.array(ci_hi6)[order6]

    #  Graph 3: graphES_moser3 
    print("    Generating graphES_moser3...")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(rt6, est6, "o-", color="navy", markersize=5, linewidth=1)
    ax.vlines(rt6, ci_lo6, ci_hi6, color="maroon", linewidth=1)
    if slope_50 is not None:
        x_line = np.linspace(-6, 21, 100)
        ax.plot(x_line, x_line * slope_50, "--", color="gray", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(-6, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser3.pdf"))
    print("    ✓ Saved graphs/graphES_moser3.pdf")
    plt.show()

# 7) Variance of the effect of having been exposed to treatment for 14 years
df_1932 = df[df["year"] == 1932].copy()
g0 = df_1932.loc[df_1932["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
g1 = df_1932.loc[df_1932["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()

# ─── 1. sdtest diffpatentswrt1918 if year==1932, by(treatmentgroup) ───
n0, n1 = len(g0), len(g1)
mean0, mean1 = g0.mean(), g1.mean()
sd0, sd1 = g0.std(ddof=1), g1.std(ddof=1)
se0, se1 = sd0 / np.sqrt(n0), sd1 / np.sqrt(n1)
ci_lo0, ci_hi0 = mean0 - 1.96 * se0, mean0 + 1.96 * se0
ci_lo1, ci_hi1 = mean1 - 1.96 * se1, mean1 + 1.96 * se1

combined = pd.concat([g0, g1])
n_c = len(combined)
mean_c = combined.mean()
sd_c = combined.std(ddof=1)
se_c = sd_c / np.sqrt(n_c)
ci_loc, ci_hic = mean_c - 1.96 * se_c, mean_c + 1.96 * se_c

print("    Variance ratio test")
print("    " + "-" * 70)
print(f"       Group |     Obs        Mean    Std. err.   Std. dev.   [95% conf. interval]")
print("    " + "-" * 70)
print(f"           0 | {n0:>7,}    {mean0:.7f}    {se0:.7f}    {sd0:.6f}    {ci_lo0:.7f}    {ci_hi0:.7f}")
print(f"           1 | {n1:>7,}    {mean1:.7f}    {se1:.7f}    {sd1:.6f}    {ci_lo1:.7f}    {ci_hi1:.7f}")
print("    " + "-" * 70)
print(f"    Combined | {n_c:>7,}    {mean_c:.7f}    {se_c:.7f}    {sd_c:.6f}    {ci_loc:.7f}    {ci_hic:.7f}")
print("    " + "-" * 70)

F_var = sd0**2 / sd1**2
df1_f, df2_f = n0 - 1, n1 - 1
p_left = stats.f.cdf(F_var, df1_f, df2_f)
p_right = 1 - p_left
p_two = 2 * min(p_left, p_right)
print(f"        ratio = sd(0) / sd(1)                          f = {F_var:.4f}")
print(f"    H0: ratio = 1                         Degrees of freedom = {df1_f}, {df2_f}")
print(f"      Pr(F < f) = {p_left:.4f}         2*Pr(F < f) = {p_two:.4f}           Pr(F > f) = {p_right:.4f}")

#  2. di r(sd_2)-r(sd_1) 
sd_diff = sd1 - sd0
print(f"\n    r(sd_2) - r(sd_1) = {sd_diff:.8f}")

#  3. scalar sd_effects=r(sd_2)-r(sd_1) 
sd_effects = sd_diff
print(f"    scalar sd_effects = {sd_effects:.8f}")

#  4. reg diffpatentswrt1918 treatmentgroup if year==1932 
if HAS_PYFIXEST:
    m7 = pf.feols("diffpatentswrt1918 ~ treatmentgroup", data=df_1932)
    print(m7.summary())
    beta_did = m7.coef()["treatmentgroup"]
else:
    import statsmodels.api as sm
    X = sm.add_constant(df_1932["treatmentgroup"])
    m7 = sm.OLS(df_1932["diffpatentswrt1918"], X).fit()
    print(m7.summary())
    beta_did = m7.params["treatmentgroup"]

#  di _b[treatmentgroup]-1.96*sd_effects, _b[treatmentgroup]+1.96*sd_effects 
ci_lo_sd = beta_did - 1.96 * sd_effects
ci_hi_sd = beta_did + 1.96 * sd_effects
print(f"\n    _b[treatmentgroup] = {beta_did:.7f}")
print(f"    CI using SD of effects: [{ci_lo_sd:.8f}, {ci_hi_sd:.8f}]")

# 8) Placebo test of variance assumptions
print("\n*** 8) Placebo test (year=1904)")

#  1. sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup) 
df_1904 = df[df["year"] == 1904].copy()
g0p = df_1904.loc[df_1904["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
g1p = df_1904.loc[df_1904["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()

def print_sdtest(g0, g1, label=""):
    n0, n1 = len(g0), len(g1)
    mean0, mean1 = g0.mean(), g1.mean()
    sd0, sd1 = g0.std(ddof=1), g1.std(ddof=1)
    se0, se1 = sd0 / np.sqrt(n0), sd1 / np.sqrt(n1)
    
    combined = pd.concat([g0, g1])
    n_c = len(combined)
    mean_c = combined.mean()
    sd_c = combined.std(ddof=1)
    se_c = sd_c / np.sqrt(n_c)
    
    # Intervalos de confianza
    from scipy import stats as st
    ci_lo0, ci_hi0 = st.t.interval(0.95, n0 - 1, loc=mean0, scale=se0)
    ci_lo1, ci_hi1 = st.t.interval(0.95, n1 - 1, loc=mean1, scale=se1)
    ci_loc, ci_hic = st.t.interval(0.95, n_c - 1, loc=mean_c, scale=se_c)
    
    print(f"    {label}")
    print("    Variance ratio test")
    print("    " + "-" * 75)
    print(f"       Group |     Obs        Mean    Std. err.   Std. dev.   [95% conf. interval]")
    print("    " + "-" * 75)
    print(f"           0 | {n0:>7,}   {mean0:>10.7f}   {se0:>10.7f}   {sd0:>10.6f}   {ci_lo0:>11.7f}   {ci_hi0:>11.7f}")
    print(f"           1 | {n1:>7,}   {mean1:>10.7f}   {se1:>10.7f}   {sd1:>10.6f}   {ci_lo1:>11.7f}   {ci_hi1:>11.7f}")
    print("    " + "-" * 75)
    print(f"    Combined | {n_c:>7,}   {mean_c:>10.7f}   {se_c:>10.7f}   {sd_c:>10.6f}   {ci_loc:>11.7f}   {ci_hic:>11.7f}")
    print("    " + "-" * 75)
    
    if sd0 > 0 and sd1 > 0:
        F_var = sd0**2 / sd1**2
        df1_f, df2_f = n0 - 1, n1 - 1
        p_left = st.f.cdf(F_var, df1_f, df2_f)
        p_right = 1 - p_left
        if F_var < 1:
            p_two = 2 * p_left
            else:
            p_two = 2 * p_right
        print(f"        ratio = sd(0) / sd(1)                          f = {F_var:.4f}")
        print(f"    H0: ratio = 1                         Degrees of freedom = {df1_f}, {df2_f}")
        print(f"      Pr(F < f) = {p_left:.4f}         2*Pr(F {'<' if F_var < 1 else '>'} f) = {p_two:.4f}           Pr(F > f) = {p_right:.4f}")
    else:
        print(f"        ratio = sd(0) / sd(1)                          f =     .")
        print(f"    H0: ratio = 1                         Degrees of freedom = {n0-1}, {n1-1}")
        print(f"      Pr(F < f) =     .         2*Pr(F > f) =     .           Pr(F > f) =     .")
    print()
    
    return sd0, sd1, F_var if (sd0 > 0 and sd1 > 0) else None

# sdtest for 1904
print_sdtest(g0p, g1p, label="sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup)")

#  2. forvalue i=1900/1939 { sdtest ... } 
print("\n    forvalue i=1900/1939 { sdtest diffpatentswrt1918, by(treatmentgroup) }")
print("=" * 80)
for yr in range(1900, 1940):
    dfy = df[df["year"] == yr]
    g0y = dfy.loc[dfy["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
    g1y = dfy.loc[dfy["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()
    if len(g0y) > 0 and len(g1y) > 0:
        print_sdtest(g0y, g1y, label=f"year == {yr}")
    else:
        print(f"    year == {yr}: No observations in one or both groups")
