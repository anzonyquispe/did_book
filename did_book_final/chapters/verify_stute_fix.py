"""
Verification script for corrected Stute test p-values.
Matches Stata's stute_test.ado: Mammen (1993) wild bootstrap + OLS re-fit.

Run this BEFORE modifying ch07.qmd.
Expected: p-values within ~±0.02 of Stata reference (different RNG).
"""
import numpy as np
import pandas as pd

# ============================================================
# Corrected functions
# ============================================================

def stute_test_cvm(Y, D, df, order=1, seed=1, brep=500):
    """
    Stute (1997) Cramer-von Mises test — matches Stata stute_test.
    Mammen (1993) wild bootstrap with OLS re-fit.
    """
    sub = df[[Y, D]].dropna()
    y = sub[Y].values.astype(float)
    d = sub[D].values.astype(float)
    n = len(y)

    # Sort by (D, Y) — matches Stata's sort(data, (2, 1))
    idx = np.lexsort((y, d))
    y, d = y[idx], d[idx]

    # Design matrix: [1, d, d^2, ..., d^order]
    X = np.column_stack([d**k for k in range(order + 1)])

    # OLS fit and residuals
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ beta

    # CvM statistic: (1/n^2) * sum(cumsum(e)^2)
    cum_e = np.cumsum(e)
    cvm = np.sum(cum_e**2) / n**2

    # Uniform draws -> Mammen two-point distribution
    rng = np.random.RandomState(seed)
    V = rng.uniform(size=(n, brep))

    c1 = (np.sqrt(5) + 1) / (2 * np.sqrt(5))   # ~0.7236
    c2 = (1 - np.sqrt(5)) / 2                    # ~-0.6180
    c3 = np.sqrt(5)                               # ~2.2361
    F = c2 + c3 * (V > c1).astype(float)
    # F takes values (1-sqrt5)/2 or (1+sqrt5)/2 — Mammen two-point

    # Wild bootstrap: Y* = X*beta + F * e, then RE-FIT OLS
    Xb = (X @ beta)[:, None]       # n x 1
    ec = e[:, None]                  # n x 1
    Y_st = Xb + F * ec              # n x brep
    b_st = np.linalg.lstsq(X, Y_st, rcond=None)[0]  # (order+1) x brep
    e_st = Y_st - X @ b_st          # n x brep — projected bootstrap residuals

    # Bootstrap CvM statistics
    bres = np.sum(np.cumsum(e_st, axis=0)**2, axis=0) / n**2  # brep vector

    # p-value: strict > , proportion (matching Stata)
    pval = np.mean(bres > cvm)

    return cvm, pval


def stute_test_panel(Y, D, G, T, df, order=1, seed=1, brep=500):
    """
    Panel Stute test with joint statistic — matches Stata stute_test.
    V matrix shared across time periods (group-level clustering).
    """
    sub = df[[Y, D, G, T]].dropna()
    times = sorted(sub[T].unique())
    GG = sub[G].nunique()

    # Single V matrix for all periods (group-level bootstrap)
    rng = np.random.RandomState(seed)
    V = rng.uniform(size=(GG, brep))

    c1 = (np.sqrt(5) + 1) / (2 * np.sqrt(5))
    c2 = (1 - np.sqrt(5)) / 2
    c3 = np.sqrt(5)

    t_tot = 0.0
    t_boot = np.zeros(brep)
    results = []

    for t in times:
        dft = sub[sub[T] == t]
        y = dft[Y].values.astype(float)
        d = dft[D].values.astype(float)
        n = len(y)
        assert n == GG, f"Unbalanced panel at T={t}: {n} obs vs {GG} groups"

        idx = np.lexsort((y, d))
        y, d = y[idx], d[idx]

        X = np.column_stack([d**k for k in range(order + 1)])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        e = y - X @ beta

        cum_e = np.cumsum(e)
        cvm = np.sum(cum_e**2) / n**2

        # Same V for all periods — rows mapped to D-sorted observations
        F = c2 + c3 * (V > c1).astype(float)

        Xb = (X @ beta)[:, None]
        ec = e[:, None]
        Y_st = Xb + F * ec
        b_st = np.linalg.lstsq(X, Y_st, rcond=None)[0]
        e_st = Y_st - X @ b_st

        bres = np.sum(np.cumsum(e_st, axis=0)**2, axis=0) / n**2
        pval = np.mean(bres > cvm)

        results.append((t, cvm, pval))
        t_tot += cvm
        t_boot += bres

    joint_pval = np.mean(t_boot > t_tot)
    return results, t_tot, joint_pval


# ============================================================
# Verification
# ============================================================

if __name__ == "__main__":
    url = ("https://raw.githubusercontent.com/Credible-Answers/did_book/main/"
           "cc_xd_didtextbook_2025_9_30/Data%20sets/Pierce%20and%20Schott%202016/"
           "pierce_schott_didtextbook.dta")
    df = pd.read_stata(url)

    stata_ref_gq4 = {2001: 0.042, 2002: 0.134, 2004: 0.478, 2005: 0.716}
    stata_ref_gq7 = {1998: 0.302, 1997: 0.508}
    stata_ref_gq9 = {2001: 0.378, 2002: 0.062, 2004: 0.384, 2005: 0.530}

    # ── GQ4: Stute Tests of Linearity (order=1) ──────────────
    print("=" * 65)
    print("GQ4: Stute Tests of Linearity (order=1)")
    print("=" * 65)
    print(f"{'Variable':<18} {'CvM stat':>12} {'p-value':>9} {'Stata ref':>11}")
    print("-" * 52)
    for yr in [2001, 2002, 2004, 2005]:
        stat, pv = stute_test_cvm(f'delta{yr}', 'ntrgap', df, order=1, seed=1)
        ref = stata_ref_gq4[yr]
        print(f"delta{yr:<13} {stat:12.7f} {pv:9.3f} {ref:11.3f}")

    # Panel joint test
    print("\nPanel Joint Test:")
    id_vars = ['indusid', 'ntrgap']
    delta_cols = [f'delta{y}' for y in [2001, 2002, 2004, 2005]]
    df_long = df[id_vars + delta_cols].melt(
        id_vars=id_vars, value_vars=delta_cols,
        var_name='year_str', value_name='delta')
    df_long['year'] = df_long['year_str'].str.replace('delta', '').astype(int)

    res, jstat, jpv = stute_test_panel(
        'delta', 'ntrgap', 'indusid', 'year', df_long,
        order=1, seed=1)
    for t, s, p in res:
        print(f"  {t:>6}  {s:12.7f}  {p:.3f}")
    print(f"  Joint: {jstat:.4f} ({jpv:.4f})   [Stata: 0.0089 (0.5400)]")

    # ── GQ7: Stute Tests of Mean Independence (order=0) ──────
    print("\n" + "=" * 65)
    print("GQ7: Stute Tests of Mean Independence (order=0)")
    print("=" * 65)
    print(f"{'Variable':<24} {'CvM stat':>12} {'p-value':>9} {'Stata ref':>11}")
    print("-" * 58)
    for yr in [1998, 1997]:
        stat, pv = stute_test_cvm(f'deltalintrend{yr}', 'ntrgap', df, order=0, seed=1)
        ref = stata_ref_gq7[yr]
        print(f"deltalintrend{yr:<11} {stat:12.7f} {pv:9.3f} {ref:11.3f}")

    # Panel joint test
    print("\nPanel Joint Test:")
    lt_cols = [f'deltalintrend{y}' for y in [1997, 1998]]
    df_lt = df[id_vars + lt_cols].melt(
        id_vars=id_vars, value_vars=lt_cols,
        var_name='year_str', value_name='deltalintrend')
    df_lt['year'] = df_lt['year_str'].str.replace('deltalintrend', '').astype(int)

    res, jstat, jpv = stute_test_panel(
        'deltalintrend', 'ntrgap', 'indusid', 'year', df_lt,
        order=0, seed=1)
    for t, s, p in res:
        print(f"  {t:>6}  {s:12.7f}  {p:.3f}")
    print(f"  Joint: {jstat:.4f} ({jpv:.4f})   [Stata: 0.0012 (0.4720)]")

    # ── GQ9: Stute Tests of Linearity with Linear Trends ─────
    print("\n" + "=" * 65)
    print("GQ9: Stute Tests of Linearity w/ Linear Trends (order=1)")
    print("=" * 65)
    print(f"{'Variable':<24} {'CvM stat':>12} {'p-value':>9} {'Stata ref':>11}")
    print("-" * 58)
    for yr in [2001, 2002, 2004, 2005]:
        stat, pv = stute_test_cvm(f'deltalintrend{yr}', 'ntrgap', df, order=1, seed=1)
        ref = stata_ref_gq9[yr]
        print(f"deltalintrend{yr:<11} {stat:12.7f} {pv:9.3f} {ref:11.3f}")

    # Panel joint test
    print("\nPanel Joint Test:")
    lt_cols = [f'deltalintrend{y}' for y in [2001, 2002, 2004, 2005]]
    df_lt2 = df[id_vars + lt_cols].melt(
        id_vars=id_vars, value_vars=lt_cols,
        var_name='year_str', value_name='deltalintrend')
    df_lt2['year'] = df_lt2['year_str'].str.replace('deltalintrend', '').astype(int)

    res, jstat, jpv = stute_test_panel(
        'deltalintrend', 'ntrgap', 'indusid', 'year', df_lt2,
        order=1, seed=1)
    for t, s, p in res:
        print(f"  {t:>6}  {s:12.7f}  {p:.3f}")
    print(f"  Joint: {jstat:.4f} ({jpv:.4f})   [Stata: 0.0119 (0.4020)]")
