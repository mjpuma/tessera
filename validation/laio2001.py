"""
Validation of soilwater.py against Laio, Porporato, Ridolfi and Rodriguez-Iturbe (2001),
"Plants in water-controlled ecosystems ... II. Probabilistic soil moisture dynamics",
Advances in Water Resources 24, 707-723.

Three tests, all for a single variety at full cover, where soilwater.py reduces exactly
to the paper's loss function (Eq. 16):

  1. Drydown: the numerical integrator against the closed-form decay of Eq. 20 (Fig. 7).
  2. Steady-state pdf: Monte Carlo simulation against the analytical pdf of Eq. 29,
     for the twelve cases of Fig. 8.
  3. Water balance: simulated partition of rainfall against Eqs. 36-44 (Fig. 10a).

Run from the repository root:  python validation/laio2001.py
Units are cm and days, as in the paper.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from soilwater import SoilParams, infiltrate, drain  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
E_MAX, E_W = 0.45, 0.01          # cm/day, as in Figs. 7, 8 and 10

# Laio et al. (2001), Table 1
SOILS = {
    "loamy sand": dict(n=0.42, beta=12.7, s_h=0.08, s_w=0.11, s_star=0.31, s_fc=0.52, k_s=100.0),
    "loam":       dict(n=0.45, beta=14.8, s_h=0.19, s_w=0.24, s_star=0.57, s_fc=0.65, k_s=20.0),
}


def soil_params(name, z_r):
    p = SOILS[name]
    return SoilParams(n=p["n"], z_r=z_r, k_s=p["k_s"], beta=p["beta"], s_h=p["s_h"], s_fc=p["s_fc"],
                      e_w=E_W, e_bare=0.0), p


# ------------------------------------------------------------------ analytical results

def rho(s, soil, p):
    """Normalized loss function, Laio et al. Eq. 16 (1/day)."""
    nzr = soil.n * soil.z_r
    eta_w, eta = E_W / nzr, E_MAX / nzr
    m = soil.k_s / (nzr * (np.exp(soil.beta * (1 - soil.s_fc)) - 1))
    out = np.zeros_like(s)
    a = (s > soil.s_h) & (s <= p["s_w"])
    out[a] = eta_w * (s[a] - soil.s_h) / (p["s_w"] - soil.s_h)
    b = (s > p["s_w"]) & (s <= p["s_star"])
    out[b] = eta_w + (eta - eta_w) * (s[b] - p["s_w"]) / (p["s_star"] - p["s_w"])
    c = (s > p["s_star"]) & (s <= soil.s_fc)
    out[c] = eta
    d = s > soil.s_fc
    out[d] = eta + m * (np.exp(soil.beta * (s[d] - soil.s_fc)) - 1)
    return out


def analytic_pdf(soil, p, lam, alpha, delta=0.0, n=40000):
    """Steady-state pdf, Eq. 29, evaluated numerically; returns grid, pdf, cdf."""
    nzr = soil.n * soil.z_r
    gamma = nzr / alpha
    lam_p = lam * np.exp(-delta / alpha)
    # grid dense near s_h, where 1/rho is large
    u = np.linspace(0, 1, n)
    s = soil.s_h + (1 - soil.s_h) * (1e-6 + (1 - 1e-6) * u ** 3)
    r = rho(s, soil, p)
    inv = 1.0 / r
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(s))])
    logp = -np.log(r) - gamma * s + lam_p * integ
    logp -= logp.max()
    pdf = np.exp(logp)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(s))])
    pdf /= cdf[-1]
    cdf /= cdf[-1]
    return s, pdf, cdf


def analytic_drydown(t, soil, p, s0=1.0):
    """Closed-form soil moisture decay in the absence of rain, Eqs. 20-21."""
    nzr = soil.n * soil.z_r
    eta_w, eta = E_W / nzr, E_MAX / nzr
    m = soil.k_s / (nzr * (np.exp(soil.beta * (1 - soil.s_fc)) - 1))
    b, sfc, sst, sw, sh = soil.beta, soil.s_fc, p["s_star"], p["s_w"], soil.s_h
    A = eta - m + m * np.exp(b * (s0 - sfc))
    t_fc = (b * (sfc - s0) + np.log(A / eta)) / (b * (m - eta))
    t_st = (sfc - sst) / eta + t_fc
    t_w = (sst - sw) / (eta - eta_w) * np.log(eta / eta_w) + t_st
    s = np.empty_like(t)
    i1 = t < t_fc
    s[i1] = s0 - np.log((A * np.exp(b * (eta - m) * t[i1]) - m * np.exp(b * (s0 - sfc))) / (eta - m)) / b
    i2 = (t >= t_fc) & (t < t_st)
    s[i2] = sfc - eta * (t[i2] - t_fc)
    i3 = (t >= t_st) & (t < t_w)
    s[i3] = sw + (sst - sw) * (eta / (eta - eta_w) * np.exp(-(eta - eta_w) / (sst - sw) * (t[i3] - t_st)) - eta_w / (eta - eta_w))
    i4 = t >= t_w
    s[i4] = sh + (sw - sh) * np.exp(-eta_w / (sw - sh) * (t[i4] - t_w))
    return s


def analytic_balance(soil, p, lam, alpha, delta):
    """Mean water balance as fractions of rainfall, Eqs. 36-44."""
    s, pdf, cdf = analytic_pdf(soil, p, lam, alpha, delta)
    nzr = soil.n * soil.z_r
    lam_p = lam * np.exp(-delta / alpha)
    r = rho(s, soil, p)
    chi = nzr * np.trapezoid(r * pdf, s)
    eta = E_MAX / nzr
    P_star = np.interp(p["s_star"], s, cdf)
    e_ns = E_MAX * (1 - P_star)
    leak = nzr * np.trapezoid(np.where(s > soil.s_fc, r - eta, 0.0) * pdf, s)
    e_s = chi - e_ns - leak
    rain = alpha * lam
    return dict(I=alpha * (lam - lam_p) / rain, Q=(alpha * lam_p - chi) / rain,
                L=leak / rain, E_ns=e_ns / rain, E_s=e_s / rain)


# ------------------------------------------------------------------ simulation

def simulate_poisson(soil, p, lam, alpha, delta=0.0, n_chains=200, n_events=3000, burn=100,
                     max_ds=0.01, seed=0):
    """Continuous-time marked Poisson rainfall. Soil moisture is sampled just before each
    storm, which by the PASTA property (Poisson arrivals see time averages) gives samples
    from the time-averaged pdf. Also returns long-run water-balance totals."""
    rng = np.random.default_rng(seed)
    cov = np.ones((n_chains, 1))
    e_max, s_star, s_w = np.array([E_MAX]), np.array([p["s_star"]]), np.array([p["s_w"]])
    s = np.full(n_chains, 0.5 * (p["s_star"] + soil.s_fc))
    samples = []
    tot = dict(rain=0.0, I=0.0, Q=0.0, ET=0.0, L=0.0, time=0.0, unstressed_time=0.0)
    for k in range(n_events):
        dt = rng.exponential(1.0 / lam, n_chains)
        depth = rng.exponential(alpha, n_chains)
        through = np.maximum(depth - delta, 0.0)
        s, q = infiltrate(s, through, soil)
        s_new, tr, ev, lk = drain(s, dt, cov, e_max, s_star, s_w, soil, max_ds=max_ds)
        if k >= burn:
            samples.append(s_new.copy())
            tot["rain"] += depth.sum(); tot["I"] += (depth - through).sum(); tot["Q"] += q.sum()
            tot["ET"] += tr.sum() + ev.sum(); tot["L"] += lk.sum(); tot["time"] += dt.sum()
        s = s_new
    samples = np.concatenate(samples)
    p_unstressed = (samples > p["s_star"]).mean()           # PASTA: share of time with s > s*
    rain = tot["rain"]
    e_ns = E_MAX * p_unstressed * tot["time"]
    return samples, dict(I=tot["I"] / rain, Q=tot["Q"] / rain, L=tot["L"] / rain,
                         E_ns=e_ns / rain, E_s=(tot["ET"] - e_ns) / rain)


def simulate_daily(soil, p, lam, alpha, n_chains=200, n_days=3000, burn=365, max_ds=0.02, seed=1):
    """The tessera scheme: the day's storms (Poisson number with mean lam) are summed and
    applied as one pulse at the start of the day, then one day of losses. Samples are taken at a random time
    within each day so that they represent the time-averaged pdf."""
    rng = np.random.default_rng(seed)
    cov = np.ones((n_chains, 1))
    e_max, s_star, s_w = np.array([E_MAX]), np.array([p["s_star"]]), np.array([p["s_w"]])
    s = np.full(n_chains, 0.5 * (p["s_star"] + soil.s_fc))
    samples = []
    for d in range(burn + n_days):
        count = rng.poisson(lam, n_chains)
        rain = np.where(count > 0, rng.gamma(np.maximum(count, 1), alpha), 0.0)
        s, _ = infiltrate(s, rain, soil)
        # sample at a random time within the day, so samples represent the time average
        u = rng.random(n_chains)
        s, *_ = drain(s, u, cov, e_max, s_star, s_w, soil, max_ds=max_ds)
        if d >= burn:
            samples.append(s.copy())
        s, *_ = drain(s, 1.0 - u, cov, e_max, s_star, s_w, soil, max_ds=max_ds)
    return np.concatenate(samples)


# ------------------------------------------------------------------ run

def ks_stat(samples, s_grid, cdf):
    x = np.sort(samples)
    ecdf = np.arange(1, len(x) + 1) / len(x)
    return np.max(np.abs(ecdf - np.interp(x, s_grid, cdf)))


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
    rows = []
    COL = {"loamy sand": "#2a78d6", "loam": "#eb6834"}

    # 1. Drydown, Fig. 7
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    t = np.linspace(0, 60, 601)
    for ax, name in zip(axs, ["loam", "loamy sand"]):
        for z_r, ls in zip((30, 60, 90), ("-", "--", "-.")):
            soil, p = soil_params(name, z_r)
            s_an = analytic_drydown(t, soil, p)
            for max_ds, step, lab in ((0.005, 0.1, "fine"), (0.02, 1.0, "daily")):
                s = np.array([1.0]); out = [1.0]; tt = [0.0]
                for k in range(int(60 / step)):
                    s, *_ = drain(s, step, np.ones((1, 1)), np.array([E_MAX]), np.array([p["s_star"]]),
                                  np.array([p["s_w"]]), soil, max_ds=max_ds)
                    out.append(s[0]); tt.append((k + 1) * step)
                err = np.max(np.abs(np.array(out) - np.interp(tt, t, s_an)))
                rows.append(dict(test="drydown (Eq. 20)", soil=name, Zr_cm=z_r, scheme=lab, max_abs_error=err))
                if lab == "daily":
                    ax.plot(tt, out, "o", ms=2.5, color=COL[name], alpha=0.7)
            ax.plot(t, s_an, ls, color="#333333", lw=1.2, label=f"Eq. 20, Zr = {z_r} cm")
        ax.set_title(f"{name.capitalize()}: drydown from saturation", loc="left", fontsize=9.5)
        ax.set_xlabel("Days without rain")
        ax.legend(frameon=False, fontsize=7.5)
    axs[0].set_ylabel("Relative soil moisture, s")
    fig.text(0.01, 0.01, "Lines: closed form (Laio et al. 2001, Eq. 20). Dots: soilwater.py with daily steps.",
             fontsize=7.5, color="#666666")
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(os.path.join(OUT, "laio_fig7_drydown.png"), dpi=180)
    plt.close(fig)

    # 2. Steady-state pdfs, Fig. 8
    fig, axs = plt.subplots(3, 2, figsize=(8.5, 8.5), sharex=True)
    for i, lam in enumerate((0.1, 0.2, 0.5)):
        for j, z_r in enumerate((30, 90)):
            ax = axs[i, j]
            for name in ("loamy sand", "loam"):
                soil, p = soil_params(name, z_r)
                s_g, pdf, cdf = analytic_pdf(soil, p, lam, 1.5)
                samp, _ = simulate_poisson(soil, p, lam, 1.5, seed=10 * i + j)
                daily = simulate_daily(soil, p, lam, 1.5, seed=100 + 10 * i + j)
                mean_an = np.trapezoid(s_g * pdf, s_g)
                sd_an = np.sqrt(np.trapezoid((s_g - mean_an) ** 2 * pdf, s_g))
                for lab, x in (("Poisson storms", samp), ("daily pulses (tessera scheme)", daily)):
                    rows.append(dict(test="steady-state pdf (Eq. 29)", soil=name, Zr_cm=z_r, lambda_per_day=lam,
                                     scheme=lab, mean_analytic=mean_an, mean_sim=x.mean(),
                                     sd_analytic=sd_an, sd_sim=x.std(), KS=ks_stat(x, s_g, cdf)))
                ax.plot(s_g, pdf, color=COL[name], lw=1.6, label=f"{name}, Eq. 29")
                h, e = np.histogram(samp, bins=np.linspace(0, 1, 101), density=True)
                ax.step(e[:-1], h, where="post", color=COL[name], lw=0.8, alpha=0.55)
            ax.set_xlim(0, 1); ax.set_ylim(0, 11)
            ax.set_title(f"λ = {lam} per day, Zr = {z_r} cm", loc="left", fontsize=9)
            if i == 2: ax.set_xlabel("Relative soil moisture, s")
            if j == 0: ax.set_ylabel("p(s)")
    axs[0, 0].legend(frameon=False, fontsize=7.5)
    fig.text(0.01, 0.005, "Smooth lines: analytical pdf (Laio et al. 2001, Eq. 29). Thin steps: soilwater.py "
             "Monte Carlo (200 chains, 2,900 storms each). α = 1.5 cm, Emax = 0.45 cm/d, Ew = 0.01 cm/d.",
             fontsize=7, color="#666666")
    fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(os.path.join(OUT, "laio_fig8_pdfs.png"), dpi=180)
    plt.close(fig)

    # 3. Water balance, Fig. 10a: shallow loam, Zr = 30 cm, alpha = 2 cm, Delta = 0.2 cm
    soil, p = soil_params("loam", 30)
    lams = np.round(np.arange(0.1, 0.71, 0.1), 2)
    an, sim = [], []
    for k, lam in enumerate(lams):
        a = analytic_balance(soil, p, lam, 2.0, 0.2)
        _, b = simulate_poisson(soil, p, lam, 2.0, delta=0.2, n_events=2000, seed=500 + k)
        an.append(a); sim.append(b)
        for comp in ("I", "Q", "L", "E_ns", "E_s"):
            rows.append(dict(test="water balance (Eqs. 36-44)", soil="loam", Zr_cm=30, lambda_per_day=lam,
                             scheme="Poisson storms", component=comp, share_analytic=a[comp], share_sim=b[comp]))
    fig, ax = plt.subplots(figsize=(6, 4))
    order = ["E_s", "E_ns", "L", "Q", "I"]
    cum_an = np.zeros(len(lams)); cum_sim = np.zeros(len(lams))
    for comp in order:
        cum_an += np.array([a[comp] for a in an]); cum_sim += np.array([b[comp] for b in sim])
        ax.plot(lams, 100 * cum_an, color="#333333", lw=1.2)
        ax.plot(lams, 100 * cum_sim, "o", color="#2a78d6", ms=4)
        ax.text(lams[-1] + 0.01, 100 * cum_an[-1] - 3, {"E_s": "<Es>", "E_ns": "<Ens>", "L": "<L>", "Q": "<Q>", "I": "<I>"}[comp],
                fontsize=8, color="#333333")
    ax.set_xlabel("λ (storms per day)"); ax.set_ylabel("Cumulative % of rainfall")
    ax.set_ylim(0, 100); ax.set_xlim(0.08, 0.78)
    ax.set_title("Water balance, loam, Zr = 30 cm (cf. Laio et al. 2001, Fig. 10a)", loc="left", fontsize=9)
    fig.text(0.01, 0.01, "Lines: Eqs. 36-44 with the analytical pdf. Dots: soilwater.py Monte Carlo.", fontsize=7.5, color="#666666")
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(os.path.join(OUT, "laio_fig10_waterbalance.png"), dpi=180)
    plt.close(fig)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "laio2001_results.csv"), index=False)
    pd.set_option("display.width", 200)
    print(df[df.test.str.startswith("drydown")].round(4).to_string(index=False))
    pdf_rows = df[df.test.str.startswith("steady")][["soil", "Zr_cm", "lambda_per_day", "scheme", "mean_analytic", "mean_sim", "sd_analytic", "sd_sim", "KS"]]
    print(pdf_rows.round(4).to_string(index=False))
    wb = df[df.test.str.startswith("water")]
    print("water balance: max |share_sim - share_analytic| =", round((wb.share_sim - wb.share_analytic).abs().max(), 4))


if __name__ == "__main__":
    main()
