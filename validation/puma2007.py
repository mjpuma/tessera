"""
Reproduction, with the tessera bucket model, of the rainfall-resolution test of
Puma, Rodriguez-Iturbe, Celia and Guswa (2007), "Implications of rainfall temporal
resolution for soil-moisture and transpiration modeling", Transport in Porous Media 68, 37-67.

Puma et al. drove a Richards-equation soil-plant model with 15-minute rainfall and with
the same rainfall aggregated to coarser intervals h (rain spread evenly over each interval),
for Bouteloua gracilis on silty loam in Colorado and Prosopis glandulosa on sandy loam in
Texas. They found that daily rainfall reproduces the 15-minute results, while 2-week and
monthly rainfall do not, and that the temporal infiltration index
    I_It = min(mu_x, K_sat) / <ET_max>_1day
falls below one at about the aggregation where errors appear.

Here the same experiment is run with the tessera point model (soilwater.py):
  * synthetic 15-minute rainfall from the rectangular pulses Poisson model with the
    parameters Puma et al. fitted to the 15-minute data (Table 2, h = 0.010 d);
  * the soil and plant parameters of their Tables 4 and 5, mapped to the loss function
    of Laio et al. (2001);
  * diurnal evapotranspiration: a 12-hour sinusoid with the daily totals of Table 4;
  * aggregation to h = 1 h, 1, 2, 4, 7, 14 and 33 days, with rain spread evenly over each
    interval and no runoff (ponded water infiltrates later), as in Puma et al.;
  * one more run with the tessera time scheme: daily totals as one pulse at the start of
    each day and constant daily evapotranspiration.

The bucket model has no vertical profile, so the levels of soil moisture differ from the
Richards-equation results. The test is whether the effect of aggregation is the same.

Run from the repository root:  python validation/puma2007.py
Units are cm and days.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from soilwater import SoilParams, drain  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
STEPS_PER_DAY = 96                      # 15-minute steps
DT = 1.0 / STEPS_PER_DAY


def s_fc_from(ks, b):
    """Field capacity as in Laio et al. (2001): where K = Ks s^(2b+3) falls to 0.05 cm/d."""
    return (0.05 / ks) ** (1.0 / (2 * b + 3))


# Puma et al. (2007), Tables 2, 4 and 5
SITES = {
    "Colorado (B. gracilis, silty loam)": dict(
        n_days=957, lam=1.369, mu_L=0.003, mu_X=51.129,
        n=0.468, z_r=30.0, k_s=33.0, b=4.98, s_h=0.137, s_w=0.288, s_star=0.454, et_max=0.37, e_w=0.01,
        puma_mean={"15 min": 0.438, "1 day": 0.436, "14 days": 0.411, "33 days": 0.394},
        puma_sd={"15 min": 0.073, "1 day": 0.072, "14 days": 0.054, "33 days": 0.042}),
    "Texas (P. glandulosa, sandy loam)": dict(
        n_days=660, lam=0.538, mu_L=0.010, mu_X=58.087,
        n=0.43, z_r=150.0, k_s=82.2, b=4.9, s_h=0.14, s_w=0.290, s_star=0.440, et_max=0.442, e_w=0.02,
        puma_mean={"15 min": 0.454, "1 day": 0.447, "14 days": 0.402, "33 days": 0.377},
        puma_sd={"15 min": 0.084, "1 day": 0.082, "14 days": 0.057, "33 days": 0.046}),
}
LEVELS = [("15 min", 1), ("1 h", 4), ("1 day", 96), ("2 days", 192), ("4 days", 384),
          ("7 days", 672), ("14 days", 1344), ("33 days", 3168)]


def rpp_rain(n_days, lam, mu_L, mu_X, n_real, seed):
    """Rectangular pulses Poisson rainfall, as depth (cm) in each 15-minute step: (n_real, n_steps)."""
    rng = np.random.default_rng(seed)
    n_steps = n_days * STEPS_PER_DAY
    rain = np.zeros((n_real, n_steps))
    for r in range(n_real):
        k = rng.poisson(lam * n_days)
        t0 = rng.uniform(0, n_days, k)
        dur = rng.exponential(mu_L, k)
        x = rng.exponential(mu_X, k)
        t1 = np.minimum(t0 + dur, n_days)
        b0 = (t0 / DT).astype(int)
        b1 = np.minimum((t1 / DT).astype(int), n_steps - 1)
        one = b0 == b1
        np.add.at(rain[r], b0[one], x[one] * (t1[one] - t0[one]))
        for i in np.nonzero(~one)[0]:
            edges = np.arange(b0[i], b1[i] + 2) * DT
            lo = np.maximum(edges[:-1], t0[i]); hi = np.minimum(edges[1:], t1[i])
            rain[r, b0[i]:b1[i] + 1] += x[i] * np.maximum(hi - lo, 0)
    return rain


def aggregate(rain, k):
    """Spread the total of each block of k steps evenly over the block."""
    n_real, n = rain.shape
    m = n // k * k
    out = rain.copy()
    blk = rain[:, :m].reshape(n_real, -1, k).sum(-1, keepdims=True) / k
    out[:, :m] = np.repeat(blk, k, axis=-1).reshape(n_real, m)
    if m < n:
        out[:, m:] = rain[:, m:].sum(1, keepdims=True) / (n - m)
    return out


def diurnal(step):
    """Atmospheric demand over the day: zero at night, a 12-hour sinusoid by day, mean 1."""
    f = (step % STEPS_PER_DAY + 0.5) / STEPS_PER_DAY
    return np.pi * np.sin(2 * np.pi * (f - 0.25)) if 0.25 <= f < 0.75 else 0.0


def run(rain, site, scheme="15 min", spin_days=33):
    """Integrate the bucket. scheme '15 min': 15-minute steps, diurnal ET, ponding.
    scheme 'tessera': daily pulse at the start of each day, constant daily ET."""
    p = site
    soil = SoilParams(n=p["n"], z_r=p["z_r"], k_s=p["k_s"], beta=2 * p["b"] + 4, s_h=p["s_h"],
                      s_fc=s_fc_from(p["k_s"], p["b"]), e_w=p["e_w"], e_bare=0.0)
    nzr = soil.nzr
    n_real, n_steps = rain.shape
    cov = np.ones((n_real, 1))
    e_max, s_star, s_w = np.array([p["et_max"]]), np.array([p["s_star"]]), np.array([p["s_w"]])
    s = np.full(n_real, p["s_star"])
    pond = np.zeros(n_real)
    keep = spin_days * STEPS_PER_DAY
    s_samples, t_ratio = [], []
    tot = dict(infil=0.0, transp=0.0)
    if scheme == "tessera":
        daily = rain.reshape(n_real, -1, STEPS_PER_DAY).sum(-1)
        for d in range(daily.shape[1]):
            room = nzr * (1 - s)
            water = daily[:, d] + pond
            inf = np.minimum(water, room); pond = water - inf
            s = s + inf / nzr
            # sample at a random time within the day so samples represent the time average
            u = (d * 0.618034) % 1.0
            s_mid, tr1, *_ = drain(s, u, cov, e_max, s_star, s_w, soil, max_ds=0.02)
            s, tr2, *_ = drain(s_mid, 1.0 - u, cov, e_max, s_star, s_w, soil, max_ds=0.02)
            if d >= spin_days:
                s_samples.append(s_mid)
                t_ratio.append((tr1 + tr2).ravel() / (p["et_max"] - p["e_w"]))
                tot["infil"] += inf.sum(); tot["transp"] += (tr1 + tr2).sum()
    else:
        for k in range(n_steps):
            water = rain[:, k] + pond
            room = nzr * (1 - s)
            inf = np.minimum(water, room); pond = water - inf
            s = s + inf / nzr
            g = diurnal(k)
            s, tr, ev, lk = drain(s, DT, cov, e_max, s_star, s_w, soil, et_scale=g, max_ds=0.01)
            if k >= keep:
                if k % 4 == 0:
                    s_samples.append(s)
                if g > 0:
                    t_ratio.append(tr.ravel() / (g * (p["et_max"] - p["e_w"]) * DT))
                tot["infil"] += inf.sum(); tot["transp"] += tr.sum()
    S = np.concatenate(s_samples)
    TR = np.concatenate(t_ratio)
    return S, TR, tot["transp"] / tot["infil"]


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
    rows, store = [], {}
    for si, (name, site) in enumerate(SITES.items()):
        rain15 = rpp_rain(site["n_days"], site["lam"], site["mu_L"], site["mu_X"], n_real=100, seed=7 + si)
        ref = None
        for lab, k in LEVELS + [("tessera daily scheme", None)]:
            r = rain15 if k in (1, None) else aggregate(rain15, k)
            S, TR, frac = run(r, site, scheme="tessera" if k is None else "15 min")
            if lab == "15 min":
                ref = S
            kk = 96 if k is None else k
            intens = aggregate(rain15, kk) if kk > 1 else rain15
            mu_x = (intens[intens > 0] / DT).mean()                       # cm/day while raining
            ks = np.sort(S); kr = np.sort(ref)
            grid = np.linspace(0, 1, 2001)
            ksd = np.max(np.abs(np.searchsorted(ks, grid) / len(ks) - np.searchsorted(kr, grid) / len(kr)))
            rows.append(dict(site=name, rainfall=lab, h_days=kk / 96, mean_S=S.mean(), sd_S=S.std(),
                             KS_vs_15min=ksd, mean_T_over_Tmax=TR.mean(), share_infiltration_transpired=frac,
                             I_It=min(mu_x, site["k_s"]) / site["et_max"], R_I=site["k_s"] / mu_x,
                             puma2007_mean_S=site["puma_mean"].get(lab, np.nan),
                             puma2007_sd_S=site["puma_sd"].get(lab, np.nan)))
            store[(name, lab)] = (S, TR)
            print(name, lab, "done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "puma2007_results.csv"), index=False)
    pd.set_option("display.width", 220)
    print(df.round(3).to_string(index=False))

    # Figure: soil moisture pdfs in the layout of Puma et al. (2007), Fig. 5
    col = {"15 min": "#2a78d6", "1 day": "#1baf7a", "tessera daily scheme": "#e87ba4",
           "14 days": "#eda100", "33 days": "#e34948"}
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.6), sharex=True)
    bins = np.linspace(0.2, 0.8, 61)
    for i, name in enumerate(SITES):
        for j, labs in enumerate((("15 min", "1 day", "tessera daily scheme"), ("15 min", "14 days", "33 days"))):
            ax = axs[i, j]
            for lab in labs:
                h, e = np.histogram(store[(name, lab)][0], bins=bins, density=True)
                ax.plot(0.5 * (e[1:] + e[:-1]), h, color=col[lab], lw=1.8 if lab == "15 min" else 1.4,
                        ls="--" if lab == "15 min" else "-", label=lab)
            ax.set_title(f"{name.split(' (')[0]}: 15-min rainfall against {'daily' if j == 0 else 'coarser'} rainfall",
                         loc="left", fontsize=9)
            ax.legend(frameon=False, fontsize=7.5)
            if i == 1: ax.set_xlabel("Root-zone relative soil moisture, s")
            if j == 0: ax.set_ylabel("Probability density")
    fig.text(0.01, 0.005, "soilwater.py bucket driven by synthetic 15-min rainfall (rectangular pulses Poisson, "
             "parameters of Puma et al. 2007, Table 2) and the same rainfall aggregated. Layout follows their Fig. 5.",
             fontsize=7, color="#666666")
    fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(os.path.join(OUT, "puma_fig5_pdfs.png"), dpi=180)
    plt.close(fig)

    # Figure: error against aggregation, and the temporal infiltration index
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for name, c in zip(SITES, ("#2a78d6", "#eb6834")):
        d = df[(df.site == name) & (df.rainfall != "tessera daily scheme")]
        short = name.split(" (")[0]
        axs[0].plot(d.h_days, d.KS_vs_15min, "o-", color=c, lw=1.6, ms=4, label=short)
        axs[1].plot(d.h_days, d.I_It, "o-", color=c, lw=1.6, ms=4, label=short)
        t = df[(df.site == name) & (df.rainfall == "tessera daily scheme")]
        axs[0].plot(t.h_days, t.KS_vs_15min, "D", color=c, ms=6, mfc="white", label=f"{short}, tessera daily scheme")
    axs[0].set_xscale("log"); axs[1].set_xscale("log"); axs[1].set_yscale("log")
    axs[1].axhline(1.0, color="#888888", lw=0.9, ls="--")
    axs[1].text(0.012, 1.15, "I_It = 1", fontsize=7.5, color="#666666")
    axs[0].set_xlabel("Rainfall aggregation h (days)"); axs[1].set_xlabel("Rainfall aggregation h (days)")
    axs[0].set_ylabel("Distance from 15-min pdf (KS)"); axs[1].set_ylabel("Temporal infiltration index I_It")
    axs[0].set_title("A. Error in the soil moisture pdf", loc="left", fontsize=9)
    axs[1].set_title("B. Puma et al. (2007) a priori index", loc="left", fontsize=9)
    axs[0].legend(frameon=False, fontsize=7); axs[1].legend(frameon=False, fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "puma_resolution.png"), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
