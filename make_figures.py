"""Figures from results.csv, written to figures/. Called at the end of run_experiments.py."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixtures_model import VARIETIES


def make_figures(run, SCEN, COLORS, df=None):
    os.makedirs("figures", exist_ok=True)
    if df is None:
        df = pd.read_csv("results.csv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#888888", "axes.labelcolor": "#333333",
                         "xtick.color": "#555555", "ytick.color": "#555555", "axes.grid": True,
                         "grid.color": "#e6e6e6", "grid.linewidth": 0.6})

    # Figure 1: the mechanism, one field in one year
    r = run(SCEN["Landrace mixture"], L=100, seed=7, keep_s=True)
    yr = 12
    d0, d1 = 365 * yr, 365 * (yr + 1)
    fig, ax = plt.subplots(figsize=(8, 3.6))
    days = np.arange(365)
    ax2 = ax.twinx()
    ax2.bar(days, r["rain"][0, d0:d1], color="#9db7d6", width=1.0, label="Rain (mm/day)")
    ax2.set_ylim(0, 4 * max(1.0, r["rain"][0, d0:d1].max())); ax2.set_ylabel("Rain (mm/day)"); ax2.grid(False); ax2.spines["top"].set_visible(False)
    ax.plot(days, r["s"][0, d0:d1], color="#333333", lw=1.6, zorder=5)
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    vcol = {"LR_early": "#2a78d6", "LR_mid": "#eb6834", "LR_late": "#1baf7a"}
    for k, (nm, c) in enumerate(vcol.items()):
        v = VARIETIES[nm]
        for pd_ in (70, 284):
            a = pd_ + int(v.stage_bounds[0] * v.maturity); b = pd_ + int(v.stage_bounds[1] * v.maturity)
            y0 = 0.92 - 0.05 * k
            ax.plot([a, min(b, 364)], [y0, y0], color=c, lw=4, solid_capstyle="round")
        ax.text(68, 0.92 - 0.05 * k, nm.replace("LR_", "").capitalize() + " landrace, flowering", color="#333333", va="center", ha="right", fontsize=8)
    ax.axhline(VARIETIES["LR_mid"].s_star, color="#999999", lw=0.8, ls="--")
    ax.text(5, VARIETIES["LR_mid"].s_star + 0.015, "s* (mid variety)", color="#666666", fontsize=8)
    ax.set_xlim(0, 364); ax.set_ylim(0, 1); ax.set_xlabel("Day of year"); ax.set_ylabel("Relative soil moisture, s")
    ax.set_title("One field, one year: soil moisture and the flowering window of each landrace", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/fig1_mechanism.png", dpi=200); plt.close(fig)

    # Figure 2: four panels
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.6))
    sc = df[(df.experiment == "scenarios_by_L") & (df.L_km == 100)]

    ax = axs[0, 0]
    for _, row in sc.iterrows():
        ax.scatter(row.p_short, row.mean_yield_t_per_ha, s=60, color=COLORS[row.scenario], edgecolor="white", lw=1.5, zorder=3)
        ax.annotate(row.scenario, (row.p_short, row.mean_yield_t_per_ha), xytext=(6, 4), textcoords="offset points", fontsize=8, color="#333333")
    ax.set_xlabel("Share of household-years below grain need")
    ax.set_ylabel("Mean yield (t/ha/yr)")
    ax.set_title("A. Average yield against shortfall risk", loc="left", fontsize=10)
    ax.set_xlim(left=0, right=sc.p_short.max() * 1.6)

    ax = axs[0, 1]
    byL = df[df.experiment == "scenarios_by_L"]
    for name, c in COLORS.items():
        d = byL[byL.scenario == name]
        ax.plot(d.L_km, d.phi, color=c, lw=2, marker="o", ms=4, label=name)
    ax.set_xscale("log"); ax.set_xlabel("Rainfall correlation length L (km, log scale)")
    ax.set_ylabel("Synchrony among farms, φ")
    ax.set_title("B. Synchrony across farms", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")

    ax = axs[1, 0]
    mt = df[df.experiment == "markets_by_theta"]
    for name, c in COLORS.items():
        d = mt[mt.scenario == name]
        ax.plot(d.theta, d.gap_buyers, color=c, lw=2, marker="o", ms=4, label=name)
    ax.set_xlabel("Market integration θ (0 isolated, 1 fully integrated)")
    ax.set_ylabel("Mean food gap, net-buyer households (t/yr)")
    ax.set_title("C. Food gap of net-buyer households", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")

    ax = axs[1, 1]
    be = df[df.experiment == "breakeven_ymax"]
    for name in ("Uniform improved", "Landrace mixture"):
        d = be[be.scenario == name]
        ax.plot(d.imp_ymax, d.p_short, color=COLORS[name], lw=2, marker="o", ms=4, label=name)
    ax.set_xlabel("Yield potential of the improved variety (t/ha per season)")
    ax.set_ylabel("Share of household-years below grain need")
    ax.set_title("D. How large must the yield gap be?", loc="left", fontsize=10)
    ax.annotate("Break-even near 2.0 t/ha", xy=(2.0, 0.49), xytext=(2.25, 0.62), fontsize=8, color="#333333",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8))
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Illustrative runs, uncalibrated parameters (64 farms, 60 years, 3 rainfall ensembles)", x=0.01, ha="left", fontsize=9, color="#666666")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig("figures/fig2_scenarios.png", dpi=200); plt.close(fig)
    print("figures written")
