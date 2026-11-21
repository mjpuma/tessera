"""
Benchmark the bucket soil model against a Richards-equation soil column (SoilWaterFlow).

Both models get the same farms, the same rainfall, the same varieties and the same
yield response; only the soil differs. The bucket has one well-mixed root zone shared by
all varieties. The column resolves water in depth, and each variety draws water through
its own roots (richards.ROOTS, illustrative). Run from the repository root:

    python validation/richards_vs_bucket.py      # about 5 minutes on 2 cores

Writes validation/richards_vs_bucket.csv and validation/richards_vs_bucket.png.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mixtures_model import (VARIETIES, Soil, Climate, Disease, make_fields, generate_rain,  # noqa: E402
                            simulate, field_yields, annual, partition)
from richards import simulate_richards, ROOTS  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
N_FIELDS, N_YEARS, SEED = 16, 11, 7          # first year is spin-up
SOIL, DIS = Soil(), Disease()
XY, AREA = make_fields(N_FIELDS)
SHARES = {
    "Uniform improved": np.tile([0, 0, 0, 1.0], (N_FIELDS, 1)),
    "Landrace mixture": np.tile([1 / 3, 1 / 3, 1 / 3, 0], (N_FIELDS, 1)),
}


def main():
    rain, _ = generate_rain(XY, N_YEARS, Climate(), seed=SEED)
    rows, ya, strace, timing = [], {}, {}, {}
    for scen, sh in SHARES.items():
        for model, fn in (("bucket", simulate), ("Richards", simulate_richards)):
            t0 = time.time()
            sim = fn(rain, sh, SOIL, VARIETIES)
            timing[model] = timing.get(model, 0.0) + time.time() - t0
            Ya = annual(field_yields(sim, sh, DIS, seed=SEED + 100))[:, 1:]
            ya[scen, model] = Ya
            strace[scen, model] = sim["s"][:, 365:]
            tot = sim["totals"]
            p = partition(Ya * AREA[:, None])
            R = 0.8 * Ya.mean()
            rows.append(dict(scenario=scen, model=model, mean_yield=Ya.mean(), CV_L=p["CV_L"],
                             CV_R=p["CV_R"], phi=p["phi"],
                             p_below_80pct_of_mean=(Ya < R).mean(),
                             worst_year_regional=(Ya * AREA[:, None]).sum(0).min() / (Ya * AREA[:, None]).sum(0).mean(),
                             transp_share=tot["transp"] / tot["rain"], evap_share=tot["evap"] / tot["rain"],
                             drain_share=tot["drain"] / tot["rain"], runoff_share=tot["runoff"] / tot["rain"],
                             wb_error_mm=(tot["rain"] - tot["runoff"] - tot["transp"] - tot["evap"]
                                          - tot["drain"] - tot["ds"]) / N_FIELDS))
            print(f"{scen:18s} {model:9s} mean yield {Ya.mean():.2f} t/ha  CV_L {p['CV_L']:.3f}  "
                  f"phi {p['phi']:.2f}  T/rain {tot['transp'] / tot['rain']:.2f}  ({time.time() - t0:.0f} s)")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "richards_vs_bucket.csv"), index=False, float_format="%.4f")
    fy = N_FIELDS * N_YEARS
    print(f"\ntime per field-year: bucket {1e3 * timing['bucket'] / (2 * fy):.1f} ms, "
          f"Richards {1e3 * timing['Richards'] / (2 * fy):.0f} ms")

    # mixture advantage in each model
    for model in ("bucket", "Richards"):
        u, m = ya["Uniform improved", model], ya["Landrace mixture", model]
        print(f"{model:9s} mixture - uniform: mean yield {m.mean() - u.mean():+.2f} t/ha, "
              f"CV_L {partition(m * AREA[:, None])['CV_L'] - partition(u * AREA[:, None])['CV_L']:+.3f}")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
    cols = {"Uniform improved": "#2a78d6", "Landrace mixture": "#eb6834"}
    for scen, c in cols.items():
        ax[0].scatter(ya[scen, "bucket"].ravel(), ya[scen, "Richards"].ravel(), s=10, color=c, alpha=0.7, label=scen)
    lim = [0, max(v.max() for v in ya.values()) * 1.05]
    ax[0].plot(lim, lim, color="0.6", lw=0.8)
    ax[0].set(xlim=lim, ylim=lim, xlabel="bucket: annual yield (t/ha)", ylabel="Richards column: annual yield (t/ha)")
    ax[0].set_title("A. Field-years, same rain", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=8)

    d = np.arange(365, 365 + 3 * 365)
    f = 0
    ax[1].plot(d / 365, strace["Landrace mixture", "bucket"][f, d - 365], color="0.4", lw=0.8, label="bucket")
    ax[1].plot(d / 365, strace["Landrace mixture", "Richards"][f, d - 365], color="#13766E", lw=0.8,
               label=f"column, top {SOIL.z_r:.0f} mm")
    ax[1].set(xlabel="year", ylabel="relative soil moisture s")
    ax[1].set_title("B. Root-zone soil moisture, one farm", loc="left", fontsize=10)
    ax[1].legend(frameon=False, fontsize=8)

    x = np.arange(2)
    for j, model in enumerate(("bucket", "Richards")):
        v = [df[(df.scenario == s) & (df.model == model)].mean_yield.item() for s in cols]
        ax[2].bar(x + (j - 0.5) * 0.38, v, 0.38, color=["0.6", "#13766E"][j], label=model)
    ax[2].set_xticks(x, list(cols))
    ax[2].set_ylabel("mean annual yield (t/ha)")
    ax[2].set_title("C. Mean yield by landscape", loc="left", fontsize=10)
    ax[2].legend(frameon=False, fontsize=8)
    fig.suptitle("Illustrative, uncalibrated parameters", fontsize=8, color="0.4", x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "richards_vs_bucket.png"), dpi=150)


if __name__ == "__main__":
    main()
