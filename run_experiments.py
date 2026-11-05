"""
Scenario experiments for the tessera prototype. Run:  python run_experiments.py

Four experiments, each averaged over three rainfall ensembles:
  E1/E2  five landscapes (how varieties are spread across farms) at five rainfall
         correlation lengths L: yields, shortfall risk and synchrony
  E3     markets at five levels of market integration theta: prices and food gaps
  E4     how large the improved variety's yield advantage must be before uniform
         stands beat the landrace mixture on shortfall risk

Writes results.csv and two figures in figures/. All parameters are illustrative.
A full run takes about ten minutes.
"""
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mixtures_model import (VARIETIES, Soil, Climate, Disease, Market, make_fields, generate_rain,
                            simulate, field_yields, annual, partition, shortfall_stats, market)

# ---------------------------------------------------------------- set-up
N_FIELDS, N_YEARS, SEEDS = 64, 61, (7, 8, 9)      # 64 farms, 61 years (first is spin-up), 3 rainfall ensembles
SOIL, DIS = Soil(), Disease()
XY, AREA = make_fields(N_FIELDS)
NAMES = list(VARIETIES)                           # LR_early, LR_mid, LR_late, IMP

# Scenario colors in fixed categorical order
COLORS = {"Uniform improved": "#2a78d6", "Landrace mixture": "#eb6834", "Mosaic of landraces": "#1baf7a",
          "Half improved, half mixture": "#eda100", "Improved added to mixture": "#e87ba4"}


def scenarios(nf=N_FIELDS, seed=3):
    """The five landscapes, each an (n_fields, 4) array of planting shares in the order
    LR_early, LR_mid, LR_late, IMP:
      Uniform improved             every farm grows only the improved variety
      Landrace mixture             every farm grows the three landraces in equal shares
      Mosaic of landraces          every farm grows one landrace, chosen at random
                                   (diversity across farms, not within them)
      Half improved, half mixture  alternate farms grow the improved variety or the mixture
      Improved added to mixture    every farm grows all four in equal shares
    """
    rng = np.random.default_rng(seed)
    uni = np.tile([0, 0, 0, 1.0], (nf, 1))
    mix = np.tile([1 / 3, 1 / 3, 1 / 3, 0], (nf, 1))
    mosaic = np.zeros((nf, 4)); mosaic[np.arange(nf), rng.integers(0, 3, nf)] = 1.0
    half = np.where((np.arange(nf) % 2 == 0)[:, None], uni, mix)
    add = np.tile([0.25, 0.25, 0.25, 0.25], (nf, 1))
    return {"Uniform improved": uni, "Landrace mixture": mix, "Mosaic of landraces": mosaic,
            "Half improved, half mixture": half, "Improved added to mixture": add}


SCEN = scenarios()


def run(shares, L=100.0, seed=7, varieties=VARIETIES, keep_s=False):
    """Simulate one landscape for one rainfall ensemble. Returns yields (t/ha) and
    production per household (t), both (n_fields, n_years)."""
    rain, _ = generate_rain(XY, N_YEARS, Climate(L=L), seed=seed)
    sim = simulate(rain, shares, SOIL, varieties)
    Y = field_yields(sim, shares, DIS, seed=seed + 100)
    Ya = annual(Y)[:, 1:]                          # drop spin-up year, t/ha
    prod = Ya * AREA[:, None]                      # tonnes per household
    out = dict(Ya=Ya, prod=prod)
    if keep_s:
        out.update(rain=rain, s=sim["s"])
    return out


def metrics(prod, mk=Market()):
    """Summary statistics for one run: mean yield, CV_L, CV_R, synchrony phi,
    shortfall risk and run length, and the share of years with a regional deficit."""
    part = partition(prod)
    sf = shortfall_stats(prod, mk.need)
    reg = prod.sum(0)
    return dict(mean_yield_t_per_ha=(prod.sum(0) / AREA.sum()).mean(), **part, **sf,
                p_regional_deficit=(reg < mk.need * len(prod)).mean())


# ---------------------------------------------------------------- experiments
rows = []          # one row of results per experiment, landscape and setting
cache = {}         # runs at L = 100 km, reused by the market experiment
# E1 and E2: every landscape at every rainfall correlation length, averaged over ensembles
for L in (10, 30, 100, 300, 1000):
    for name, sh in SCEN.items():
        ms = []
        for sd in SEEDS:
            r = run(sh, L=L, seed=sd)
            if L == 100:
                cache[(name, sd)] = r
            ms.append(metrics(r["prod"]))
        m = {k: np.mean([x[k] for x in ms]) for k in ms[0]}
        rows.append(dict(experiment="scenarios_by_L", scenario=name, L_km=L, **m))
    print("done L", L)

# E3: markets across integration theta (L = 100 km). Net buyers are households that
# on average harvest less than they need.
for theta in (0.0, 0.25, 0.5, 0.75, 1.0):
    mk = Market(theta=theta)
    for name in SCEN:
        vals = []
        for sd in SEEDS:
            r = cache[(name, sd)]
            o = market(r["prod"], AREA, mk, seed=sd + 200)
            buyers = ~o["seller"]
            vals.append(dict(price_cv=o["price"].std() / o["price"].mean(),
                             price_max_ratio=o["price"].max() / mk.p_bar,
                             gap_buyers=o["gap"][buyers].mean(),
                             gap_sellers=o["gap"][~buyers].mean(),
                             welfare_sd_sellers=o["welfare"][~buyers].std(),
                             share_buyers=buyers.mean()))
        m = {k: np.mean([v[k] for v in vals]) for k in vals[0]}
        rows.append(dict(experiment="markets_by_theta", scenario=name, theta=theta, L_km=100, **m))
print("done markets")

# E4: break-even yield potential of the improved variety. The improved variety's y_max is
# varied from 1.6 to 3.0 t/ha per season while the landraces stay fixed.
for ymax in np.arange(1.6, 3.01, 0.2):
    v = copy.deepcopy(VARIETIES); v["IMP"].y_max = float(ymax)
    for name in ("Uniform improved", "Landrace mixture"):
        ms = [metrics(run(SCEN[name], L=100, seed=sd, varieties=v)["prod"]) for sd in SEEDS]
        m = {k: np.mean([x[k] for x in ms]) for k in ms[0]}
        rows.append(dict(experiment="breakeven_ymax", scenario=name, imp_ymax=round(ymax, 1), L_km=100, **m))
print("done break-even")

df = pd.DataFrame(rows)
df.to_csv("results.csv", index=False)

# ---------------------------------------------------------------- figures
from make_figures import make_figures
make_figures(run, SCEN, COLORS, df)
