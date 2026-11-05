"""
Scenario experiments for the mixtures prototype. Run:  python run_experiments.py
Writes results.csv and two figures in figures/. All parameters are illustrative.
"""
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mixtures_model import (VARIETIES, Soil, Climate, Disease, Market, make_fields, generate_rain,
                            simulate, field_yields, annual, partition, shortfall_stats, market)

N_FIELDS, N_YEARS, SEEDS = 64, 61, (7, 8, 9)      # first year is spin-up
SOIL, DIS = Soil(), Disease()
XY, AREA = make_fields(N_FIELDS)
NAMES = list(VARIETIES)                           # LR_early, LR_mid, LR_late, IMP

# Scenario colors in fixed categorical order
COLORS = {"Uniform improved": "#2a78d6", "Landrace mixture": "#eb6834", "Mosaic of landraces": "#1baf7a",
          "Half improved, half mixture": "#eda100", "Improved added to mixture": "#e87ba4"}


def scenarios(nf=N_FIELDS, seed=3):
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
    part = partition(prod)
    sf = shortfall_stats(prod, mk.need)
    reg = prod.sum(0)
    return dict(mean_yield_t_per_ha=(prod.sum(0) / AREA.sum()).mean(), **part, **sf,
                p_regional_deficit=(reg < mk.need * len(prod)).mean())


rows = []
cache = {}
# E1 and E2: scenarios across rainfall correlation lengths
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

# E3: markets across integration theta (L = 100 km)
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

# E4: break-even yield potential of the improved variety
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
