"""
Richards-equation soil water for tessera, using SoilWaterFlow
(https://github.com/mjpuma/SoilWaterFlow).

`simulate_richards` is a drop-in alternative to `mixtures_model.simulate`: it takes the
same rainfall, planting shares, soil and varieties and returns the same dictionary, so
`field_yields`, `annual`, `partition` and `market` work unchanged. The difference is the
soil: instead of one well-mixed root zone shared by all varieties, each field is a soil
column (2 m by default) in which water infiltrates, redistributes and drains, and each
variety takes up water through its own root profile. Varieties then compete only for
the water their roots can reach.

It takes about 0.6 s per field-year on one core, against a few milliseconds for the
bucket, so it is meant for benchmarks and for a few calibrated sites rather than the
full scenario runs. See validation/richards_vs_bucket.py.

Root traits are new parameters: tessera's varieties have none. The defaults in ROOTS
are ILLUSTRATIVE and uncalibrated, like every other parameter in tessera.
"""
from dataclasses import dataclass

import numpy as np

from soilwaterflow.model import Soil as ColumnSoil, Variety as ColumnVariety, simulate_fields

from mixtures_model import VARIETIES, Soil, _crop_calendar, stage_yields


@dataclass
class Roots:
    depth: float                 # rooting depth, mm
    shape: tuple = (1.0, 2.0)    # beta(a, b) root density over depth / rooting depth; a < b: shallow-heavy


# Illustrative: longer-season landraces root deeper; the improved variety is shallower.
ROOTS = {
    "LR_early": Roots(900.0, (1.0, 2.2)),
    "LR_mid": Roots(1100.0, (1.1, 1.9)),
    "LR_late": Roots(1300.0, (1.2, 1.7)),
    "IMP": Roots(800.0, (1.0, 2.5)),
}


def column_soil(soil: Soil, air_entry=150.0):
    """A Brooks-Corey soil matching tessera's bucket soil where the two models overlap:
    porosity n, saturated conductivity K_s, and Brooks-Corey exponent b from the leakage
    exponent beta = 2b + 4 (Laio et al. 2001). Bare-soil evaporation E_b acts over the top
    150 mm, at its full rate above field capacity and falling linearly to zero in dry soil.
    The air-entry head is not a bucket parameter; it is illustrative."""
    b = (soil.beta - 4.0) / 2.0
    cs = ColumnSoil(porosity=soil.n, air_entry=air_entry, b=b, k_sat=soil.k_s, s_res=0.0,
                    e_soil=soil.e_bare, z_evap=150.0)
    cs.h_evap = float(cs.head(soil.s_fc))
    return cs


def column_varieties(varieties, soil: Soil, cs: ColumnSoil, roots=ROOTS):
    """Map tessera varieties to SoilWaterFlow varieties. Potential transpiration is
    E_max - E_w, as in the bucket; s* and s_w are converted to heads with the soil's
    retention curve."""
    out = []
    for name, v in varieties.items():
        r = roots[name]
        out.append(ColumnVariety.from_saturations(
            name, t_max=v.e_max - soil.e_w, s_star=v.s_star, s_w=v.s_w, soil=cs,
            root_depth=r.depth, root_shape=tuple(r.shape)))
    return out


def simulate_richards(rain, shares, soil: Soil = None, varieties=VARIETIES, roots=ROOTS,
                      s0=0.3, depth=2000.0, dz=25.0, dt_max=0.05, storm_hours=4.0, n_workers=None):
    """Run the soil column for every field and compute yields by variety.

    rain, shares, soil, varieties : as in mixtures_model.simulate
    roots       : rooting depth and shape of each variety (ROOTS by default)
    depth, dz   : soil column depth and layer thickness, mm
    storm_hours : each day's rain falls at a constant rate over this many hours
                  (the bucket applies it as a single pulse)
    Returns the same dictionary as mixtures_model.simulate. `s` is the mean saturation of
    the top soil.z_r mm, the bucket's root zone, so the two can be compared directly.
    Extra key `results`: the SoilWaterFlow result of each field (daily fluxes, uptake by depth).
    """
    soil = soil or Soil()
    names = list(varieties)
    vs = [varieties[k] for k in names]
    nf, n_days = rain.shape
    n_years = n_days // 365
    nv = len(vs)

    cal = [_crop_calendar(v, n_years) for v in vs]
    n_crops = min(c[3] for c in cal)
    crop_idx = np.stack([c[0] for c in cal])     # (nv, n_days)
    stage_idx = np.stack([c[1] for c in cal])
    kc_all = np.stack([c[2] for c in cal])
    active = crop_idx >= 0

    # canopy cover of each variety on each field and day: planting share x cover x growing
    cover = shares[:, None, :] * (kc_all * active).T[None, :, :]     # (nf, n_days, nv)

    cs = column_soil(soil)
    cvs = column_varieties(varieties, soil, cs, roots)
    res = simulate_fields(cs, cvs, rain, cover=cover, n_workers=n_workers, depth=depth, dz=dz,
                          s0=s0, dt_max=dt_max, storm_hours=storm_hours, report_depth=soil.z_r)

    # transpiration by field, variety, crop and growth stage, as in mixtures_model.simulate
    t_act = np.zeros((nf, nv, n_crops, 3))
    t_pot = np.zeros((nf, nv, n_crops, 3))
    act = np.stack([r.transpiration for r in res])         # (nf, n_days, nv)
    pot = np.stack([r.transpiration_pot for r in res])
    for j in range(nv):
        days = np.nonzero(active[j] & (crop_idx[j] < n_crops))[0]
        k, st = crop_idx[j, days], stage_idx[j, days]
        for f in range(nf):
            np.add.at(t_act[f, j], (k, st), act[f, days, j])
            np.add.at(t_pot[f, j], (k, st), pot[f, days, j])

    totals = dict(rain=rain.sum(), runoff=sum(r.runoff.sum() for r in res),
                  transp=float(act.sum()), evap=sum(r.evaporation.sum() for r in res),
                  drain=sum(r.drainage.sum() for r in res),
                  ds=sum(r.storage[-1] - r.storage0 for r in res))
    s_trace = np.stack([r.s_top for r in res]).astype(np.float32)
    return dict(y_var=stage_yields(t_act, t_pot, vs), s=s_trace, totals=totals, names=names,
                n_crops=n_crops, results=res)
