"""
Prototype model: crop mixtures, yield stability and regional food supply.

Implements appendix sections A1 to A3 of "Crop Mixtures and Regional Food Supply":
  A1  daily stochastic soil water balance per field, shared by all varieties in the
      field, with Poisson storms whose rate varies by season, ENSO and place;
      stage-based yield response (Doorenbos and Kassam); correlated pest and
      disease losses that depend on shared varieties.
  A2  Wang and Loreau (2014) partition of regional variability:
      CV_R^2 = phi * CV_L^2.
  A3  grain storage bucket, local price with market integration theta,
      household net position and food gap.

All parameter values are ILLUSTRATIVE (semi-arid East Africa, sorghum-like
varieties). Nothing here is calibrated. Units: mm, days, t/ha, USD/t.
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.stats import norm

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

@dataclass
class Variety:
    name: str
    maturity: int        # days from planting to harvest
    y_max: float         # t/ha with no water stress
    s_star: float        # soil moisture where transpiration starts to fall
    s_w: float           # wilting point
    t_max: float         # mm/day maximum transpiration at full cover
    ky: tuple = (0.20, 0.55, 0.45)          # yield response: vegetative, flowering, grain fill
    stage_bounds: tuple = (0.45, 0.65)      # stage boundaries as fraction of season


VARIETIES = {
    "LR_early": Variety("LR_early", 85, 1.4, 0.40, 0.18, 3.5),
    "LR_mid":   Variety("LR_mid", 105, 1.7, 0.45, 0.20, 4.0),
    "LR_late":  Variety("LR_late", 125, 2.0, 0.50, 0.22, 4.5),
    "IMP":      Variety("IMP", 100, 2.8, 0.50, 0.22, 5.0),
}


@dataclass
class Soil:
    n: float = 0.43       # porosity
    z_r: float = 500.0    # root zone depth, mm
    k_s: float = 40.0     # saturated drainage rate, mm/day
    c: float = 11.0       # drainage exponent
    e_max: float = 1.0    # maximum soil evaporation, mm/day

    @property
    def nzr(self):
        return self.n * self.z_r


@dataclass
class Climate:
    # daily storm probability outside / inside the two rainy seasons
    p_base: float = 0.04
    p_long: float = 0.30      # long rains, days 60-151
    p_short: float = 0.28     # short rains, days 274-365
    h: float = 9.0            # mean storm depth, mm
    enso_ar: float = 0.5      # ENSO lag-1 autocorrelation (annual)
    load_long: float = -0.10  # ENSO loading on long-rains storm rate
    load_short: float = 0.35  # ENSO loading on short-rains storm rate
    sigma_eta: float = 0.35   # spread of seasonal storm-rate anomalies
    L: float = 100.0          # km, correlation length of seasonal anomalies
    L_day: float = 20.0       # km, correlation length of daily storm occurrence


@dataclass
class Disease:
    d0: float = 0.08          # mean loss fraction
    sd: float = 0.06
    rho: float = 0.6          # correlation scale for fields with identical varieties
    d_max: float = 0.6


LONG_RAINS = (60, 152)
SHORT_RAINS = (274, 366)
PLANT_DAYS = (70, 284)        # planting day of year for each season


# --------------------------------------------------------------------------
# Landscape and rainfall
# --------------------------------------------------------------------------

def make_fields(n_fields=64, extent_km=150.0, seed=1):
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, extent_km, size=(n_fields, 2))
    area = rng.uniform(0.5, 1.5, size=n_fields)       # ha, so some households are net buyers
    return xy, area


def _chol_exp(xy, length):
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    cov = np.exp(-d / length) + 1e-9 * np.eye(len(xy))
    return np.linalg.cholesky(cov)


def generate_rain(xy, n_years, clim: Climate, seed=7):
    """Daily rainfall (n_fields, n_days). The same seed gives the same underlying
    random numbers for any L, so runs at different L are directly comparable."""
    rng = np.random.default_rng(seed)
    nf = len(xy)
    n_days = 365 * n_years
    doy = np.arange(n_days) % 365
    year = np.arange(n_days) // 365

    # ENSO index, AR(1), unit variance
    e = np.zeros(n_years)
    eps = rng.standard_normal(n_years)
    e[0] = eps[0]
    for y in range(1, n_years):
        e[y] = clim.enso_ar * e[y - 1] + np.sqrt(1 - clim.enso_ar ** 2) * eps[y]

    # seasonal anomalies, spatially correlated with length L
    chol_L = _chol_exp(xy, clim.L)
    eta_long = clim.sigma_eta * (chol_L @ rng.standard_normal((nf, n_years)))
    eta_short = clim.sigma_eta * (chol_L @ rng.standard_normal((nf, n_years)))

    base = np.full(n_days, clim.p_base)
    in_long = (doy >= LONG_RAINS[0]) & (doy < LONG_RAINS[1])
    in_short = (doy >= SHORT_RAINS[0]) & (doy < SHORT_RAINS[1])
    p = np.tile(base, (nf, 1))
    p_long = clim.p_long * np.exp(clim.load_long * e[None, :] + eta_long)    # (nf, years)
    p_short = clim.p_short * np.exp(clim.load_short * e[None, :] + eta_short)
    p[:, in_long] = p_long[:, year[in_long]]
    p[:, in_short] = p_short[:, year[in_short]]
    p = np.clip(p, 0.0, 0.95)

    # daily occurrence: Gaussian copula with short correlation length
    chol_d = _chol_exp(xy, clim.L_day)
    z = chol_d @ rng.standard_normal((nf, n_days))
    occur = z < norm.ppf(p)
    depth = rng.exponential(clim.h, size=(nf, n_days))
    return np.where(occur, depth, 0.0), e


# --------------------------------------------------------------------------
# A1: soil water balance and yields
# --------------------------------------------------------------------------

def _crop_calendar(var: Variety, n_years):
    """For each day: crop index (or -1), stage (0,1,2) and cover factor kc."""
    n_days = 365 * n_years
    crop = np.full(n_days, -1)
    stage = np.zeros(n_days, dtype=int)
    kc = np.zeros(n_days)
    k = 0
    for y in range(n_years):
        for pd in PLANT_DAYS:
            start = 365 * y + pd
            end = start + var.maturity
            if end > n_days:
                continue
            idx = np.arange(start, end)
            f = (idx - start) / var.maturity
            crop[idx] = k
            stage[idx] = np.where(f < var.stage_bounds[0], 0, np.where(f < var.stage_bounds[1], 1, 2))
            # canopy cover: ramp up to 40 percent of season, decline after 80 percent
            kc[idx] = np.clip(np.minimum(0.3 + 0.7 * f / 0.4, 1.0 - 0.4 * np.clip((f - 0.8) / 0.2, 0, 1)), 0, 1)
            k += 1
    return crop, stage, kc, k


def simulate(rain, shares, soil: Soil, varieties=VARIETIES, s0=0.3):
    """Run the daily water balance.

    rain   : (n_fields, n_days) mm
    shares : (n_fields, n_varieties) planting shares c_v, rows sum to 1
    Returns per-variety per-season yields before disease, soil moisture,
    and water-balance totals.
    """
    names = list(varieties)
    vs = [varieties[k] for k in names]
    nf, n_days = rain.shape
    n_years = n_days // 365
    nv = len(vs)

    cal = [_crop_calendar(v, n_years) for v in vs]
    n_crops = min(c[3] for c in cal)
    s_star = np.array([v.s_star for v in vs])
    s_w = np.array([v.s_w for v in vs])
    t_max = np.array([v.t_max for v in vs])

    t_act = np.zeros((nf, nv, n_crops, 3))
    t_pot = np.zeros((nf, nv, n_crops, 3))
    s = np.full(nf, s0)
    s_trace = np.zeros((nf, n_days), dtype=np.float32)
    totals = dict(rain=rain.sum(), runoff=0.0, transp=0.0, evap=0.0, drain=0.0, ds=0.0)
    nzr = soil.nzr

    crop_idx = np.stack([c[0] for c in cal])     # (nv, n_days)
    stage_idx = np.stack([c[1] for c in cal])
    kc_all = np.stack([c[2] for c in cal])

    s_start = s.copy()
    for t in range(n_days):
        # rain pulse with saturation excess runoff
        room = nzr * (1.0 - s)
        infil = np.minimum(rain[:, t], room)
        totals["runoff"] += (rain[:, t] - infil).sum()
        s = s + infil / nzr

        # transpiration by variety (shared bucket)
        f = np.clip((s[:, None] - s_w[None, :]) / (s_star - s_w)[None, :], 0.0, 1.0)   # (nf, nv)
        active = crop_idx[:, t] >= 0
        pot = shares * (t_max * kc_all[:, t] * active)[None, :]                       # (nf, nv)
        act = pot * f
        evap = soil.e_max * s
        drain = soil.k_s * s ** soil.c
        loss = act.sum(1) + evap + drain
        # do not remove more water than is there
        scale = np.minimum(1.0, s * nzr / np.maximum(loss, 1e-12))
        act *= scale[:, None]
        evap *= scale
        drain *= scale
        s = s - (act.sum(1) + evap + drain) / nzr
        s = np.clip(s, 0.0, 1.0)

        totals["transp"] += act.sum()
        totals["evap"] += evap.sum()
        totals["drain"] += drain.sum()
        for j in np.nonzero(active)[0]:
            k = crop_idx[j, t]
            if k < n_crops:
                st = stage_idx[j, t]
                t_act[:, j, k, st] += act[:, j]
                t_pot[:, j, k, st] += pot[:, j]
        s_trace[:, t] = s

    totals["ds"] = ((s - s_start) * nzr).sum()

    # Doorenbos and Kassam yield by stage, per unit area of each variety
    ky = np.array([v.ky for v in vs])                         # (nv, 3)
    y_max = np.array([v.y_max for v in vs])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(t_pot > 0, t_act / t_pot, 1.0)
    factor = np.clip(1.0 - ky[None, :, None, :] * (1.0 - ratio), 0.0, 1.0).prod(-1)   # (nf, nv, crops)
    y_var = y_max[None, :, None] * factor
    return dict(y_var=y_var, s=s_trace, totals=totals, names=names, n_crops=n_crops)


def jaccard(shares):
    present = shares > 0
    inter = present.astype(float) @ present.T.astype(float)
    union = present.sum(1)[:, None] + present.sum(1)[None, :] - inter
    return inter / union


def disease_losses(shares, n_crops, dis: Disease, seed=11):
    rng = np.random.default_rng(seed)
    J = jaccard(shares)
    nf = len(J)
    cov = dis.rho * J + (1 - dis.rho) * np.eye(nf)
    chol = np.linalg.cholesky(cov + 1e-9 * np.eye(nf))
    z = chol @ rng.standard_normal((nf, n_crops))
    return np.clip(dis.d0 + dis.sd * z, 0.0, dis.d_max)


def field_yields(sim, shares, dis: Disease, seed=11):
    """Season yields Y_{i,t} = (1 - D) * sum_v c_v Y_v  (t/ha), shape (n_fields, n_crops)."""
    y_mix = (shares[:, :, None] * sim["y_var"]).sum(1)
    D = disease_losses(shares, sim["n_crops"], dis, seed)
    return (1.0 - D) * y_mix


def annual(season_yields):
    """Sum the two seasons planted in each year: (n_fields, n_years)."""
    nf, nc = season_yields.shape
    ny = nc // 2
    return season_yields[:, : 2 * ny].reshape(nf, ny, 2).sum(-1)


# --------------------------------------------------------------------------
# A2: regional synchrony
# --------------------------------------------------------------------------

def partition(production):
    """production: (n_fields, n_years) in tonnes. Returns CV_L, CV_R, phi."""
    m = production.mean(1)
    cov = np.cov(production)
    sd = np.sqrt(np.diag(cov))
    cv_l = sd.sum() / m.sum()
    cv_r = np.sqrt(cov.sum()) / m.sum()
    phi = cov.sum() / sd.sum() ** 2
    return dict(CV_L=cv_l, CV_R=cv_r, phi=phi)


def shortfall_stats(yields_annual, R):
    short = yields_annual < R
    p = short.mean()
    S = np.maximum(R - yields_annual, 0).mean()
    # mean length of runs of consecutive shortfall years
    runs = []
    for row in short:
        n = 0
        for v in row:
            if v:
                n += 1
            elif n:
                runs.append(n)
                n = 0
        if n:
            runs.append(n)
    return dict(p_short=p, exp_short=S, mean_run=np.mean(runs) if runs else 0.0)


# --------------------------------------------------------------------------
# A3: storage, prices, households
# --------------------------------------------------------------------------

@dataclass
class Market:
    theta: float = 0.2        # market integration, 0 isolated, 1 fully integrated
    eps: float = 0.3          # price elasticity of demand (absolute)
    p_bar: float = 250.0      # normal local price, USD/t
    sigma_w: float = 0.15     # variability of import parity price (log)
    delta: float = 0.10       # storage loss per year
    s_max_frac: float = 0.25  # storage capacity as share of normal consumption
    need: float = 1.8         # household grain need, t/yr
    cash: float = 150.0       # household cash available for grain, USD/yr
    c_bar: float = None       # normal regional demand, t/yr (None: mean regional harvest)


def market(production, area, mk: Market, seed=21):
    """production: (n_fields, n_years) tonnes. Returns price series and household outcomes."""
    rng = np.random.default_rng(seed)
    nf, ny = production.shape
    # Normal regional demand is set to the landscape's mean harvest, so that p_bar is
    # the price in a normal year and prices respond to departures from normal.
    c_bar = production.sum(0).mean() if mk.c_bar is None else mk.c_bar
    s_max = mk.s_max_frac * c_bar
    p_w = mk.p_bar * np.exp(mk.sigma_w * rng.standard_normal(ny))
    P = production.sum(0)
    S = 0.0
    Q = np.zeros(ny)
    for y in range(ny):
        A = P[y] + (1 - mk.delta) * S
        S = min(s_max, max(A - c_bar, 0.0))
        Q[y] = A - S
    lnp = (1 - mk.theta) * (np.log(mk.p_bar) - np.log(Q / c_bar) / mk.eps) + mk.theta * np.log(p_w)
    price = np.exp(lnp)

    net = production - mk.need                              # tonnes, per household
    seller = net.mean(1) > 0
    gap = np.maximum(mk.need - production - mk.cash / price[None, :], 0.0)
    welfare = net * (price[None, :] - mk.p_bar)             # first-order, USD
    return dict(price=price, Q=Q, P=P, seller=seller, gap=gap, welfare=welfare)
