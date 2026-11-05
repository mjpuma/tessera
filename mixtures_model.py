"""
tessera model: crop mixtures, yield stability and regional food supply.

The model follows the effect of crop mixtures from single fields up to a region and its
grain market. It has three parts, matching appendix sections A1 to A3 of
"Crop Mixtures and Regional Food Supply":

  A1  Soil water and yields. Each field has a daily stochastic soil water balance
      (soilwater.py; the loss function of Laio et al. 2001, weighted by the canopy cover
      of each variety). All varieties in a field share the same soil water. Storms arrive
      as a Poisson process whose rate varies by season, with ENSO and from place to place.
      Yields fall with water stress by growth stage (Doorenbos and Kassam 1979), and pest
      and disease losses are more correlated between fields that grow the same varieties.
  A2  From fields to the region. Regional variability is split into the variability of
      individual farms and their synchrony (Wang and Loreau 2014): CV_R^2 = phi * CV_L^2.
  A3  Markets, storage and households. A regional grain store, a local price that depends
      on local supply and outside prices (market integration theta), and each household's
      food gap.

A typical run (see run_experiments.py):

    xy, area = make_fields()                       # where the farms are and how big
    rain, enso = generate_rain(xy, n_years, Climate())
    sim = simulate(rain, shares, Soil())           # soil water and yields by variety
    Y = annual(field_yields(sim, shares, Disease()))   # t/ha per farm per year
    partition(Y * area[:, None])                   # CV_L, CV_R and synchrony phi
    market(Y * area[:, None], area, Market())      # prices and household food gaps

All parameter values are ILLUSTRATIVE (semi-arid East Africa, sorghum-like varieties).
Nothing here is calibrated. Units: mm, days, t/ha, USD/t.
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.stats import norm, poisson
from soilwater import SoilParams, infiltrate, drain

# ==========================================================================
# Parameters
# ==========================================================================

@dataclass
class Variety:
    """Traits of one crop variety. Soil moisture values are relative (0 to 1)."""
    name: str
    maturity: int        # days from planting to harvest
    y_max: float         # yield with no water stress, t/ha per season
    s_star: float        # soil moisture below which transpiration starts to fall (stomata close)
    s_w: float           # wilting point: transpiration stops
    e_max: float         # maximum evapotranspiration at full cover, mm/day (Laio et al. E_max)
    ky: tuple = (0.20, 0.55, 0.45)          # yield response factors: vegetative, flowering, grain fill
    stage_bounds: tuple = (0.45, 0.65)      # where each stage ends, as a fraction of the season


# Three landraces with different maturity and one improved variety. The early landrace
# tolerates drought better (lower s_star and s_w); the improved variety yields more when
# water is plentiful but is less tolerant.
VARIETIES = {
    #                    name       days  y_max  s_star  s_w   e_max
    "LR_early": Variety("LR_early", 85,   1.4,   0.40,   0.18, 3.5),
    "LR_mid":   Variety("LR_mid",   105,  1.7,   0.45,   0.20, 4.0),
    "LR_late":  Variety("LR_late",  125,  2.0,   0.50,   0.22, 4.5),
    "IMP":      Variety("IMP",      100,  2.8,   0.50,   0.22, 5.0),
}


# Soil parameters live in soilwater.py. The defaults there are a sandy loam after
# Laio et al. (2001), Table 1, converted to mm and days.
Soil = SoilParams


@dataclass
class Climate:
    """Rainfall generator settings, loosely based on bimodal semi-arid East Africa."""
    # storm arrival rate (storms per day) outside / inside the two rainy seasons
    p_base: float = 0.04
    p_long: float = 0.30      # long rains, days 60-151
    p_short: float = 0.28     # short rains, days 274-365
    h: float = 9.0            # mean storm depth, mm
    enso_ar: float = 0.5      # ENSO lag-1 autocorrelation (annual)
    load_long: float = -0.10  # ENSO loading on long-rains storm rate (El Nino slightly drier)
    load_short: float = 0.35  # ENSO loading on short-rains storm rate (El Nino wetter)
    sigma_eta: float = 0.35   # spread of seasonal storm-rate anomalies between places and years
    L: float = 100.0          # km, correlation length of seasonal anomalies (key parameter)
    L_day: float = 20.0       # km, correlation length of individual storms


@dataclass
class Disease:
    """Pest and disease losses as a fraction of yield."""
    d0: float = 0.08          # mean loss fraction
    sd: float = 0.06          # spread of losses
    rho: float = 0.6          # correlation between two fields with identical varieties
    d_max: float = 0.6        # largest possible loss


# Calendar: day-of-year windows for the two rainy seasons and the planting day of each
LONG_RAINS = (60, 152)
SHORT_RAINS = (274, 366)
PLANT_DAYS = (70, 284)


# ==========================================================================
# Landscape and rainfall
# ==========================================================================

def make_fields(n_fields=64, extent_km=150.0, seed=1):
    """Place farms at random in a square region. Returns coordinates (km) and farm areas (ha).
    Farm areas vary from 0.5 to 1.5 ha so that some households grow less than they need."""
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, extent_km, size=(n_fields, 2))
    area = rng.uniform(0.5, 1.5, size=n_fields)
    return xy, area


def _chol_exp(xy, length):
    """Cholesky factor of an exponential correlation matrix, corr = exp(-distance / length).
    Multiplying independent normal draws by it gives spatially correlated draws."""
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    cov = np.exp(-d / length) + 1e-9 * np.eye(len(xy))     # tiny jitter keeps it positive definite
    return np.linalg.cholesky(cov)


def generate_rain(xy, n_years, clim: Climate, seed=7):
    """Daily rainfall for every farm: array (n_fields, n_days) in mm, and the ENSO index.

    Rainfall is built in three layers:
      1. an ENSO index for each year;
      2. a storm rate for each farm, season and year, which depends on the season, ENSO
         and a random anomaly that is spatially correlated over length L;
      3. daily storm counts and depths, correlated between nearby farms over L_day.
    The same seed gives the same underlying random numbers for any L, so runs at
    different L are directly comparable."""
    rng = np.random.default_rng(seed)
    nf = len(xy)
    n_days = 365 * n_years
    doy = np.arange(n_days) % 365          # day of year for each simulation day
    year = np.arange(n_days) // 365        # year index for each simulation day

    # 1. ENSO index: first-order autoregressive, unit variance
    e = np.zeros(n_years)
    eps = rng.standard_normal(n_years)
    e[0] = eps[0]
    for y in range(1, n_years):
        e[y] = clim.enso_ar * e[y - 1] + np.sqrt(1 - clim.enso_ar ** 2) * eps[y]

    # 2a. Seasonal anomalies for each farm and year, correlated across farms over length L
    chol_L = _chol_exp(xy, clim.L)
    eta_long = clim.sigma_eta * (chol_L @ rng.standard_normal((nf, n_years)))
    eta_short = clim.sigma_eta * (chol_L @ rng.standard_normal((nf, n_years)))

    # 2b. Storm rate p (storms per day) for each farm and day:
    #     base rate outside the rains; in the rains, seasonal rate x exp(ENSO + anomaly)
    base = np.full(n_days, clim.p_base)
    in_long = (doy >= LONG_RAINS[0]) & (doy < LONG_RAINS[1])
    in_short = (doy >= SHORT_RAINS[0]) & (doy < SHORT_RAINS[1])
    p = np.tile(base, (nf, 1))
    p_long = clim.p_long * np.exp(clim.load_long * e[None, :] + eta_long)    # (nf, years)
    p_short = clim.p_short * np.exp(clim.load_short * e[None, :] + eta_short)
    p[:, in_long] = p_long[:, year[in_long]]
    p[:, in_short] = p_short[:, year[in_short]]
    p = np.clip(p, 0.0, 2.0)

    # 3. Daily storm counts: Poisson with rate p. Correlated normal draws are turned into
    #    Poisson counts (a Gaussian copula), so nearby farms tend to get rain on the same days.
    #    Each storm's depth is exponential with mean h, so the day's total is a gamma draw.
    #    This is the daily total of the marked Poisson rainfall of Laio et al. (2001).
    chol_d = _chol_exp(xy, clim.L_day)
    z = chol_d @ rng.standard_normal((nf, n_days))
    count = poisson.ppf(norm.cdf(z), p)
    depth = rng.gamma(np.maximum(count, 1.0), clim.h)
    return np.where(count > 0, depth, 0.0), e


# ==========================================================================
# A1: soil water balance and yields
# ==========================================================================

def _crop_calendar(var: Variety, n_years):
    """Daily calendar for one variety, planted at the start of both seasons every year.

    Returns, for each simulation day: which crop (season) is growing (-1 if none),
    its growth stage (0 vegetative, 1 flowering, 2 grain fill), its canopy cover kc
    (0 to 1), and the total number of crops."""
    n_days = 365 * n_years
    crop = np.full(n_days, -1)
    stage = np.zeros(n_days, dtype=int)
    kc = np.zeros(n_days)
    k = 0
    for y in range(n_years):
        for pd in PLANT_DAYS:
            start = 365 * y + pd
            end = start + var.maturity
            if end > n_days:              # skip a crop that would not finish before the run ends
                continue
            idx = np.arange(start, end)
            f = (idx - start) / var.maturity          # fraction of the season elapsed
            crop[idx] = k
            stage[idx] = np.where(f < var.stage_bounds[0], 0, np.where(f < var.stage_bounds[1], 1, 2))
            # canopy cover: rises from 0.3 to 1 over the first 40 percent of the season,
            # stays at 1, then falls to 0.6 over the last 20 percent
            kc[idx] = np.clip(np.minimum(0.3 + 0.7 * f / 0.4, 1.0 - 0.4 * np.clip((f - 0.8) / 0.2, 0, 1)), 0, 1)
            k += 1
    return crop, stage, kc, k


def simulate(rain, shares, soil: Soil, varieties=VARIETIES, s0=0.3, max_ds=0.02):
    """Run the daily soil water balance for every field and compute yields by variety.

    rain   : (n_fields, n_days) mm, applied as a pulse at the start of each day
    shares : (n_fields, n_varieties) planting share of each variety, rows sum to 1
             (for example [0, 0, 0, 1] for a uniform improved stand)
    Returns a dict with:
      y_var  : (n_fields, n_varieties, n_crops) yield of each variety per hectare it occupies,
               before pest and disease losses
      s      : (n_fields, n_days) soil moisture at the end of each day
      totals : water balance totals over all fields (mm), for checking that water is conserved
    """
    names = list(varieties)
    vs = [varieties[k] for k in names]
    nf, n_days = rain.shape
    n_years = n_days // 365
    nv = len(vs)

    # Crop calendars and variety traits as arrays, one entry per variety
    cal = [_crop_calendar(v, n_years) for v in vs]
    n_crops = min(c[3] for c in cal)
    s_star = np.array([v.s_star for v in vs])
    s_w = np.array([v.s_w for v in vs])
    e_max = np.array([v.e_max for v in vs])

    # Accumulators: actual and unstressed (potential) transpiration, summed by
    # field, variety, crop and growth stage. Their ratio drives yield.
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
        # 1. Rain enters the soil; anything beyond saturation runs off
        s, runoff = infiltrate(s, rain[:, t], soil)
        totals["runoff"] += runoff.sum()

        # 2. Canopy cover of each variety today: planting share x cover x (growing or not)
        active = crop_idx[:, t] >= 0
        cov = shares * (kc_all[:, t] * active)[None, :]                        # (nf, nv)

        # 3. One day of losses (transpiration, evaporation, leakage) from the shared soil water
        s, act, evap, leak = drain(s, 1.0, cov, e_max, s_star, s_w, soil, max_ds=max_ds)
        pot = cov * (e_max - soil.e_w)          # what each variety would transpire with no stress

        # 4. Book-keeping: water balance totals and transpiration by crop and stage
        totals["transp"] += act.sum()
        totals["evap"] += evap.sum()
        totals["drain"] += leak.sum()
        for j in np.nonzero(active)[0]:
            k = crop_idx[j, t]
            if k < n_crops:
                st = stage_idx[j, t]
                t_act[:, j, k, st] += act[:, j]
                t_pot[:, j, k, st] += pot[:, j]
        s_trace[:, t] = s

    totals["ds"] = ((s - s_start) * nzr).sum()       # change in stored water

    # Yield (Doorenbos and Kassam 1979): in each stage, yield falls by K_y times the
    # relative transpiration deficit (1 - actual / potential); stage effects multiply.
    ky = np.array([v.ky for v in vs])                         # (nv, 3)
    y_max = np.array([v.y_max for v in vs])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(t_pot > 0, t_act / t_pot, 1.0)
    factor = np.clip(1.0 - ky[None, :, None, :] * (1.0 - ratio), 0.0, 1.0).prod(-1)   # (nf, nv, crops)
    y_var = y_max[None, :, None] * factor
    return dict(y_var=y_var, s=s_trace, totals=totals, names=names, n_crops=n_crops)


def jaccard(shares):
    """Share of varieties two fields have in common (intersection over union), for every
    pair of fields. 1 means identical variety sets, 0 means nothing in common."""
    present = shares > 0
    inter = present.astype(float) @ present.T.astype(float)
    union = present.sum(1)[:, None] + present.sum(1)[None, :] - inter
    return inter / union


def disease_losses(shares, n_crops, dis: Disease, seed=11):
    """Pest and disease loss fraction for each field and season. Losses in two fields are
    correlated in proportion to the varieties they share: corr = rho x Jaccard similarity."""
    rng = np.random.default_rng(seed)
    J = jaccard(shares)
    nf = len(J)
    cov = dis.rho * J + (1 - dis.rho) * np.eye(nf)
    chol = np.linalg.cholesky(cov + 1e-9 * np.eye(nf))
    z = chol @ rng.standard_normal((nf, n_crops))
    return np.clip(dis.d0 + dis.sd * z, 0.0, dis.d_max)


def field_yields(sim, shares, dis: Disease, seed=11):
    """Field yield each season, t/ha: share-weighted variety yields, minus pest and disease
    losses. Y_{i,t} = (1 - D_{i,t}) * sum_v c_v Y_v. Shape (n_fields, n_crops)."""
    y_mix = (shares[:, :, None] * sim["y_var"]).sum(1)
    D = disease_losses(shares, sim["n_crops"], dis, seed)
    return (1.0 - D) * y_mix


def annual(season_yields):
    """Add up the two seasons in each year: (n_fields, n_years)."""
    nf, nc = season_yields.shape
    ny = nc // 2
    return season_yields[:, : 2 * ny].reshape(nf, ny, 2).sum(-1)


# ==========================================================================
# A2: regional synchrony and shortfall statistics
# ==========================================================================

def partition(production):
    """Split regional variability into farm variability and synchrony (Wang and Loreau 2014).

    production : (n_fields, n_years) in tonnes
    Returns CV_L (typical farm's coefficient of variation, weighted by size),
    CV_R (coefficient of variation of regional total) and phi (synchrony, 0 to 1),
    which satisfy CV_R^2 = phi * CV_L^2."""
    m = production.mean(1)                 # mean production of each farm
    cov = np.cov(production)               # covariance between every pair of farms
    sd = np.sqrt(np.diag(cov))             # standard deviation of each farm
    cv_l = sd.sum() / m.sum()
    cv_r = np.sqrt(cov.sum()) / m.sum()
    phi = cov.sum() / sd.sum() ** 2
    return dict(CV_L=cv_l, CV_R=cv_r, phi=phi)


def shortfall_stats(yields_annual, R):
    """Shortfall measures for production below a need R (same units as yields_annual):
    p_short   share of farm-years in shortfall
    exp_short mean shortfall size, counting zero for years without one
    mean_run  mean number of consecutive years in shortfall"""
    short = yields_annual < R
    p = short.mean()
    S = np.maximum(R - yields_annual, 0).mean()
    # mean length of runs of consecutive shortfall years, farm by farm
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


# ==========================================================================
# A3: storage, prices, households
# ==========================================================================

@dataclass
class Market:
    """Regional grain market and household settings."""
    theta: float = 0.2        # market integration, 0 isolated (local harvest sets price), 1 fully integrated
    eps: float = 0.3          # price elasticity of demand (absolute); low means prices swing a lot
    p_bar: float = 250.0      # local price in a normal year, USD/t
    sigma_w: float = 0.15     # year-to-year variability of the outside (import parity) price, log scale
    delta: float = 0.10       # share of stored grain lost per year
    s_max_frac: float = 0.25  # storage capacity as a share of normal demand
    need: float = 1.8         # household grain need, t/yr
    cash: float = 150.0       # household cash available for buying grain, USD/yr
    c_bar: float = None       # normal regional demand, t/yr (None: mean regional harvest)


def market(production, area, mk: Market, seed=21):
    """Storage, prices and household food gaps, year by year.

    production : (n_fields, n_years) tonnes; each field is one household
    Returns the price series and, for each household, whether it is a net seller,
    its food gap (t) and the first-order welfare effect of price changes (USD)."""
    rng = np.random.default_rng(seed)
    nf, ny = production.shape
    # Normal regional demand is set to the landscape's mean harvest, so that p_bar is
    # the price in a normal year and prices respond to departures from normal.
    c_bar = production.sum(0).mean() if mk.c_bar is None else mk.c_bar
    s_max = mk.s_max_frac * c_bar
    p_w = mk.p_bar * np.exp(mk.sigma_w * rng.standard_normal(ny))     # outside price each year

    # Storage: grain above normal demand goes into store (up to capacity); in poor years
    # the store is drawn down. Q is the grain actually available for consumption each year.
    P = production.sum(0)
    S = 0.0
    Q = np.zeros(ny)
    for y in range(ny):
        A = P[y] + (1 - mk.delta) * S          # harvest plus what is left in store
        S = min(s_max, max(A - c_bar, 0.0))    # store any surplus, up to capacity
        Q[y] = A - S

    # Price: a weighted mix (in logs) of the local supply-demand price and the outside price
    lnp = (1 - mk.theta) * (np.log(mk.p_bar) - np.log(Q / c_bar) / mk.eps) + mk.theta * np.log(p_w)
    price = np.exp(lnp)

    # Households: net position, food gap and welfare effect of the price (Deaton 1989)
    net = production - mk.need                              # surplus (+) or deficit (-), t
    seller = net.mean(1) > 0                                 # net seller on average
    gap = np.maximum(mk.need - production - mk.cash / price[None, :], 0.0)   # need not grown or bought
    welfare = net * (price[None, :] - mk.p_bar)             # gain (+) or loss (-) from the price, USD
    return dict(price=price, Q=Q, P=P, seller=seller, gap=gap, welfare=welfare)
