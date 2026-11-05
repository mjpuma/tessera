"""
Point soil water balance shared by the tessera model and the validation scripts.

The loss function follows Laio et al. (2001), Eq. 16, with the evapotranspiration
of each variety weighted by its canopy cover:

    nZr ds/dt = I(t) - sum_v cov_v [ Ew r_v(s) + (Emax_v - Ew) f_v(s) ]
                     - (1 - sum_v cov_v) Eb b(s)
                     - L(s)

    r_v(s) = clip((s - s_h) / (s_w,v - s_h), 0, 1)       evaporation and cuticular losses
    f_v(s) = clip((s - s_w,v) / (s*_v - s_w,v), 0, 1)    stomatal control of transpiration
    b(s)   = clip((s - s_h) / (s_fc - s_h), 0, 1)        bare-soil evaporation
    L(s)   = Ks (exp(beta (s - s_fc)) - 1) / (exp(beta (1 - s_fc)) - 1)   for s > s_fc

For a single variety at full cover (cov = 1) this is exactly Laio et al. (2001),
Eq. 16, with E_w, E_max, s_h, s_w, s*, s_fc, K_s and beta as in that paper.

Rain arrives as a pulse at the start of each time step. Losses are integrated with
adaptive second-order (midpoint) sub-steps so that soil moisture changes by at most
about `max_ds` per sub-step, which keeps the scheme accurate during fast drainage after
large storms.
Units: mm and days (or any consistent pair).
"""
from dataclasses import dataclass
import numpy as np


# ==========================================================================
# Soil parameters
# ==========================================================================

@dataclass
class SoilParams:
    """Soil properties. Defaults: sandy loam after Laio et al. (2001), Table 1, in mm and days.
    Soil moisture s is relative: 0 is oven dry, 1 is saturated."""
    n: float = 0.43          # porosity
    z_r: float = 500.0       # root zone depth (mm)
    k_s: float = 800.0       # saturated hydraulic conductivity (mm/day)
    beta: float = 13.8       # leakage exponent (Laio et al. 2001, beta = 2b + 4)
    s_h: float = 0.14        # hygroscopic point
    s_fc: float = 0.56       # field capacity
    e_w: float = 0.1         # evaporation at the wilting point (mm/day)
    e_bare: float = 1.5      # bare-soil evaporation at field capacity (mm/day)

    @property
    def nzr(self):
        """Water the root zone holds when saturated (mm): porosity x root depth.
        Dividing a water depth by nzr converts it to a change in s."""
        return self.n * self.z_r


# ==========================================================================
# Loss functions
# ==========================================================================


def leakage(s, soil: SoilParams):
    """Drainage below the root zone (mm/day). Zero at or below field capacity, rising
    exponentially to K_s at saturation (Laio et al. 2001, Eq. 14)."""
    out = np.zeros_like(s)
    m = s > soil.s_fc
    out[m] = soil.k_s * (np.exp(soil.beta * (s[m] - soil.s_fc)) - 1.0) / (np.exp(soil.beta * (1.0 - soil.s_fc)) - 1.0)
    return out


def loss_rates(s, cov, e_max, s_star, s_w, soil: SoilParams, et_scale=1.0):
    """Instantaneous loss rates.

    s      : (nf,) relative soil moisture
    cov    : (nf, nv) canopy cover of each variety (planting share x cover x active)
    e_max, s_star, s_w : (nv,) variety parameters
    et_scale : scalar or (nf,) multiplier on atmospheric demand (e.g. diurnal cycle)
    Returns transpiration (nf, nv), evaporation (nf,), leakage (nf,).
    """
    sc = s[:, None]                     # column vector so it broadcasts across varieties
    # f: stomatal control, 1 above s*, falling linearly to 0 at the wilting point
    f = np.clip((sc - s_w) / (s_star - s_w), 0.0, 1.0)
    # r: evaporation and cuticular losses, 1 above the wilting point, 0 at the hygroscopic point
    r = np.clip((sc - soil.s_h) / (s_w - soil.s_h), 0.0, 1.0)
    et_scale = np.asarray(et_scale)
    k = et_scale[:, None] if et_scale.ndim == 1 else et_scale
    # Transpiration of each variety: its cover x (E_max - E_w) x stomatal control
    transp = k * cov * (e_max - soil.e_w) * f
    # Uncovered ground evaporates as bare soil, rising linearly from s_h to field capacity
    bare = np.clip(1.0 - cov.sum(1), 0.0, 1.0)
    b = np.clip((s - soil.s_h) / (soil.s_fc - soil.s_h), 0.0, 1.0)
    # Total evaporation: under the canopy (E_w term) plus bare soil
    evap = (k * cov * soil.e_w * r).sum(1) + (et_scale * bare * soil.e_bare * b)
    return transp, evap, leakage(s, soil)


# ==========================================================================
# Time stepping
# ==========================================================================


def infiltrate(s, rain, soil: SoilParams):
    """Add a rain pulse; water beyond saturation runs off. Returns (s, runoff)."""
    room = soil.nzr * (1.0 - s)             # space left in the root zone (mm)
    infil = np.minimum(rain, room)
    return s + infil / soil.nzr, rain - infil


def drain(s, dt, cov, e_max, s_star, s_w, soil: SoilParams, et_scale=1.0, max_ds=0.005):
    """Integrate losses over an interval dt (scalar or per field) with adaptive explicit sub-steps.

    Returns new s and the water removed over dt by each component:
    transpiration (nf, nv), evaporation (nf,), leakage (nf,).

    Each field takes its own sub-steps: short ones when losses are fast (just after a large
    storm, when leakage is strong) and few or one when the soil is dry. The loop runs
    until every field has covered the whole interval.
    """
    nzr = soil.nzr
    s = s.copy()
    # time still to integrate for each field, and running totals of each loss
    remaining = np.broadcast_to(np.asarray(dt, dtype=float), s.shape).copy()
    t_acc = np.zeros(cov.shape)
    e_acc = np.zeros(s.shape)
    l_acc = np.zeros(s.shape)
    active = remaining > 0
    while active.any():
        # second-order (midpoint) step: choose h from the rate at s, then use the rate at s(t + h/2)
        tr, ev, lk = loss_rates(s, cov, e_max, s_star, s_w, soil, et_scale)
        rate = tr.sum(1) + ev + lk                          # total loss rate, mm/day
        # step length: small enough that s changes by at most max_ds, and no longer than what is left
        h = np.where(rate > 0, np.minimum(remaining, max_ds * nzr / np.maximum(rate, 1e-12)), remaining)
        h = np.where(active, h, 0.0)                        # fields already finished take no step
        s_mid = np.maximum(s - 0.5 * rate * h / nzr, soil.s_h)
        tr, ev, lk = loss_rates(s_mid, cov, e_max, s_star, s_w, soil, et_scale)
        rate = tr.sum(1) + ev + lk
        # never remove more water than is present above the hygroscopic point
        avail = np.maximum(s - soil.s_h, 0.0) * nzr
        scale = np.where(rate * h > avail, avail / np.maximum(rate * h, 1e-12), 1.0)
        # record the water removed by each component, then update soil moisture
        t_acc += tr * (h * scale)[:, None]
        e_acc += ev * h * scale
        l_acc += lk * h * scale
        s = s - rate * h * scale / nzr
        remaining = remaining - h
        active = remaining > 1e-12
    return s, t_acc, e_acc, l_acc
