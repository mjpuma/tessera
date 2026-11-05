# tessera

Do crop variety mixtures stabilize food supply beyond the field?

`tessera` is a prototype model that follows the effect of crop mixtures from a single field up to regional food supply and local grain markets. A tessera is one tile in a mosaic: each farm is a tile, and the question is whether the whole picture holds together in a bad year.

The model was developed as a possible modeling component for a proposed study of sorghum and maize variety mixtures in East and West Africa. **All parameter values are illustrative and uncalibrated.** Results show how the model behaves, not findings about real farming systems.

## The model

The model has three linked parts. Field $i$ is a farm, $t$ is a year and $v$ is a variety.

### 1. Soil water and yields

Each field has a daily stochastic soil water balance in the tradition of Rodriguez-Iturbe et al. (1999) and Laio et al. (2001), in `soilwater.py`:

$$nZ_r \frac{ds}{d\tau} = I(\tau) - \sum_v c_v\left[E_w\,r_v(s) + (E_{max,v} - E_w)\,f_v(s)\right] - (1 - C)\,E_b\,b(s) - \Lambda(s)$$

- $s$ is relative soil moisture, from 0 (dry) to 1 (saturated), and $\tau$ is the day.
- $n$ is soil porosity and $Z_r$ is root zone depth (mm), so $nZ_r$ is the water the root zone can hold.
- $I$ is rain that infiltrates. Each day's rain is the total of a Poisson number of storms with exponentially distributed depths, applied at the start of the day. Rain beyond saturation runs off.
- $c_v$ is the canopy cover of variety $v$: its planting share times its cover at that point in the season. $C = \sum_v c_v$. All varieties in a field draw on the same soil water.
- $f_v(s)$ is the stomatal control of transpiration. It is 1 above the variety's stress threshold $s^*_v$, falls linearly to 0 at its wilting point $s_{w,v}$, and $E_{max,v}$ is its maximum evapotranspiration.
- $r_v(s)$ ramps evaporation and cuticular losses $E_w$ from 0 at the hygroscopic point $s_h$ to 1 at $s_{w,v}$.
- $E_b\,b(s)$ is evaporation from bare soil, rising linearly from $s_h$ to field capacity $s_{fc}$.
- $\Lambda(s) = K_s\left(e^{\beta(s - s_{fc})} - 1\right) / \left(e^{\beta(1 - s_{fc})} - 1\right)$ for $s > s_{fc}$ is leakage below the root zone.

For a single variety at full cover ($c_v = 1$) this is exactly the loss function of Laio et al. (2001), Eq. 16. Losses are integrated with adaptive second-order sub-steps. The scheme is checked against the published analytical results in [`validation/`](validation/).

The storm rate varies by season, by year with ENSO, and from place to place. Year-to-year anomalies are correlated between fields according to $e^{-d_{ij}/L}$, where $d_{ij}$ is the distance between fields and $L$ is the **rainfall correlation length** (km). A short $L$ means neighboring farms often have different rainfall.

Each variety's yield falls with its transpiration deficit in each growth stage (vegetative, flowering, grain fill), weighted by yield response factors $K_y$ (Doorenbos and Kassam 1979). Varieties differ in maturity, so their flowering windows fall on different days. Pest and disease losses are more correlated between fields that grow the same varieties.

### 2. From fields to the region

Regional production is the sum over farms. Its variability is split into how variable individual farms are and how much they move together, using the identity of Wang and Loreau (2014):

$$CV_R^2 = \varphi \, CV_L^2$$

where, with $m_i$ the mean production of farm $i$ and $\sigma_{ij}$ the covariance of production between farms $i$ and $j$ (so $\sigma_i = \sqrt{\sigma_{ii}}$ is farm $i$'s standard deviation):

$$CV_L = \frac{\sum_i \sigma_i}{\sum_i m_i}, \qquad CV_R = \frac{\sqrt{\sum_i \sum_j \sigma_{ij}}}{\sum_i m_i}, \qquad \varphi = \frac{\sum_i \sum_j \sigma_{ij}}{\left(\sum_i \sigma_i\right)^2}$$

- $CV_L$ is **local variability**: the coefficient of variation of a typical farm's production from year to year, weighted by farm size.
- $CV_R$ is **regional variability**: the coefficient of variation of total regional production.
- $\varphi$ is **synchrony**, between 0 and 1. It is 1 when all farms rise and fall together and close to 0 when their good and bad years are independent.

Regional supply can therefore be steady for two reasons: farms are individually steady (low $CV_L$), or their bad years do not line up (low $\varphi$).

### 3. Markets, storage and households

Regional grain goes into a storage bucket. Each year, grain above normal demand is stored up to a capacity $S_{max}$, stored grain loses a share $\delta$ per year, and stocks are drawn down in poor years. The local price $p_t$ then depends on the grain available, $Q_t$:

$$\ln p_t = (1-\theta)\left[\ln \bar{p} - \frac{1}{\varepsilon}\ln\frac{Q_t}{\bar{C}}\right] + \theta \ln p^w_t$$

- $\bar{p}$ is the price in a normal year and $\bar{C}$ is normal regional demand, set to the mean regional harvest.
- $\varepsilon$ is the price elasticity of demand. A low value means prices swing strongly with the harvest.
- $p^w_t$ is the outside (import parity) price.
- $\theta$ is **market integration**, from 0 to 1. At 0 the local price is set entirely by the local harvest; at 1 it is set entirely by outside prices.

For each household, the **food gap** is the grain it needs but can neither grow nor buy:

$$G_{h,t} = \max\left(R_h - y_{h,t} - \frac{I_{h,t}}{p_t},\ 0\right)$$

where $R_h$ is the household's grain need, $y_{h,t}$ its harvest and $I_{h,t}$ the cash it can spend on grain. A **net buyer** is a household that on average harvests less grain than it needs; a **net seller** harvests more. A price rise helps net sellers and hurts net buyers (Deaton 1989).

A fourth part, seed exchange networks and the spread of new varieties, is described in the design but not yet implemented.

### Outputs

- **Shortfall risk:** the share of household-years in which a household's harvest is below its grain need.
- **Mean run length:** the average number of consecutive years in shortfall. The model has memory through soil moisture, ENSO and grain stocks, so bad years can cluster.
- **Regional deficit:** the share of years in which total regional production is below total household need.

## Example output

![Soil moisture and flowering windows](figures/fig1_mechanism.png)

![Scenario results](figures/fig2_scenarios.png)

These runs use 64 farms over 60 years with 3 rainfall ensembles, and five landscapes: all improved variety, all landrace mixture, a mosaic of single landraces, half improved and half mixture, and the improved variety added to the mixture.

## Validation

The soil water model is checked against published results in [`validation/`](validation/):

- **Laio et al. (2001).** For a single variety at full cover the model is their model, and it reproduces their closed-form drydown (Eq. 20), steady-state soil moisture pdfs (Eq. 29) and water balance (Eqs. 36-44).
- **Puma et al. (2007).** Repeating their rainfall-resolution experiment with the bucket shows that daily rainfall, and the daily scheme used here, give the same root-zone soil moisture distribution as 15-minute rainfall, while 2-week and monthly rainfall do not.

The bucket does not capture the effects of the vertical soil profile. These are discussed under limitations in the validation notes.

![Laio et al. (2001) pdfs reproduced](validation/laio_fig8_pdfs.png)

![Puma et al. (2007) rainfall resolution reproduced](validation/puma_fig5_pdfs.png)

## Running it

```
pip install -r requirements.txt
python run_experiments.py
```

This writes `results.csv` and the figures in `figures/`. A full run takes about two minutes on a laptop.

## Files

| File | Contents |
|---|---|
| `soilwater.py` | Point soil water balance (loss function of Laio et al. 2001) and its integrator |
| `mixtures_model.py` | The model: rainfall generator, yields, disease, synchrony, markets |
| `validation/` | Checks of the soil water model against Laio et al. (2001) and Puma et al. (2007) |
| `run_experiments.py` | Scenario experiments across rainfall correlation length, market integration and yield potential |
| `make_figures.py` | Figures from `results.csv` |
| `results.csv` | Results of the example runs |

## Where to change things

| What | Where |
|---|---|
| Variety traits (yield potential, stress threshold, wilting point, maturity, yield response) | `VARIETIES` in `mixtures_model.py` |
| Soil | `Soil` |
| Rainfall climatology, ENSO loading, correlation lengths | `Climate` |
| Pest and disease losses | `Disease` |
| Storage, demand elasticity, market integration, household need and cash | `Market` |
| Landscapes | `scenarios()` in `run_experiments.py` |

## Known simplifications

- One well-mixed root zone. Vertical redistribution and differences in rooting depth between varieties are not represented. A one-dimensional soil water flow model such as [SoilWaterFlow](https://github.com/mjpuma/SoilWaterFlow) could be used as a benchmark at a few sites (see [`validation/`](validation/)).
- Fixed planting dates. In practice farmers wait for the rains to start.
- Varieties in a mixture interact only through shared soil water. There is no light competition and no disease dilution within the field.
- One regional market with a simple storage rule. The competitive storage model would be the standard alternative.
- Seed networks are not yet implemented.

## References

Deaton, A. (1989). Rice prices and income distribution in Thailand: a non-parametric analysis. *Economic Journal*, 99(395), 1–37.

Doorenbos, J., and Kassam, A. H. (1979). *Yield Response to Water*. FAO Irrigation and Drainage Paper 33. Rome: FAO.

Laio, F., Porporato, A., Ridolfi, L., and Rodriguez-Iturbe, I. (2001). Plants in water-controlled ecosystems: active role in hydrologic processes and response to water stress. II. Probabilistic soil moisture dynamics. *Advances in Water Resources*, 24(7), 707–723.

Puma, M. J., Rodriguez-Iturbe, I., Celia, M. A., and Guswa, A. J. (2007). Implications of rainfall temporal resolution for soil-moisture and transpiration modeling. *Transport in Porous Media*, 68, 37–67.

Rodriguez-Iturbe, I., Porporato, A., Ridolfi, L., Isham, V., and Cox, D. R. (1999). Probabilistic modelling of water balance at a point: the role of climate, soil and vegetation. *Proceedings of the Royal Society A*, 455(1990), 3789–3805.

Wang, S., and Loreau, M. (2014). Ecosystem stability in space: α, β and γ variability. *Ecology Letters*, 17(8), 891–901.

## Contact

Michael J. Puma, Center for Climate Systems Research, Columbia University.
