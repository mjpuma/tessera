# tessera

Do crop variety mixtures stabilize food supply beyond the field?

`tessera` is a prototype model that follows the effect of crop mixtures from a single field up to regional food supply and local grain markets. A tessera is one tile in a mosaic: each farm is a tile, and the question is whether the whole picture holds together in a bad year.

The model was developed as a possible modeling component for a proposed study of sorghum and maize variety mixtures in East and West Africa. **All parameter values are illustrative and uncalibrated.** Results show how the model behaves, not findings about real farming systems.

## The model

The model has three linked parts. Field $i$ is a farm, $t$ is a year and $v$ is a variety.

### 1. Soil water and yields

Each field has a daily stochastic soil water balance in the tradition of Rodriguez-Iturbe et al. (1999) and Laio et al. (2001):

$$nZ_r \frac{ds}{d\tau} = I(\tau) - \sum_v c_v\,T_v(s) - E(s) - \Lambda(s)$$

- $s$ is relative soil moisture, from 0 (dry) to 1 (saturated), and $\tau$ is the day.
- $n$ is soil porosity and $Z_r$ is root zone depth (mm), so $nZ_r$ is the water the root zone can hold.
- $I$ is rain that infiltrates. Storms arrive as a Poisson process with exponentially distributed depths, and rain beyond saturation runs off.
- $c_v$ is the planting share of variety $v$ in the field (1 for a uniform stand). All varieties in a field draw on the same soil water.
- $T_v$ is transpiration by variety $v$. It is at its maximum when $s$ is above the variety's stress threshold $s^*_v$, falls linearly below it, and stops at the wilting point $s_{w,v}$.
- $E$ is soil evaporation and $\Lambda$ is leakage (drainage below the root zone).

The daily chance of a storm varies by season, by year with ENSO, and from place to place. Year-to-year anomalies are correlated between fields according to $e^{-d_{ij}/L}$, where $d_{ij}$ is the distance between fields and $L$ is the **rainfall correlation length** (km). A short $L$ means neighboring farms often have different rainfall.

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

## Running it

```
pip install -r requirements.txt
python run_experiments.py
```

This writes `results.csv` and the figures in `figures/`. A full run takes about two minutes on a laptop.

## Files

| File | Contents |
|---|---|
| `mixtures_model.py` | The model: rainfall generator, soil water balance, yields, disease, synchrony, markets |
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

- One root zone, so differences in rooting depth between varieties are not represented.
- Fixed planting dates. In practice farmers wait for the rains to start.
- Varieties in a mixture interact only through shared soil water. There is no light competition and no disease dilution within the field.
- One regional market with a simple storage rule. The competitive storage model would be the standard alternative.
- Seed networks are not yet implemented.

## References

Deaton, A. (1989). Rice prices and income distribution in Thailand: a non-parametric analysis. *Economic Journal*, 99(395), 1–37.

Doorenbos, J., and Kassam, A. H. (1979). *Yield Response to Water*. FAO Irrigation and Drainage Paper 33. Rome: FAO.

Laio, F., Porporato, A., Ridolfi, L., and Rodriguez-Iturbe, I. (2001). Plants in water-controlled ecosystems: active role in hydrologic processes and response to water stress. II. Probabilistic soil moisture dynamics. *Advances in Water Resources*, 24(7), 707–723.

Rodriguez-Iturbe, I., Porporato, A., Ridolfi, L., Isham, V., and Cox, D. R. (1999). Probabilistic modelling of water balance at a point: the role of climate, soil and vegetation. *Proceedings of the Royal Society A*, 455(1990), 3789–3805.

Wang, S., and Loreau, M. (2014). Ecosystem stability in space: α, β and γ variability. *Ecology Letters*, 17(8), 891–901.

## Contact

Michael J. Puma, Center for Climate Systems Research, Columbia University.
