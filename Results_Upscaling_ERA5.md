# ERA5-Land Monthly Upscaling of Upland Methane Exchange

## Methods

### Analysis goal and study domain

We developed a spatial upscaling workflow to estimate non-wetland upland methane (CH4) exchange from 2000 through 2025 using ERA5-Land climate predictors, dynamic land-cover information, monthly inundation screening, and flux behavior learned from NEON tower-gradient observations. The analysis was designed around a two-stage ecological premise: upland ecosystems should not be assigned a fixed source or sink rate directly from ecosystem type alone, because a given upland location can behave as either a weak CH4 sink or weak CH4 source depending on climate, soil moisture, and site context. Therefore, the workflow first estimated the probability that a grid-cell month was a weak source, then estimated flux magnitude conditionally for weak-sink and weak-source states.

Wetland and inundated ecosystems were excluded because the objective was to isolate non-wetland upland exchange rather than reproduce the global wetland CH4 source. This exclusion was necessary because inundated and wetland ecosystems have fundamentally different transport pathways, substrate availability, redox structure, and source magnitudes than aerated upland soils. Including them would have allowed wetland-like emission rates to dominate the upscaled estimate and would have made the result inappropriate for comparison with terrestrial upland soil uptake in the Global Methane Budget.

**Table 1. Summary of major modeling stages and why each was needed.**

| Stage | Response or output | Main predictors or inputs | Why this step was necessary |
|---|---|---|---|
| Upland filtering | Non-wetland upland training and prediction domain | NEON ecosystem labels, MODIS upland classes, WAD2M inundation fraction | Prevents wetland or inundated methane-emission processes from being mixed with aerated upland soil exchange. |
| Monthly aggregation | Site-month CH4 budget and source state | ERA5-gapfilled 30-minute NEON CH4 fluxes | Reduces half-hourly noise and aligns the response with monthly ERA5-Land spatial predictors. |
| Class-probability model | `P(weak source)` for each site-month or grid-cell month | EcoType, ERA5-Land air temperature, ERA5-Land soil moisture, MAP, MAT | Separates the question of whether a location is likely to be a source from the question of how large the flux should be. |
| Probability calibration | Calibrated monthly source probability | Raw model probability | Provides a check that probability values used in expected-flux calculations track observed source frequency. |
| Conditional magnitude models | Weak-sink and weak-source flux magnitude | ERA5-Land air temperature, ERA5-Land soil moisture, MAP, MAT, calibrated source probability, EcoType, site random effect | Allows sink and source flux rates to respond to conditions instead of using fixed ecosystem lookup rates. |
| Continuous expected flux | Monthly expected CH4 flux | `P(source)`, modeled source flux, modeled sink flux | Avoids abrupt source/sink threshold behavior and lets near-threshold cells contribute intermediate fluxes. |
| Chamber/process constraint | Soft-constrained expected flux | Chamber/process sink rates and upland chamber source bound | Keeps non-wetland upland rates within plausible empirical/process-informed ranges without fully discarding model variation. |

### NEON training data and monthly flux aggregation

The model was trained from the ERA5-gapfilled NEON half-hourly CH4 flux product (`NEON_ERA5_gapfilled_30min.csv.gz`). Fluxes were first restricted to upland NEON sites by joining with the site behavior table and excluding ecosystem labels containing wetland, inundated, flooded, marsh, swamp, bog, fen, lake, or rice terms. This filtering ensured that both the class model and magnitude model were learned from the same non-wetland upland domain used for spatial prediction.

Half-hourly ERA5-gapfilled CH4 fluxes were aggregated to the site-month scale. For each `SITE_ID`, ecosystem type, year, and month, the workflow calculated monthly CH4 budget as the sum of half-hourly gapfilled fluxes:

```r
monthly_budget_mgC_m2 = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE)
```

This monthly scale was used because the scientific question concerns seasonal-to-interannual shifts between weak uptake and weak emission, and because monthly aggregation reduces the influence of individual half-hourly noise, short storage artifacts, and uneven observation density. A month was classified as a weak-source month when its total monthly CH4 exchange was positive:

```r
weak_source = as.integer(monthly_budget_mgC_m2 > 0)
```

The same monthly aggregation was also converted to g C m-2 month-1 for flux magnitude modeling:

```r
monthly_flux_gC_m2_month = monthly_budget_mgC_m2 / 1000
```

### Monthly weak-source probability model

The first stage estimated the probability that an upland site-month was a weak source. The monthly class model used a binomial generalized linear model with a quasibinomial error structure:

```r
weak_source ~ EcoType +
  scale(mean_ERA5_Tair_C) +
  scale(mean_ERA5_VSWC) +
  scale(MAP) +
  scale(MAT)
```

`mean_ERA5_Tair_C` and `mean_ERA5_VSWC` were monthly ERA5-derived conditions. `MAP` and `MAT` represented site-level mean annual precipitation and mean annual temperature, respectively. Ecosystem type was included because croplands, forests, grasslands, and shrublands differ in soil structure, plant inputs, management, and methane oxidation/emission potential. Season was intentionally excluded from the class-probability model because fixed northern-hemisphere seasonal labels do not transfer cleanly to a global prediction domain. Monthly temporal structure is therefore represented through actual monthly temperature and soil moisture rather than through a categorical season term.

The class model was fit with inverse class-frequency weights so that weak-source and weak-sink months contributed equally to model fitting. This balancing step was necessary because the NEON monthly training data were source-heavy. Without class balancing, the model could learn the sampling imbalance in the training set rather than the environmental conditions associated with source behavior. Class weights were calculated as:

```r
class_balance_weight = n_total / (number_of_classes * n_class)
```

where `n_class` was the number of training months in each sink/source class.

The predicted source probability was then calibrated with a second binomial model using the logit of the raw source probability as the only predictor. This calibration step was included so that the probability used in the continuous expected-flux calculation represented an observed monthly source frequency rather than only the linear predictor scale of the balanced classifier. In the current no-season run, calibration did not materially change the probability values, but it provides a formal checkpoint for future probability-model revisions.

**Table 2. Class-probability model design.**

| Feature | Implementation |
|---|---|
| Training unit | NEON upland site-month |
| Response | `weak_source = 1` when monthly ERA5-gapfilled CH4 budget was positive |
| Model family | Quasibinomial GLM |
| Predictors | EcoType, monthly ERA5-Land air temperature, monthly ERA5-Land soil moisture, MAP, MAT |
| Class balancing | Inverse class-frequency weights so weak-source and weak-sink months had equal total weight |
| Season term | Excluded |
| Default source threshold | `P(source) >= 0.80` for hard-threshold summaries |
| Spatial prediction unit | ERA5-Land grid-cell month from 2000-2025 |

### Conditional flux magnitude models

The second stage estimated flux magnitude separately for weak-sink and weak-source months. This separation was necessary because the mechanisms and magnitude distributions for uptake and emission are not symmetric. Upland CH4 uptake is constrained by diffusion, methanotrophic activity, and soil aeration, whereas weak upland emissions can occur under wetter microsites, pulses, or storage/transport conditions. A single linear flux model would tend to average these states and could produce unstable near-zero predictions.

For each state, the model used log-absolute monthly flux magnitude:

```r
log_abs_flux = log(max(abs(monthly_flux_gC_m2_month), 1e-6))
```

The log transformation was used because monthly CH4 flux magnitudes are right-skewed and include many near-zero values. Modeling log magnitude stabilizes variance and reduces the leverage of large source months. After prediction, the sign was restored by assigning negative values to weak-sink predictions and positive values to weak-source predictions.

Separate hierarchical linear mixed-effects models were fit for weak-sink and weak-source months:

```r
log_abs_flux ~ z_Tair + z_VSWC + z_MAP + z_MAT + source_probability +
  (1 + z_Tair + z_VSWC | EcoType) + (1 | SITE_ID)
```

where `z_Tair`, `z_VSWC`, `z_MAP`, and `z_MAT` were standardized predictors using the means and standard deviations of the NEON monthly training data. The model allowed ecosystem-specific intercepts and, where supported by the data, ecosystem-specific temperature and soil moisture slopes. This partial-pooling structure was necessary because fitting fully separate models by ecosystem type would be unstable for underrepresented classes, while complete pooling would assume that croplands, forests, grasslands, and shrublands respond identically to climate. A site random effect was included during training to absorb persistent site-level differences that should not be misattributed to climate predictors. During spatial prediction, site-level effects were not used; predictions were made with new spatial locations treated as unseen sites.

If the random-slope model was singular or failed, the workflow fell back to a simpler mixed model with ecosystem-specific intercepts and site random effects. If that also failed, the workflow used a fixed-effect ecosystem model. This fallback logic was necessary because the sink training data were much smaller than the source training data, and mixed-effect random-slope models can be overparameterized when some ecosystem-by-state combinations have few observations.

**Table 3. Conditional magnitude model design.**

| Feature | Weak-sink model | Weak-source model |
|---|---|---|
| Training subset | Site-months with negative monthly flux | Site-months with positive monthly flux |
| Response scale | `log(abs(monthly_flux_gC_m2_month))` | `log(abs(monthly_flux_gC_m2_month))` |
| Sign restoration | Predicted magnitude multiplied by -1 | Predicted magnitude retained as positive |
| Fixed predictors | Standardized Tair, VSWC, MAP, MAT, calibrated source probability | Standardized Tair, VSWC, MAP, MAT, calibrated source probability |
| Random effects | EcoType partial pooling and site intercept | EcoType partial pooling and site intercept |
| Spatial prediction | Population-level prediction for unseen spatial cells | Population-level prediction for unseen spatial cells |

### Continuous expected-flux calculation

The primary continuous upscaling estimate avoided hard source/sink classification. For every grid-cell month, expected flux was calculated as:

```r
expected_flux =
  P(source) * source_flux +
  (1 - P(source)) * sink_flux
```

This step was necessary because many upland pixels are near the source/sink decision boundary. A hard threshold forces each pixel-month into a full source or full sink rate, producing abrupt spatial and temporal discontinuities. The continuous expected-flux formulation allows near-threshold locations to contribute near-zero exchange and treats source behavior probabilistically rather than categorically.

### Chamber/process constraints on magnitude

Two continuous magnitude scenarios were evaluated. The first, condition-only magnitude, used the hierarchical magnitude model predictions directly. The second, chamber-constrained magnitude, applied soft shrinkage toward independent chamber/process constraints. The constrained scenario shrank sink rates toward chamber/process sink rates and shrank high source rates toward the maximum positive non-wetland upland chamber source bound available in the local chamber reference table. FLUXNET-CH4 source classes were not used for source bounds because the relevant FLUXNET source ecosystems are inundated or partially inundated and are not comparable to the non-wetland upland domain analyzed here.

The constraint was applied as soft shrinkage rather than hard clipping. For sink fluxes, predicted rates were blended with chamber/process sink rates using a constraint weight of 0.65. For source fluxes, only predictions above the upland chamber source bound were shrunk toward that bound. Soft shrinkage was used because hard caps can create artificial plateaus and discard environmental signal, whereas shrinkage preserves conditional variation while preventing implausibly large non-wetland upland source rates.

### Spatial predictors and land-surface filtering

Spatial prediction used monthly ERA5-Land grids for 2000-2025. ERA5-Land 2 m temperature was converted from Kelvin to degrees Celsius. ERA5-Land volumetric soil water layers 1 and 2 were averaged to represent monthly near-surface soil moisture. ERA5-Land total precipitation was converted to monthly millimeters and used to derive long-term MAP. Long-term MAT was calculated from the monthly temperature fields.

Dynamic land cover was incorporated from processed MODIS MCD12C1 annual land-cover rasters. MODIS classes were mapped to four broad upland ecosystem types: cropland, forest, grassland, and shrubland. When a prediction year was outside the available MODIS range, the nearest available processed MODIS year was used. This dynamic land-cover step was necessary because broad ecosystem area changes through time and because cropland extent cannot be represented reliably from the ecoregion fallback alone.

Monthly inundation was screened with WAD2M inundation fraction. Cells with inundation fraction greater than 0.05 were excluded for that month, and remaining cells were area-weighted by `1 - inundation_fraction`. This step was necessary to maintain the upland/non-inundated domain even in locations where a nominally upland land-cover class included seasonally inundated fractions.

### Area weighting and conversion to Tg CH4 yr-1

For each ERA5-Land grid cell, cell area was calculated in million hectares. Monthly carbon flux in g C m-2 month-1 was converted to Tg CH4 by multiplying by area and converting carbon mass to methane mass:

```r
Tg CH4 = flux_gC_m2_month * area_mha * 0.0133333333333333
```

This conversion includes the relationship between m2 and Mha, grams and teragrams, and molecular conversion from C to CH4. Monthly cell contributions were summed by year to obtain annual net exchange in Tg CH4 yr-1. Positive values indicate net emission to the atmosphere, and negative values indicate net uptake.

### Threshold-based comparison workflow

The continuous expected-flux estimate was evaluated alongside the earlier hard-threshold workflow. In the threshold workflow, a grid-cell month was assigned to weak source if its source probability exceeded a specified threshold and weak sink otherwise. The default threshold was 0.80. Sensitivity was evaluated from 0.50 to 0.95 in increments of 0.05. This threshold analysis was retained because it directly tests the decision rule implied by the weak-source probability model and identifies the threshold at which the upscaled budget changes sign.

### Class-change diagnostics

To determine whether spatial locations were stable or changed source/sink state through time, the workflow summarized class changes in two ways. Monthly switching counted cells that changed between weak sink and weak source at least once across all 312 months from 2000-2025. Annual switching first classified each cell-year as weak source if at least six months were assigned weak source, then counted cells whose annual dominant class changed at least once across years. These diagnostics were necessary because a monthly probability model can produce seasonal switching even when the annual dominant class is stable; reporting both scales separates within-year seasonality from interannual shifts.

## Results

### Source probability and spatial class dynamics

The no-season monthly source-probability model predicted a more sink-dominated upland domain than the earlier seasonal class model. At the default weak-source probability threshold of 0.80, mean monthly weak-source area from 2000-2025 was 1,081 Mha, ranging from 1,041 to 1,164 Mha among years. Mean monthly weak-sink area was 6,780 Mha, ranging from 6,608 to 7,095 Mha. These results indicate that removing season from the class-probability model substantially reduced the area assigned to weak-source behavior under the hard-threshold decision rule.

The class-probability model had moderate ranking strength but was intentionally conservative at the default 0.80 source threshold. Across 1,764 NEON upland site-months, 1,482 months were positive and 282 were negative. The model AUC was 0.683, indicating that the model ranked a randomly selected source month above a randomly selected sink month about 68% of the time. The Brier score was 0.230, Tjur's R2 was 0.081, and deviance explained was 5.93%. These values indicate that the no-season climate-and-ecosystem model contains useful but incomplete information about monthly source state. At the 0.80 hard threshold, training-set source sensitivity was only 2.16%, while sink specificity was 98.6%; this is expected because the 0.80 threshold was chosen as a strong-evidence source rule rather than as an accuracy-optimizing threshold.

**Table 4. No-season monthly source-probability model strength.**

| Diagnostic | Value | Interpretation |
|---|---:|---|
| Site-months | 1,764 | NEON upland monthly training observations |
| Sink months | 282 | Monthly budget <= 0 |
| Source months | 1,482 | Monthly budget > 0 |
| Source fraction | 0.840 | Training data are source-heavy |
| AUC | 0.683 | Moderate ability to rank source months above sink months |
| Brier score | 0.230 | Mean squared probability error |
| Tjur's R2 | 0.081 | Mean probability separation between source and sink months |
| Deviance explained | 5.93% | Limited but nonzero explanatory strength |
| Accuracy at 0.80 threshold | 17.6% | Low because threshold strongly defaults to sink |
| Source sensitivity at 0.80 | 2.16% | Few training source months exceed the strong-evidence threshold |
| Sink specificity at 0.80 | 98.6% | Most sink months are correctly kept as sink |
| Source precision at 0.80 | 88.9% | Months classified as source are usually observed source months |

Class switching also became much less common after removing season. Across the 53,708 modeled upland grid cells, 3,916 cells changed between weak sink and weak source at least once from 2000-2025. This represented 7.29% of cells and 1,117 Mha, or 9.38% of the modeled upland area. Annual dominant-class switching was smaller still: 1,892 cells changed annual dominant class at least once, representing 3.52% of cells and 545 Mha, or 4.57% of modeled upland area. This reduction indicates that the previous season term was driving much of the monthly class switching, whereas the no-season model assigns more stable source/sink behavior based on climate conditions and ecosystem type.

Grasslands accounted for the largest area with annual class switching, with 223 Mha changing annual dominant class across 2000-2025. Forests contributed 198 Mha, croplands 118 Mha, and shrublands 6.65 Mha. In 2025, the mapped annual dominant class included 8,223 Mha of weak-sink area and 1,116 Mha of weak-source area. Croplands contained the largest absolute weak-source area in 2025 (467 Mha), followed by forests (375 Mha), grasslands (272 Mha), and shrublands (2.31 Mha).

**Table 5. Spatial class dynamics under the no-season probability model.**

| Scale | Cells changing class | Percent of cells | Area changing class (Mha) | Percent of modeled area |
|---|---:|---:|---:|---:|
| Monthly class over 2000-2025 | 3,916 | 7.29% | 1,117 | 9.38% |
| Annual dominant class over 2000-2025 | 1,892 | 3.52% | 545 | 4.57% |

### Threshold sensitivity of hard-classified estimates

The hard-threshold workflow showed that the global estimate was sensitive to the probability threshold used to assign weak-source behavior. After removing season from the class-probability model, both fixed-rate hard-threshold scenarios were net sinks at the default 0.80 threshold. Under the NEON ERA5 fixed-rate scenario, mean annual net exchange was -5.40 Tg CH4 yr-1. Under the chamber/process sink plus NEON source-rate scenario, mean annual net exchange was -15.85 Tg CH4 yr-1.

Increasing the source-probability threshold made the estimate progressively more sink-like by requiring stronger evidence before a cell-month could be assigned to the weak-source class. The chamber/process sink plus NEON source-rate scenario crossed from net source to net sink at a threshold of 0.75, where the mean annual estimate was -8.73 Tg CH4 yr-1. The NEON ERA5 fixed-rate scenario crossed from net source to net sink at the default threshold of 0.80, where the mean annual estimate was -5.40 Tg CH4 yr-1. This threshold behavior indicates that the no-season model is substantially more conservative about source classification than the previous seasonal probability model.

**Table 6. Annual hard-threshold estimates at the default 0.80 source threshold.**

| Rate scenario | Minimum (Tg CH4 yr-1) | Mean (Tg CH4 yr-1) | Maximum (Tg CH4 yr-1) | SD (Tg CH4 yr-1) |
|---|---:|---:|---:|---:|
| NEON ERA5 rates | -5.89 | -5.40 | -4.65 | 0.293 |
| Chamber/process sink + NEON source | -16.86 | -15.85 | -14.74 | 0.539 |

### Magnitude model performance

The hierarchical log-magnitude models were refit using the no-season source probabilities as predictors. The weak-sink model was trained on 282 sink months and had a mean observed flux of -0.0173 g C m-2 month-1. Its fitted mean was -0.0115 g C m-2 month-1, with RMSE of 0.0139 g C m-2 month-1, MAE of 0.0100 g C m-2 month-1, and observed-fitted correlation of 0.452.

The weak-source model was trained on 1,482 source months and had a mean observed flux of +0.0439 g C m-2 month-1. Its fitted mean was +0.0405 g C m-2 month-1, with RMSE of 0.0178 g C m-2 month-1, MAE of 0.0134 g C m-2 month-1, and observed-fitted correlation of 0.848. Source-magnitude skill remained high, but sink-magnitude skill declined relative to the seasonal-probability run because the source-probability predictor now contains less monthly timing information.

**Table 7. Conditional magnitude model strength.**

| Magnitude model | Training months | Observed mean (g C m-2 month-1) | Fitted mean (g C m-2 month-1) | RMSE | MAE | Observed-fitted correlation |
|---|---:|---:|---:|---:|---:|---:|
| Weak sink | 282 | -0.0173 | -0.0115 | 0.0139 | 0.0100 | 0.452 |
| Weak source | 1,482 | +0.0439 | +0.0405 | 0.0178 | 0.0134 | 0.848 |

**Table 8. Ecosystem-specific magnitude model diagnostics.**

| Magnitude model | EcoType | Training months | Observed mean | Fitted mean | RMSE | MAE | Correlation |
|---|---|---:|---:|---:|---:|---:|---:|
| Weak sink | Cropland | 4 | -0.0073 | -0.0059 | 0.0039 | 0.0037 | 0.974 |
| Weak sink | Forest | 169 | -0.0173 | -0.0118 | 0.0138 | 0.0099 | 0.461 |
| Weak sink | Grassland | 69 | -0.0198 | -0.0114 | 0.0172 | 0.0133 | 0.462 |
| Weak sink | Shrubland | 40 | -0.0137 | -0.0112 | 0.0071 | 0.0054 | 0.364 |
| Weak source | Cropland | 80 | +0.0373 | +0.0351 | 0.0147 | 0.0126 | 0.321 |
| Weak source | Forest | 851 | +0.0434 | +0.0399 | 0.0182 | 0.0136 | 0.844 |
| Weak source | Grassland | 363 | +0.0573 | +0.0541 | 0.0206 | 0.0158 | 0.828 |
| Weak source | Shrubland | 188 | +0.0227 | +0.0194 | 0.0101 | 0.0083 | 0.743 |

The probability calibration diagnostic showed that observed source frequency generally increased across predicted source-probability bins, although the no-season model compressed the probability range. The lowest probability bin had mean predicted probability of 0.273 and observed source frequency of 0.667, whereas the highest probability bin had mean predicted probability of 0.779 and observed source frequency of 0.938. This pattern indicates that the model still ranks source-prone months, but excluding season reduces separation among months. It also shows that even low-probability bins contain many observed source months in the NEON training data, which remains important for interpreting positive continuous expected-flux estimates.

**Table 9. Probability calibration by decile.**

| Probability bin | n | Mean predicted P(source) | Observed source fraction |
|---:|---:|---:|---:|
| 1 | 177 | 0.273 | 0.667 |
| 2 | 177 | 0.382 | 0.774 |
| 3 | 177 | 0.438 | 0.842 |
| 4 | 177 | 0.488 | 0.627 |
| 5 | 176 | 0.526 | 0.898 |
| 6 | 176 | 0.550 | 0.858 |
| 7 | 176 | 0.577 | 0.915 |
| 8 | 176 | 0.609 | 0.966 |
| 9 | 176 | 0.657 | 0.920 |
| 10 | 176 | 0.779 | 0.938 |

### Continuous expected-flux estimates

The continuous expected-flux approach reduced abrupt threshold behavior by weighting weak-source and weak-sink magnitude predictions by monthly source probability. Unlike the hard-threshold workflow, the continuous expected-flux workflow remained positive after removing season because source probabilities still contributed continuously even when they were below the 0.80 hard-classification threshold. Across 2000-2025, the condition-only continuous estimate averaged +30.3 Tg CH4 yr-1, with annual values ranging from +29.8 to +31.0 Tg CH4 yr-1.

The chamber-constrained continuous magnitude scenario produced a much smaller positive estimate. Across 2000-2025, the constrained continuous estimate averaged +5.97 Tg CH4 yr-1, with annual values ranging from +5.52 to +6.34 Tg CH4 yr-1. The constrained estimate was therefore much closer to zero than the unconstrained estimate, but it still did not reproduce the Global Methane Budget terrestrial soil sink reference of approximately -35 Tg CH4 yr-1.

Ecosystem-level continuous components showed why the constrained expected-flux estimate remained positive. In the chamber-constrained scenario, mean monthly contributions were positive for croplands (+0.229 Tg CH4 month-1), forests (+0.158 Tg CH4 month-1), and grasslands (+0.118 Tg CH4 month-1), while shrublands were slightly negative (-0.0076 Tg CH4 month-1). Mean source probabilities were 0.752 for croplands, 0.651 for forests, 0.580 for grasslands, and 0.580 for shrublands. The chamber constraints reduced source magnitudes substantially, but the continuous expected-flux calculation still allows sub-threshold source probability to contribute positive flux.

**Table 10. Annual continuous expected-flux estimates.**

| Magnitude scenario | Minimum (Tg CH4 yr-1) | Mean (Tg CH4 yr-1) | Maximum (Tg CH4 yr-1) | SD (Tg CH4 yr-1) |
|---|---:|---:|---:|---:|
| Condition-only magnitude | +29.82 | +30.27 | +30.96 | 0.350 |
| Chamber-constrained magnitude | +5.52 | +5.97 | +6.34 | 0.218 |

**Table 11. Mean monthly continuous components by ecosystem type.**

| Scenario | EcoType | Monthly contribution (Tg CH4 month-1) | Mean P(source) | Mean sink flux | Mean source flux | Mean constrained sink flux | Mean constrained source flux |
|---|---|---:|---:|---:|---:|---:|---:|
| Chamber-constrained | Cropland | +0.229 | 0.752 | -0.0059 | +0.0530 | -0.0043 | +0.0214 |
| Chamber-constrained | Forest | +0.158 | 0.651 | -0.0082 | +0.0505 | -0.0167 | +0.0205 |
| Chamber-constrained | Grassland | +0.118 | 0.580 | -0.0088 | +0.0391 | -0.0190 | +0.0165 |
| Chamber-constrained | Shrubland | -0.0076 | 0.580 | -0.0094 | +0.0250 | -0.0173 | +0.0116 |
| Condition-only | Cropland | +0.588 | 0.752 | -0.0059 | +0.0530 | -0.0043 | +0.0214 |
| Condition-only | Forest | +0.612 | 0.651 | -0.0082 | +0.0505 | -0.0167 | +0.0205 |
| Condition-only | Grassland | +1.109 | 0.580 | -0.0088 | +0.0391 | -0.0190 | +0.0165 |
| Condition-only | Shrubland | +0.214 | 0.580 | -0.0094 | +0.0250 | -0.0173 | +0.0116 |

### Comparison with the terrestrial soil sink reference

All spatial scenarios remained more positive than the Global Methane Budget terrestrial soil sink reference, but removing season changed the sign of the hard-threshold estimates. The fixed-rate hard-threshold estimates averaged -5.40 Tg CH4 yr-1 for NEON ERA5 rates and -15.85 Tg CH4 yr-1 for the chamber/process sink plus NEON source scenario. The continuous condition-only magnitude estimate averaged +30.3 Tg CH4 yr-1, while the continuous chamber-constrained estimate averaged +5.97 Tg CH4 yr-1.

These results suggest that the current ERA5-Land monthly framework improves the treatment of conditional source/sink behavior and flux magnitude, but it does not yet reconcile NEON-derived upland exchange with the global terrestrial soil sink. The no-season class-probability model produces a more conservative hard-threshold classification and a net sink under the default threshold. However, the continuous expected-flux model remains positive because sub-threshold source probabilities still contribute positive source flux across very large upland areas. The remaining uncertainty is therefore centered on how source probability should be interpreted for expected-flux estimation, not only on the binary class threshold.

## Figure Captions

**Figure 1. Monthly spatial upscaling components by ecosystem type and fixed-rate scenario.** Monthly CH4 exchange components from the hard-threshold upscaling workflow, averaged across 2000-2025. Bars show monthly Tg CH4 contributions by ecosystem type and assigned class using the default weak-source probability threshold of 0.80. Positive values indicate net CH4 emission and negative values indicate net uptake. This figure illustrates how the binary source/sink decision and fixed rate table distribute monthly exchange among croplands, forests, grasslands, and shrublands. File: `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/FIGURES/monthly_spatial_components.png`.

**Figure 2. Hard-threshold ERA5-Land upland CH4 upscaling diagnostics.** Multipanel summary of the hard-threshold source/sink workflow using the no-season class-probability model. (A) Mean annual net CH4 exchange from the fixed-rate hard-threshold workflow for the NEON ERA5 rate scenario and the chamber/process sink plus NEON source-rate scenario; the horizontal reference line and band indicate the Global Methane Budget terrestrial soil sink comparison value. (B) Annual hard-threshold net CH4 exchange from 2000-2025; after removing season from the class-probability model, both fixed-rate hard-threshold scenarios are net sinks at the default 0.80 source-probability threshold. (C) Inset within panel D showing mean monthly upland area assigned to weak-source and weak-sink classes in each year; the no-season class-probability model assigns most upland area to weak sink, with mean weak-source area of 1,081 Mha and mean weak-sink area of 6,780 Mha. (D) Enlarged map of dominant annual source/sink class for upland grid cells in 2025, where a cell is classified as weak source when at least six months exceed the 0.80 source-probability threshold. (E) Source-threshold sensitivity from 0.50 to 0.95; the chamber/process sink plus NEON source-rate scenario crosses from source to sink at 0.75, whereas the NEON ERA5 fixed-rate scenario crosses at 0.80. File: `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/FIGURES/ERA5_hard_threshold_multipanel.png`.

**Figure 3. Continuous expected-flux model diagnostics.** Multipanel summary of the continuous expected-flux workflow and model checks. (A) Annual net CH4 exchange from the continuous expected-flux approach, where expected flux is calculated as `P(source) * source_flux + (1 - P(source)) * sink_flux`; the condition-only magnitude scenario averages +30.3 Tg CH4 yr-1, whereas the chamber-constrained magnitude scenario averages +5.97 Tg CH4 yr-1. (B) Annual budget comparison for the continuous scenarios, with bars showing 2000-2025 means, whiskers showing the annual minimum-to-maximum range, and the Global Methane Budget terrestrial soil sink reference shown as the blue dashed line/band. (C) Observed versus fitted monthly CH4 flux for the separate weak-sink and weak-source log-magnitude models; the weak-sink model had RMSE of 0.0139 g C m-2 month-1 and observed-fitted correlation of 0.452, while the weak-source model had RMSE of 0.0178 g C m-2 month-1 and observed-fitted correlation of 0.848. (D) Binned calibration of the monthly source-probability model against observed source frequency in the NEON monthly training data; observed source frequency generally increases with predicted probability, but the no-season model compresses predictions into a narrower probability range and even the lowest predicted-probability bin contains many source months. File: `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/FIGURES/ERA5_continuous_expected_flux_multipanel.png`.

## Output Files

Key tables generated by the workflow include:

- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/annual_expected_flux_2000_2025.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/monthly_expected_flux_components.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/class_probability_model_skill.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/magnitude_model_skill.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/probability_calibration_skill.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/cell_class_change_totals_2000_2025.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/annual_cell_class_change_summary_2000_2025.csv`
- `/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/OUTPUT/spatial_threshold_sensitivity.csv`
