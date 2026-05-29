# ERA5 Methods And Results

## Methods

We used `NEON.ERA5.HalfHourlyGapfill.R` to estimate continuous annual CH4 budgets from ERA5-forced half-hourly gap filling. Hourly ERA5 point covariates were obtained for each NEON tower coordinate through the Open-Meteo Archive API and cached locally. The workflow extracted 2 m air temperature and 0-7 cm soil volumetric water content, then linearly interpolated hourly values to a 30-minute grid for each site-year represented in the NEON flux observations.

Observed total CH4 fluxes with matching ERA5 covariates were used to train a generalized additive model. The response was total CH4 flux in mg C m-2 per 30 min. Predictors included cyclic smooths for hour of day and day of year, smooth terms for ERA5 air temperature and ERA5 soil moisture, a tensor interaction between ERA5 temperature and moisture, season, ecosystem type, and a site random effect. The fitted model was evaluated using observed-versus-fitted fluxes, residuals, and population-level partial effects for the ERA5 covariates.

Annual ERA5 budgets retained observed half-hour fluxes where available and filled missing half-hours with ERA5-driven model predictions. Half-hour fluxes were summed within each site-year and converted to g C m-2 yr-1. Site-level ERA5 annual behavior was classified from the fraction of site-years with positive annual budgets: consistent source sites had at least 75% positive years, consistent sink sites had no more than 25% positive years, and all remaining sites were classified as fluctuating. ERA5 annual budgets were compared with model-standardized 30-minute annual budgets and with lookup-filled scaled annual budgets. Figures summarized annual scaled versus ERA5 flux magnitudes, diel flux patterns, diel source probability, site counts, ecosystem-type composition, and the spatial distribution of ERA5 annual categories.

## Results

The ERA5 half-hourly model used 131,360 observed flux records with matched ERA5 covariates. Model RMSE was 0.221 mg C m-2 per 30 min, mean absolute error was 0.0909 mg C m-2 per 30 min, and the observed-fitted correlation was 0.142. Across 39 sites, mean ERA5 annual budgets ranged from -0.257 to 1.57 g C m-2 yr-1, with a median of 0.345 g C m-2 yr-1.

ERA5 annual budgets were strongly correlated with model-standardized 30-minute annual budgets (Spearman rho = 0.948). The mean ERA5-minus-model-standardized difference was -0.0505 g C m-2 yr-1 and RMSE was 0.117 g C m-2 yr-1. Annual budget sign agreed for 37 of 39 sites. Relative to lookup-filled scaled annual categories, ERA5 annual categories agreed for 34 of 39 sites.

ERA5 annual behavior classified 32 sites as consistent sources, six as consistent sinks, and one as fluctuating. ERA5 consistent sinks included three forests, one grassland, one shrubland, and one wetland. The single fluctuating site was a forest. ERA5 consistent sources included 18 forests, eight grasslands, four shrublands, and two croplands.

Scaled and ERA5 annual budgets showed the same broad class structure but differed in magnitude. Among ERA5 consistent source sites, the median scaled annual budget was 0.122 g C m-2 yr-1, while the median ERA5 annual budget was 0.422 g C m-2 yr-1. Among ERA5 consistent sink sites, the median scaled annual budget was -0.0304 g C m-2 yr-1, while the median ERA5 annual budget was -0.157 g C m-2 yr-1. Thus, ERA5 gap filling generally preserved site ordering and annual source/sink behavior, but it amplified inferred annual source or sink magnitude for several sites.

