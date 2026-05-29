# NEON Driver Scale Analysis

## Design
- Source/sink classes are derived in this script from balanced 30-minute months, sampling-adjusted daily fluxes, and ERA5 annual budgets.
- Source probability and flux magnitude are modeled separately at 30-minute, daily, and annual scales.
- 30-minute and daily temporal drivers are centered within site before modeling. This estimates within-site temporal effects separately from among-site differences.
- Annual models use ERA5 annual anomalies within site plus ERA5/site mean and soil/ecosystem attributes. They are intentionally simpler linear models, with site-block bootstrap handling repeated-site dependence.
- Variable importance uncertainty uses site-block bootstrap weights, not row-level resampling.
- Bootstrap replicates: 1.
- For final inference, rerun with a larger bootstrap count, for example `NEON_DRIVER_BOOT_N=100 Rscript workflows/NEON.DriverScale.Analysis.R`.

## Model-Ready Tables
- 30-minute: 131360 rows, 39 sites.
- Daily: 18242 rows, 39 sites.
- Annual ERA5: 150 rows, 39 sites.

## Outputs
- `OUTPUT/30min_site_behavior.csv`
- `OUTPUT/30min_halfhour_balancing_bins.csv`
- `OUTPUT/30min_total_vs_gradient_behavior_comparison.csv`
- `OUTPUT/NEON_driver_scale_source_sink_summary.csv`
- `OUTPUT/NEON_driver_scale_source_sink_counts.csv`
- `OUTPUT/NEON_driver_scale_class_change_30min_to_daily.csv`
- `OUTPUT/NEON_driver_scale_class_change_daily_to_annual.csv`
- `OUTPUT/NEON_driver_scale_class_change_30min_to_annual.csv`
- `OUTPUT/NEON_driver_scale_30min_driver_data.csv`
- `OUTPUT/NEON_driver_scale_daily_driver_data.csv`
- `OUTPUT/NEON_driver_scale_annual_era5_driver_data.csv`
- `OUTPUT/NEON_driver_scale_model_summary.csv`
- `OUTPUT/NEON_driver_scale_variable_importance_bootstrap.csv`
- `OUTPUT/NEON_driver_scale_variable_importance_summary.csv`
- `OUTPUT/NEON_driver_scale_partial_effects.csv`
- `OUTPUT/NEON_driver_scale_cross_scale_summary.csv`
- `FIGURES/NEON_driver_scale_source_sink_counts.png`
- `FIGURES/NEON_driver_scale_source_sink_transitions.png`
- `FIGURES/NEON_driver_scale_variable_importance.png`
- `FIGURES/NEON_driver_scale_partial_effects.png`
