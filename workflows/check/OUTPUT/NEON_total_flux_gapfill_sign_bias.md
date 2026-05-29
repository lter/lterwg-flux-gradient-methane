# NEON Total-Flux Gapfill Sign Bias Check

## Interpretation
Gapfilling slightly shifts the annual budgets in the source-positive direction. The shift is modest relative to the fact that most observed-only scaled annual budgets are already positive, but it can change sign in low-to-moderate coverage site-years.

The sign conclusion should therefore be reported with coverage diagnostics: gapfilling is not the sole reason the budgets are source-like, but it does increase the number of source-classified site-years.

## Summary
- Interpretable site-years: n = 46, observed source = 41, gapfilled source = 43, sign flips = 4 (sink->source 3, source->sink 1), median budget shift = 0.0667 g C m-2 yr-1, median interval positive-fraction shift = 0.119
- All other observed site-years: n = 74, observed source = 60, gapfilled source = 61, sign flips = 15 (sink->source 8, source->sink 7), median budget shift = 0.0561 g C m-2 yr-1, median interval positive-fraction shift = 0.152

## Interpretable Site-Years With Sign Flips
- DCFS 2021: observed = 0.135, gapfilled = -0.0866, coverage = 5.3%, observed source -> gapfilled sink
- KONZ 2022: observed = -0.0785, gapfilled = 0.131, coverage = 33.0%, observed sink -> gapfilled source
- SERC 2022: observed = -0.0427, gapfilled = 0.139, coverage = 5.1%, observed sink -> gapfilled source
- STEI 2022: observed = -0.145, gapfilled = 0.00234, coverage = 6.7%, observed sink -> gapfilled source

## Outputs
- `OUTPUT/NEON_total_flux_gapfill_sign_bias_diagnostic.csv`
- `OUTPUT/NEON_total_flux_gapfill_sign_bias_summary.csv`
- `FIGURES/NEON_total_flux_gapfill_sign_bias.png`
