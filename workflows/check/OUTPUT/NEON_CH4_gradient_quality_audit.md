# NEON CH4 Gradient Quality Audit

## What This Checks
This audit separates broad retained-gradient behavior from sampled raw filter behavior. The all-site retained checks use `SITE_DATA_FILTERED_CH4.Rdata`; the raw-versus-retained checks load sampled site `*_Evaluation.Rdata` and `*_FILTER_AA_AW.Rdata` files so `dConc`, `mean_A`, `mean_B`, and dConc SNR are available.

Raw dConc filtering diagnostics were run for sampled sites: KONZ. Set `NEON_GRADIENT_AUDIT_MAX_RAW_SITES` to broaden this sample, but note that the raw evaluation files are large.

## Cross-Gradient Flags
- MBR / not flagged: n = 436904, positive FG = 49.9%, median FG = -4.83e-05, sites = 47
- AE / missing: n = 555509, positive FG = 42%, median FG = -0.000168, sites = 47
- WP / missing: n = 547149, positive FG = 42.1%, median FG = -0.00545, sites = 47

## dConc Filtering Bias
- AE: evaluable n = 253158, retained n = 35544, negative dConc evaluable = 46.7%, retained = 47.8%, retention negative/positive ratio = 1.05, retained median |dConc| = 0.00446, retained median dConcSNR = 3.43
- MBR: evaluable n = 485554, retained n = 44854, negative dConc evaluable = 46.6%, retained = 50.4%, retention negative/positive ratio = 1.16, retained median |dConc| = 0.00482, retained median dConcSNR = 3.55
- WP: evaluable n = 253158, retained n = 35564, negative dConc evaluable = 46.7%, retained = 47.9%, retention negative/positive ratio = 1.05, retained median |dConc| = 0.00446, retained median dConcSNR = 3.43

A retention negative/positive ratio above 1 means the filter retained negative/source-like CH4 gradients more readily than positive/sink-like gradients in the sampled raw files. A value below 1 means the opposite.

## Largest Retained Negative Concentration Offsets
- KONZ / MBR / 4_2: n = 6689, negative dConc = 58.4%, median dConc = -0.00348, median |dConc| = 0.00683, median SNR = 4.41
- KONZ / MBR / 3_1: n = 8048, negative dConc = 56.4%, median dConc = -0.00331, median |dConc| = 0.00658, median SNR = 4.31
- KONZ / MBR / 4_1: n = 1454, negative dConc = 53.7%, median dConc = -0.00313, median |dConc| = 0.00641, median SNR = 3.73
- KONZ / MBR / 3_2: n = 453, negative dConc = 53.4%, median dConc = -0.00287, median |dConc| = 0.0068, median SNR = 3.7
- KONZ / MBR / 4_3: n = 468, negative dConc = 53.2%, median dConc = -0.0028, median |dConc| = 0.00676, median SNR = 3.73
- KONZ / MBR / 2_1: n = 585, negative dConc = 50.6%, median dConc = -0.00204, median |dConc| = 0.00668, median SNR = 3.84
- KONZ / WP / 4_2: n = 9317, negative dConc = 48.5%, median dConc = 0.00234, median |dConc| = 0.0059, median SNR = 4.18
- KONZ / AE / 4_2: n = 9300, negative dConc = 48.5%, median dConc = 0.00235, median |dConc| = 0.00588, median SNR = 4.18
- KONZ / AE / 3_2: n = 925, negative dConc = 47.1%, median dConc = 0.00267, median |dConc| = 0.0056, median SNR = 3.59
- KONZ / WP / 3_2: n = 930, negative dConc = 46.8%, median dConc = 0.00271, median |dConc| = 0.00557, median SNR = 3.59

## Outputs
- `OUTPUT/NEON_CH4_retained_gradient_quality_by_site.csv`
- `OUTPUT/NEON_CH4_cross_gradient_flag_summary.csv`
- `OUTPUT/NEON_CH4_sensor_pair_behavior_summary.csv`
- `OUTPUT/NEON_CH4_dConc_filter_bias_by_pair.csv`
- `OUTPUT/NEON_CH4_dConc_filter_bias_overall.csv`
- `OUTPUT/NEON_CH4_concentration_offset_summary.csv`
- `OUTPUT/NEON_CH4_gradient_quality_raw_load_issues.csv`
- `FIGURES/NEON_CH4_cross_gradient_flag_audit.png`
- `FIGURES/NEON_CH4_dConc_filter_bias_audit.png`
- `FIGURES/NEON_CH4_sensor_pair_behavior_audit.png`
