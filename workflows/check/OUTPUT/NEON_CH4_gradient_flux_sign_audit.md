# NEON CH4 Gradient-Flux Sign Convention Audit

## Interpretation
`dConc` is `concentration_A - concentration_B`. For rows where tower position A is above B, Fickian upward-positive flux should have `sign(FG_mean) = -sign(dConc)`: higher CH4 aloft than below implies downward uptake, and higher CH4 below than aloft implies upward emission.

The filtered CH4 data overwhelmingly follow that opposite-sign relationship for A-above-B rows. That means the current `FG_mean` sign convention is consistent with upward-positive flux rather than obviously inverted.

The positive annual/source behavior is therefore not explained by a simple sign-convention inversion in `FG_mean`. It is driven by many retained observations having negative CH4 gradients (`dConc < 0`), which translate to positive upward flux under the current convention.

## Overall Sign Relationship
- AE / A above B: n = 119951, opposite sign = 100%, same sign = 0%, median FG = -0.000119, median dConc = 0.00261
- MBR / A above B: n = 103165, opposite sign = 78.3%, same sign = 21.7%, median FG = 0.000931, median dConc = 0.0022
- WP / A above B: n = 118222, opposite sign = 100%, same sign = 0%, median FG = -0.00361, median dConc = 0.00259

## Total-Flux Sites With Highest Positive Fraction
- NOGP: n = 9164, positive total flux = 81%, median total = 0.461, median gradient = 0.462, median storage = -0.000229
- MOAB: n = 10008, positive total flux = 77.4%, median total = 0.219, median gradient = 0.219, median storage = -9.19e-05
- LENO: n = 350, positive total flux = 76.3%, median total = 0.16, median gradient = 0.151, median storage = 0.00204
- OAES: n = 9494, positive total flux = 76.1%, median total = 0.447, median gradient = 0.447, median storage = -0.000351
- RMNP: n = 128, positive total flux = 75%, median total = 0.103, median gradient = 0.102, median storage = 0.000628
- JORN: n = 3904, positive total flux = 73.4%, median total = 0.145, median gradient = 0.144, median storage = -0.000252
- STER: n = 8521, positive total flux = 73.3%, median total = 0.233, median gradient = 0.233, median storage = -9.61e-05
- KONZ: n = 16097, positive total flux = 72%, median total = 0.294, median gradient = 0.293, median storage = -1.31e-05
- SCBI: n = 2215, positive total flux = 70.4%, median total = 0.0154, median gradient = 0.0153, median storage = 0.000549
- GRSM: n = 2700, positive total flux = 70.3%, median total = 0.265, median gradient = 0.263, median storage = 0.000331

## Outputs
- `OUTPUT/NEON_CH4_gradient_flux_sign_audit_by_site.csv`
- `OUTPUT/NEON_CH4_gradient_flux_sign_audit_overall.csv`
- `OUTPUT/NEON_total_flux_sign_summary_by_site.csv`
- `FIGURES/NEON_CH4_gradient_flux_sign_audit.png`
