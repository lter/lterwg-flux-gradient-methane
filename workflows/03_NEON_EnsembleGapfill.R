# Ensemble NEON CH4 gradient fluxes:
# - Sites with CO2 RSHPs: equal-weight mean of the CH4 fluxes from the pairs that are reliable for CO2
#   (RSHP.CH4 from step 02).
# - Sites without a CO2 RSHP (CH4.fallback.sites): no-EC rule of essd-2026-413 Supplement S2. Reference = mean of
#   the above-canopy (AA) MBR pairs in each half-hour; AE/WP pairs with CCC >= 0.5 against it (>= 30 shared
#   half-hours) are rescaled to it by reduced-major-axis, F' = (F - a) / b, and averaged with the AA MBR fluxes.
library(tidyverse)

# Import the RSHP information (NEON_RSHP_EVAL, CH4.fallback.sites):
load(fs::path(localdir.ch4 ,paste0("NEON_RSHP_CCC.Rdata")) )
if (!exists("CH4.fallback.sites")) CH4.fallback.sites <- setdiff(site.list, NEON_RSHP_EVAL$Site[NEON_RSHP_EVAL$RSHP.CH4 == '1'])

# Import the filtered data:
load(fs::path(localdir.ch4, paste0("/NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata")))

ccc.min <- 0.5   # CCC against the AA MBR reference (fallback sites)
n.min   <- 30    # shared half-hours needed to screen a pair

lccc <- function(x, y) { mx <- mean(x); my <- mean(y)
  2 * mean((x - mx) * (y - my)) / (mean((x - mx)^2) + mean((y - my)^2) + (mx - my)^2) }

# No-EC (S2) selection for one site: returns the selected rows with FG_mean rescaled, plus a pair table
s2_select <- function(d) {
  d   <- d %>% filter(!is.na(FG_mean))
  aa  <- d %>% filter(Approach == "MBR", Canopy_L1 == "AA")
  if (!nrow(aa)) return(list(data = d[0, ], pairs = tibble()))
  REF <- aa %>% reframe(.by = time.rounded, REF = mean(FG_mean))
  out <- list(aa); info <- list()
  cand <- d %>% filter(Approach != "MBR") %>% distinct(Approach, dLevelsAminusB)
  for (i in seq_len(nrow(cand))) {
    x <- d %>% filter(Approach == cand$Approach[i], dLevelsAminusB == cand$dLevelsAminusB[i])
    G <- x %>% reframe(.by = time.rounded, GF = mean(FG_mean)) %>% inner_join(REF, by = "time.rounded")
    cr <- if (nrow(G) >= n.min) lccc(G$GF, G$REF) else NA_real_
    b  <- if (nrow(G) >= 3) sign(cor(G$GF, G$REF)) * sd(G$GF) / sd(G$REF) else NA_real_
    a  <- if (nrow(G) >= 3) mean(G$GF) - b * mean(G$REF) else NA_real_
    keep <- is.finite(cr) && cr >= ccc.min && is.finite(b) && b != 0
    info[[i]] <- tibble(Approach = cand$Approach[i], dLevelsAminusB = cand$dLevelsAminusB[i],
                        n_ref = nrow(G), CCC_ref = cr, scale_b = b, scale_a = a, selected = keep)
    if (keep) out[[length(out) + 1]] <- x %>% mutate(FG_mean = (FG_mean - a) / b)
  }
  list(data = bind_rows(out), pairs = bind_rows(info))
}

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE <- list()
CH4_PAIR_SELECTION <- list()

for( site in site.list){
  print( site)
  if (is.null(SITE_DATA_FILTERED[[site]])) { message("  no filtered CH4 data"); next }

  d <- SITE_DATA_FILTERED[[site]] %>% filter(gas == "CH4") %>% distinct() %>%
    mutate( time.rounded = timeEndA.local %>% round_date( unit = "30 minutes") )

  if (site %in% CH4.fallback.sites) {
    s2 <- s2_select(d)
    sel <- s2$data
    rule <- "S2: AA MBR reference + AE/WP screened and rescaled (no CO2 RSHP)"
    CH4_PAIR_SELECTION[[site]] <- bind_rows(
      sel %>% filter(Approach == "MBR") %>% distinct(Approach, dLevelsAminusB) %>% mutate(selected = TRUE),
      s2$pairs) %>% mutate(Site = site, rule = rule)
  } else {
    pairs <- NEON_RSHP_EVAL %>% filter( Site == site, RSHP.CH4 == '1') %>%
      mutate(Approach = as.character(Approach)) %>% distinct(Approach, dLevelsAminusB)
    sel <- d %>% mutate(Approach = as.character(Approach)) %>% semi_join(pairs, by = c('Approach', 'dLevelsAminusB'))
    rule <- "CO2 RSHP pairs"
    CH4_PAIR_SELECTION[[site]] <- pairs %>% mutate(selected = TRUE, Site = site, rule = rule)
  }

  ens <- sel %>% filter(!is.na(FG_mean)) %>%
    reframe( .by = time.rounded,
             FG_ENSEMBLE_RSHP = mean(FG_mean),
             FG_ENSEMBLE_NPAIRS = n_distinct(paste(Approach, dLevelsAminusB)))

  # All half-hours with CH4 gradient data; the ensemble is NA where no selected pair has a value
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- d %>%
    mutate( count = case_when( is.na(FG_mean) == FALSE ~ 1, TRUE ~ 0)) %>%
    reframe( .by= c(time.rounded),
             count = sum(count),
             PAR = mean(PAR, na.rm=T),
             Tair_C = mean(Tair_C, na.rm=T)) %>%
    left_join(ens, by = "time.rounded") %>%
    mutate(CH4_PAIR_RULE = rule,
           month = format(time.rounded,'%m') %>% as.numeric,
           season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(time.rounded,'%H'))
}

CH4_PAIR_SELECTION <- bind_rows(CH4_PAIR_SELECTION)
print(CH4_PAIR_SELECTION %>% filter(selected) %>% count(rule, name = "pairs"))

# ENSEMBLE:
save( SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE, CH4_PAIR_SELECTION,
      file=fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")))
write.csv(CH4_PAIR_SELECTION, fs::path(localdir.ch4, "CH4_PAIR_SELECTION.csv"), row.names = FALSE)

if (!exists("upload.drive") || isTRUE(upload.drive))
  googledrive::drive_upload(media = fs::path(localdir.ch4, paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")), overwrite = T, path = drive_url)
