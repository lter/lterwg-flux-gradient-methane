# RSHP: ####
library( ggplot2)

source(fs::path(DirRepo.ch4, 'functions/calc_CCC_SHP.R'))

load(fs::path(localdir.ch4, paste0("/SITEval_DATA_FILTERED_CH4.Rdata")) )

canopy.info <- read.csv( file.path(paste(localdir, "Val_canopy.csv", sep="/"))) %>% rename(Site = SITE_ID)

# Use the RSHP model to Filter THE DATA: #####

# Create Datasets with only the data in it from RSHP: 
site.list <- SITEval_DATA_FILTERED %>% names()

# Compare FG to EC: ####
source(fs::path(DirRepo.eval, 'functions/calc.One2One.CCC_testing.R'))

# Canopy-structure-based candidate pair restriction (paper Sect. X):
#
#   Short canopy (hc ≲ 5 m): lower sensor at canopy top ±2 m, upper 5–10 m higher.
#   Tall canopy  (hc ≳ 15 m): AW config — zA in trunk space (hc−15 to −20 m),
#                              zB above canopy (hc+10 to +15 m).
#
# SE-Sto  (hc ≈ 0.25 m, tower = 3 m): 5–10 m separation is physically impossible;
#           the Canopy_L1 AA/AW filter is sufficient — all 9 AW+AA pairs retained.
# SE-Svb  (hc ≈ 23.5 m, 16 levels): ideal AW pairs (lower within canopy) lack CH4
#           concentration measurements; only pair 13_8 (85/35 m, AA) has data.
# US-Uaf  (hc ≈ 5 m, levels: 1, 2, 4, 8 m): Level 4 (8 m) as upper sensor gives
#           4–7 m separations — the closest achievable to the 5–10 m rule.
#           Drops small-separation pairs 2_1 (1 m), 3_1 (3 m), 3_2 (2 m).
#
# NOTE: CHM in Val_canopy.csv appears incorrect for SE-Svb (3 instead of ~23.5 m)
# and US-Uaf (18 instead of ~5 m); Canopy_L1 labels were assigned from the actual
# canopy heights and are correct — only the CHM column is wrong.
canopy_struct_allowed <- tribble(
  ~Site,    ~dLevelsAminusB,
  # SE-Sto — tower height precludes 5-10 m separation; keep all AA/AW pairs
  "SE-Sto", "3_1",
  "SE-Sto", "3_2",
  "SE-Sto", "4_1",
  "SE-Sto", "4_2",
  "SE-Sto", "4_3",
  "SE-Sto", "5_1",
  "SE-Sto", "5_2",
  "SE-Sto", "5_3",
  "SE-Sto", "5_4",
  # SE-Svb — only pair with CH4 measurements (ideal AW pairs lack data)
  "SE-Svb", "13_8",
  # US-Uaf — Level 4 (8 m) as upper sensor; excludes small-separation pairs
  "US-Uaf", "4_1",
  "US-Uaf", "4_2",
  "US-Uaf", "4_3"
)

# Add canopy information and apply both the AA/AW filter and the
# canopy-structure candidate restriction:
SITEval_DATA_FILTEREDc <- list()
for( site in site.list){
  print(site)
  SITEval_DATA_FILTEREDc[[site]] <- SITEval_DATA_FILTERED[[site]] %>%
    mutate(Site = site) %>%
    left_join(canopy.info, by = c('Site', 'dLevelsAminusB')) %>%
    filter(Canopy_L1 != "WW", Canopy_L1 != "WA") %>%
    inner_join(
      canopy_struct_allowed %>% filter(Site == site) %>%
        dplyr::select(dLevelsAminusB),
      by = 'dLevelsAminusB'
    )
}

# Data now only contains AW and AA:

CCC_VAL <- ccc.val(sites.tibble = SITEval_DATA_FILTEREDc)

CCC_VAL %>% ggplot(aes(x = CCC, y = Site, col = Approach)) +
  geom_point() + theme_bw() +
  facet_wrap(~ Season + gas, nrow = 3) + xlim(-1, 1)

CCC_VAL$CCC %>% summary

# Save the files:
save(CCC_VAL, file = fs::path(localdir.ch4, "/CCC_CH4.Rdata"))

# Add RSHP information to the fluxes — tiered filter: ####
#
# For validation towers the RSHP is built directly from CH4 CCC (no CO2/H2O
# reference as in NEON) using three tiers per site × gas × season:
#
#   Tier 1 — standard   : CCC > 0.1  (well-performing pairs)
#   Tier 2 — best_avail : no Tier-1 pair exists; use all pairs with CCC > 0
#   Tier 3 — unreliable : no pair has CCC > 0; use the single best-CCC pair
#                         and flag it so downstream analyses can caveat it
#
# Fix: Season is now derived from the observation timestamp BEFORE the CCC
# join, and the join includes Season as a key.  The previous Season-free join
# created 4× duplicates (one per seasonal CCC row) which inflated row counts.

CCC_THRESHOLD <- 0.1   # Tier-1 threshold

# Build a compact selection table: for each site×gas×season, decide which
# pairs pass and assign their tier.
rshp_pairs <- CCC_VAL %>%
  filter(gas == "CH4") %>%
  group_by(Site, gas, Season) %>%
  mutate(
    tier = case_when(
      CCC > CCC_THRESHOLD ~ "standard",
      CCC > 0             ~ "best_avail",
      TRUE                ~ "unreliable"
    ),
    best_CCC = max(CCC, na.rm = TRUE)
  ) %>%
  # Keep Tier-1 pairs if any exist; otherwise keep Tier-2; otherwise the
  # single best-CCC pair for Tier-3.
  mutate(
    has_tier1 = any(tier == "standard"),
    has_tier2 = any(tier == "best_avail")
  ) %>%
  filter(
    (has_tier1  &  tier == "standard")   |
    (!has_tier1 &  has_tier2  &  tier == "best_avail") |
    (!has_tier1 & !has_tier2  &  CCC == best_CCC)
  ) %>%
  mutate(
    RSHP_tier = case_when(
      has_tier1  ~ "standard",
      has_tier2  ~ "best_avail",
      TRUE       ~ "unreliable"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(Site, gas, Season, Approach, dLevelsAminusB, CCC, RSHP_tier)

message("RSHP pair selection summary:")
rshp_pairs %>%
  count(Site, Season, RSHP_tier) %>%
  print(n = 50)

# Save the selection table for reporting
save(rshp_pairs,
     file = fs::path(localdir.ch4, "/VAL_RSHP_pairs.Rdata"))

SITEval_DATA_FILTERED_RSHPc <- list()

for (site in site.list) {

  # Derive season from timestamp first, then join CCC on matching season
  SITEval_DATA_FILTERED_RSHPc[[site]] <- SITEval_DATA_FILTEREDc[[site]] %>%
    mutate(
      month  = format(timeEndA.local, '%m') %>% as.numeric(),
      Season = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5)  ~ "Spring",
        month %in% c(6, 7, 8)  ~ "Summer",
        TRUE                   ~ "Autumn"),
      hour   = format(timeEndA.local, '%H')
    ) %>%
    # Season-aware join — no more 4× duplication
    inner_join(
      rshp_pairs %>% filter(Site == site) %>%
        dplyr::select(gas, Season, Approach, dLevelsAminusB, CCC, RSHP_tier),
      by = c('gas', 'Season', 'Approach', 'dLevelsAminusB')
    )
}

# For validation sites there is no separate RSHP model — both paths use the
# tiered CCC filter above.
SITEval_DATA_FILTERED_RSHP <- SITEval_DATA_FILTERED_RSHPc

#Save files:
  save( SITEval_DATA_FILTERED_RSHPc, SITEval_DATA_FILTERED_RSHP,
        file= fs::path(localdir.ch4, paste0("/SITEval_DATA_FILTERED_RSHP.Rdata")))

# Harmonize data and test result:  ####

# Explore 2 options: CCC > 0 and RSHP
harmonize_val <- function( tibble){
  
  site.list <- names(tibble )
  harmonized <- list()
  for( i in site.list){
    print(i)
    
    df <- tibble[[i]]
    
    harmonized[[i]] <-   df %>% 
      mutate( time.rounded = round_date(timeEndA.local, unit = "30 minutes") ) %>% 
      reframe(.by= c(time.rounded, gas),
                                        FG_mean = median(FG_mean, na.rm=T),
                                        EC_mean = mean(EC_mean, na.rm=T),
                                        Tair_C = mean(Tair_C, na.rm=T),
                                        PAR = mean(PAR, na.rm=T))  %>% 
      mutate( Month = format(time.rounded , '%m') %>% as.numeric) %>% mutate(
        Season = case_when(
          Month %in% c(12, 1, 2) ~ "Winter",
          Month %in% 3:5 ~ "Spring",
          Month %in% 6:8 ~ "Summer",
          Month %in% 9:11 ~ "Autumn")) 
    
  }
  return(harmonized )
  
}

# Both objects use the tiered CCC filter; harmonise both for downstream compatibility
SITEval_DATA_FILTERED_RSHP_H  <- harmonize_val(tibble = SITEval_DATA_FILTERED_RSHP)
SITEval_DATA_FILTERED_RSHPc_H <- harmonize_val(tibble = SITEval_DATA_FILTERED_RSHPc)

#Save files: 
save( SITEval_DATA_FILTERED_RSHP_H, SITEval_DATA_FILTERED_RSHPc_H,
      file= fs::path(localdir.ch4, paste0("/SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata")))

# CCC for Harmonized data:

harmonized_CCC <- function( tibble){
  
  CCCval_Harmonized <- data.frame()
  
  for( i in names( tibble)){
    print(i)
    
    seasons <- tibble[[i]] $Season %>% unique
    
    for(s in seasons){
      
      data <- tibble[[i]] %>% filter( Season == s) 
      GAS <- data$gas %>% unique
      
      for( g in GAS ){
        
        df.summary <-linear_terms2('FG_mean', 'EC_mean' , df = tibble[[i]] %>% filter(Season == s, gas == g), site=i) %>% mutate( Season = s, gas = g)
        
        CCCval_Harmonized  <- rbind(CCCval_Harmonized, df.summary)
      }
    }}
  
  return(CCCval_Harmonized )
}


CCC_RSHP_H <- harmonized_CCC(tibble=SITEval_DATA_FILTERED_RSHP_H )
CCC_RSHPc_H <- harmonized_CCC(tibble= SITEval_DATA_FILTERED_RSHPc_H )

# Data VIZ:
CCC_RSHP_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)

CCC_RSHPc_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)