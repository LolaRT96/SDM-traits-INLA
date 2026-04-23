# ==============================================================================
# Title: Data Preparation for Empirical SDM – Merluccius merluccius with FishBase Traits
# Author: M. Grazia Pennino
# Date:   2025-07-08
# Description:
#   - Loads ICES DATRAS haul header (HH) and catch-at-age (CA) data
#   - Builds presence/absence for Merluccius merluccius across all hauls
#   - Calculates per‐haul mean length from CA length classes
#   - Joins environmental covariates from HH
#   - Adds selected functional traits from FishBase
#   - Saves ready‐to‐model dataset
# ==============================================================================

library(dplyr)
library(icesDatras)
library(rfishbase)   # For FishBase trait data

# --- 1. Load DATRAS data ------------------------------------------------------
hh_all <- readRDS("C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/data/hh_all.rds")    # all hauls, all species
ca_all <- readRDS("C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/data/ca_all.rds")    # all catch records

# --- 2. Identify Merluccius merluccius records --------------------------------
# Use c = 126484 for M. merluccius
target_spec <- 126484 # M. merluccius
# target_spec <- 126439 #micromessistius poutaasou
# target_spec <- 127146 #L. boscii

# presence_df: one row per haul where merluccius was caught
presence_df <- ca_all %>%
  filter(SpecCode == target_spec) %>%
  distinct(Year, Survey, StNo, HaulNo) %>%
  mutate(presence = 1L)

# --- 2.5. Compute haul‐level mean length from CA -------------------------------

# CA has columns LngtClass (length bin, cm) and CANoAtLngt (# fish at that length)

mean_length_df <- ca_all %>%
  filter(SpecCode == target_spec, !is.na(LngtClass), CANoAtLngt > 0) %>%
  group_by(Year, Survey, StNo, HaulNo) %>%
  summarise(
    mean_length_mm = sum(LngtClass * CANoAtLngt, na.rm = TRUE) /
      sum(CANoAtLngt, na.rm = TRUE),
    .groups = "drop"
  )

# --- 3. Build full haul list & assign 0/1 -------------------------------------

all_hauls <- hh_all %>%
  distinct(Year, Survey, StNo, HaulNo)

sdm_base <- all_hauls %>%
  left_join(
    presence_df %>% mutate(Year = as.integer(Year)), 
    by = c("Year","Survey","StNo","HaulNo")
  ) %>%
  mutate(presence = if_else(is.na(presence), 0L, 1L))

glimpse(sdm_base)

# --- 4. Attach environmental covariates from HH -------------------------------

covars <- hh_all %>%
  select(Year, Survey, StNo, HaulNo,
         ShootLat, ShootLong,
         Depth,            # bottom depth
         BotTemp,          # bottom temperature
         BotSal) %>%       # bottom salinity
  distinct()

# Definimos las columnas clave
key_cols <- c("Year", "Survey", "StNo", "HaulNo")

sdm_data <- sdm_base %>%
  # Antes del join convertimos solo las columnas clave
  left_join(
    covars %>% mutate(
      Year = as.integer(Year),
      across(c("Survey","StNo","HaulNo"), as.character)
    ),
    by = key_cols
  ) %>%
  left_join(
    mean_length_df %>% mutate(
      Year = as.integer(Year),
      across(c("Survey","StNo","HaulNo"), as.character)
    ),
    by = key_cols
  ) %>%
  filter(!is.na(Depth))

# sdm_data <- sdm_base %>%
#   left_join(covars,       by = c("Year","Survey","StNo","HaulNo")) %>%
#   left_join(mean_length_df, by = c("Year","Survey","StNo","HaulNo")) %>%
#   filter(!is.na(Depth))   # drop any hauls lacking core covariates

# --- 5. Add FishBase functional traits ----------------------------------------

fb_traits <- species("Merluccius merluccius",
                     fields = c("Species", "Length", "Weight",
                                "LongevityWild", "Vulnerability",
                                "DepthRangeShallow", "DepthRangeDeep",
                                "DemersPelag"))

sdm_data <- sdm_data %>%
  mutate(
    FB_max_length_cm   = fb_traits$Length[1],       # FishBase maximum length
    Weight_g           = fb_traits$Weight[1],
    LongevityWild_yrs  = fb_traits$LongevityWild[1],
    Vulnerability_idx  = fb_traits$Vulnerability[1],
    DepthRangeMin_m    = fb_traits$DepthRangeShallow[1],
    DepthRangeMax_m    = fb_traits$DepthRangeDeep[1],
    DemersPelag        = fb_traits$DemersPelag[1]
  )

# Note: `mean_length_cm` is your *observed* average on each haul,
#       `FB_max_length_cm` is the species‐level trait from FishBase.

# --- 6. Save final dataset ----------------------------------------------------
saveRDS(sdm_data, file = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/data/sdm_data_merluc.rds")

# --- 7. Quick summary ----------------------------------------------------------
message("Final SDM dataset for M. merluccius:")
glimpse(sdm_data)
print(table(sdm_data$presence))

