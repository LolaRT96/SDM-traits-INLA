
# -----------------------------
# Title: Download and Harmonize ICES DATRAS HH and CA Data
# Authors: M.Grazia Pennino & M. D. Riesgo
# Date: 2025-07-04
#
# Description:
# This script downloads haul header (HH) and catch-at-age (CA) data from the ICES DATRAS database
# for selected surveys and years, harmonizes column types to avoid binding issues,
# and saves the combined datasets as RDS files for downstream analyses.
# -----------------------------

# Create data directory if it does not exist
if (!dir.exists("data")) {
  dir.create("data")
}

# -----------------------------
# STEP 0. Load required packages
# -----------------------------
library(dplyr)
library(purrr)
library(icesDatras)
library(icesVocab)

# -----------------------------
# STEP 1. Define surveys and years
# -----------------------------
surveys <- c("SP-NORTH", "PT-IBTS")  
# Additional surveys like "SP-ARSA" can be added later
years <- 2010:2022
quarter <- 1:4  # All quarters

# -----------------------------
# STEP 2. Safe HH downloader with harmonized column types
# -----------------------------

safe_get_HH_year <- possibly(function(survey, yr) {
  df <- getDATRAS("HH", survey, year = yr, quarters = quarter)
  
  # Harmonize ID columns to character to avoid conflicts when binding data frames
  id_cols <- c("StNo", "HaulNo", "Gear", "DoorType", "GearEx",
               "Ship", "StatRec", "LngtCode", "DepthStratum", "HydroStNo")
  df <- df %>%
    mutate(across(any_of(id_cols), as.character))
  
  return(df)
}, otherwise = data.frame())

# -----------------------------
# STEP 3. Download and combine all HH data
# -----------------------------

hh_all <- map_df(surveys, function(survey) {
  map_df(years, function(yr) {
    message("Downloading HH data for ", survey, " - Year ", yr)
    safe_get_HH_year(survey, yr)
  })
})

# -----------------------------
# STEP 4. Safe CA downloader with harmonized column types
# -----------------------------
safe_get_CA_year <- possibly(function(survey, yr) {
  df <- getDATRAS("CA", survey, year = yr, quarters = quarter)
  
  # Convert specific columns to character to avoid type mismatch issues
  df <- df %>%
    mutate(across(everything(), as.character))
  
  return(df)
}, otherwise = data.frame())

# -----------------------------
# STEP 5. Download and combine all CA data
# -----------------------------

ca_all <- map_df(surveys, function(survey) {
  map_df(years, function(yr) {
    message("Downloading CA data for ", survey, " - Year ", yr)
    safe_get_CA_year(survey, yr)
  })
})

# -----------------------------
# STEP 6. Save combined datasets to disk
# -----------------------------
# saveRDS(hh_all, file = "~/SDMs_Traits/data/hh_all.rds")
# saveRDS(ca_all, file = "~/SDMs_Traits/data/ca_all.rds")

message("✅ DATRAS HH and CA data successfully downloaded and saved.")

summary(hh_all)
summary(ca_all)

# -----------------------------
# STEP 7. Identify empirical case study species
# -----------------------------

ca_all <- ca_all %>%
  mutate(
    CANoAtLngt = as.numeric(CANoAtLngt),
    IndWgt     = as.numeric(IndWgt),
    LngtClass  = as.numeric(LngtClass)
  )


# Find all species (SpecCode) with ≥1 catch record in every year
species_full_years <- ca_all %>%
  filter(Year %in% years) %>%
  group_by(SpecCode) %>%
  summarize(
    years_present = n_distinct(Year),
    total_catch   = sum(CANoAtLngt, na.rm = TRUE),
    records       = n()
  ) %>%
  ungroup() %>%
  filter(years_present == length(years)) %>%
  arrange(desc(total_catch))

# Show top 10 species based on total catch
head(species_full_years, 10)

# -----------------------------
# STEP 8. Translate SpecCode to scientific names using icesVocab
# -----------------------------

# Define a function that maps ICES SpecCode to scientific names using the ICES SpecWoRMS vocabulary

translate_spec_code_ices <- function(spec_codes) {
  species_vocab <- icesVocab::getCodeList("SpecWoRMS")
  
  translated <- species_vocab %>%
    filter(Key %in% as.character(spec_codes)) %>%
    rename(SpecCode = Key,
           scientific_name = Description) %>%
    select(SpecCode, scientific_name)
  
  return(translated)
}

# Apply the function to top species

top_spec_codes <- head(species_full_years$SpecCode, 10)
translated_species <- translate_spec_code_ices(top_spec_codes)

# Print the translated table
print(translated_species)

# -----------------------------
# (Optional) Join back to original table and export
# -----------------------------

species_full_years$SpecCode <- as.character(species_full_years$SpecCode)
translated_species$SpecCode <- as.character(translated_species$SpecCode)

species_full_years_named <- species_full_years %>%
  left_join(translated_species, by = "SpecCode")

# View with names
head(species_full_years_named)

# Save to file
# write.csv(species_full_years_named, 
#           ~/SDMs_Traits/data/species_full_years_named.csv", 
#           row.names = FALSE)


