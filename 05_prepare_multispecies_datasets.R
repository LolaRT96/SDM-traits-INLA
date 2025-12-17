# ======================================================================
# Title: Multispecies Spatial SDM with Trait-Modulated Response (INLA)
# Description: Prepares multi-species ICES DATRAS data using SpecCode translation,
#              calculates mean observed length per haul and species,
#              and prepares data for hierarchical SDM modeling.
# ======================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(INLA)
library(icesVocab)
library(sp)
library(tibble)

# --- 1. Load DATRAS data -------------------------------------------------------

hh_all <- readRDS("data/hh_all.rds")
ca_all <- readRDS("data/ca_all.rds")

# --- 2. Identify top species by total catch ------------------------------------

top_species_table <- ca_all %>%
  group_by(SpecCode) %>%
  summarise(total_catch = sum(CANoAtLngt, na.rm = TRUE)) %>%
  slice_max(total_catch, n = 10)


top_species_table <- top_species_table %>%
  mutate(SpecCode = as.character(SpecCode))

# --- 3. Translate SpecCode to scientific name using ICES vocabulary ------------

translate_spec_code_ices <- function(spec_codes) {
  species_vocab <- icesVocab::getCodeList("SpecWoRMS")
  
  # Asegurar que los códigos de entrada son character
  spec_codes <- as.character(spec_codes)
  
  translated <- species_vocab %>%
    mutate(Key = as.character(Key)) %>%        # <- asegurar tipo
    filter(Key %in% spec_codes) %>%
    rename(SpecCode = Key, scientific_name = Description) %>%
    mutate(SpecCode = as.character(SpecCode)) %>%   # <- asegurar salida character
    dplyr::select(SpecCode, scientific_name)
  
  return(translated)
}

translated_species <- translate_spec_code_ices(top_species_table$SpecCode)

# Merge names
top_species_named <- top_species_table %>%
  left_join(translated_species, by = "SpecCode") %>%
  filter(!is.na(scientific_name))  # keep only those with a name

# --- 4. Subset data for selected species ---------------------------------------

ca_sel <- ca_all %>%
  filter(SpecCode %in% top_species_named$SpecCode)

hauls_sel <- ca_sel %>%
  distinct(Year, Survey, StNo, HaulNo)

hh_sel <- hh_all %>%
  semi_join(hauls_sel, by = c("Year", "Survey", "StNo", "HaulNo"))

# --- 5. Build presence-absence table -------------------------------------------

presence_df <- ca_sel %>%
  mutate(scientific_name = top_species_named$scientific_name[match(SpecCode, top_species_named$SpecCode)]) %>%
  group_by(Year, Survey, StNo, HaulNo, scientific_name) %>%
  summarise(presence = 1, .groups = "drop")

# All combinations of hauls x species
all_combinations <- expand.grid(
  Year = unique(hh_sel$Year),
  Survey = unique(hh_sel$Survey),
  StNo = unique(hh_sel$StNo),
  HaulNo = unique(hh_sel$HaulNo),
  scientific_name = unique(top_species_named$scientific_name),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

all_pa <- all_combinations %>%
  left_join(presence_df, by = c("Year", "Survey", "StNo", "HaulNo", "scientific_name")) %>%
  mutate(presence = if_else(is.na(presence), 0, presence))

# --- 6. Calculate mean length per haul and species -----------------------------

# Note: We use LngtClass and CANoAtLngt to calculate weighted mean
mean_lengths <- ca_sel %>%
  mutate(scientific_name = top_species_named$scientific_name[match(SpecCode, top_species_named$SpecCode)]) %>%
  group_by(Year, Survey, StNo, HaulNo, scientific_name) %>%
  summarise(
    mean_length = weighted.mean(LngtClass, CANoAtLngt, na.rm = TRUE),
    .groups = "drop"
  )

# --- 7. Add environmental covariates ------------------------------------------

covars <- hh_sel %>%
  dplyr::select(Year, Survey, StNo, HaulNo, ShootLat, ShootLong, Depth, BotTemp) %>%
  distinct()

# --- 8. Join all tables and finalize dataset -----------------------------------

sdm_data_multi <- all_pa %>%
  left_join(covars, by = c("Year", "Survey", "StNo", "HaulNo")) %>%
  left_join(mean_lengths, by = c("Year", "Survey", "StNo", "HaulNo", "scientific_name")) %>%
  filter(!is.na(Depth), !is.na(BotTemp)) %>%  # ¡NO filtramos por mean_length!
  mutate(
    Species = factor(scientific_name),
    # Imputamos la media general para ausencias (podrías usar otra estrategia)
    mean_length = ifelse(is.na(mean_length), mean(mean_length, na.rm = TRUE), mean_length),
    mean_length_scaled = scale(mean_length)
  )

summary(sdm_data_multi)

# --- 9. Export final dataset ---------------------------------------------------

saveRDS(sdm_data_multi, file = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/data/sdm_multispecies_clean.rds")

message("✅ Multispecies data prepared with scaled mean observed length per haul and species.")

# --- 10. Exploratory analyis ---------------------------------------------------



library(ggpubr)

p1 <- ggplot(sdm_data_multi, aes(x = as.factor(presence), y = BotTemp, fill = as.factor(presence))) +
  geom_boxplot() +
  labs(x = "Presence", y = "Bottom temperature", title = "Temperature by presence") +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "grey80", "1" = "tomato"))

p2 <- ggplot(sdm_data_multi, aes(x = as.factor(presence), y = Depth, fill = as.factor(presence))) +
  geom_boxplot() +
  labs(x = "Presence", y = "Depth", title = "Depth by presence") +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "grey80", "1" = "tomato"))

(p1 | p2)


# ggplot(sdm_data_multi, aes(x = mean_length, fill = Species)) +
#   geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
#   labs(title = "Distribution of mean length per species", x = "Mean length", y = "Frequency") +
#   facet_wrap(~ Species, scales = "free_y") +
#   theme_minimal() +
#   theme(legend.position = "none")

library(ggridges)

ggplot(sdm_data_multi, aes(x = mean_length, y = Species, fill = Species)) +
  geom_density_ridges(alpha = 0.8) +
  labs(title = "Distribution of mean length per species",
       x = "Mean length", y = "Species") +
  xlim(150, 300) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(sdm_data_multi, aes(x = mean_length, y = Species, fill = Species)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  scale_fill_paletteer_d("colorBlindness::ModifiedSpectralScheme11Steps") +
  labs(
    title = "Distribution of mean length per species",
    x = "Mean length",
    y = "Species"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


#NO ENTIENDO QUE SE QUIERE REPRESENTAR EN ESTOS PLOTS (?¿?)

# Filtrar solo presencias
sdm_data_pres <- sdm_data_multi %>%
  filter(presence == 1)

ggplot(sdm_data_multi_pres, aes(x = BotTemp, y = mean_length, color = as.factor(presence))) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
  facet_wrap(~ Species, scales = "free") +
  labs(x = "Bottom temperature", y = "Mean length", color = "Presence",
       title = "Trait vs environment colored by presence") +
  theme_minimal()

ggplot(sdm_data_multi, aes(x = BotTemp, y = mean_length, color = as.factor(presence))) +
  geom_point(alpha = 0.4) +
  facet_wrap(~ Species, scales = "free") +
  labs(x = "Bottom temperature", y = "Mean length", color = "Presence",
       title = "Trait vs environment colored by presence") +
  theme_minimal()
