
# ==============================================================================
# Title: Figure S1 – Map of haul locations (EVHOE and SP-NORTH)
# Author: M. Grazia Pennino
# Date:   2025-07-08
#
# Description:
#   This script generates a supplementary figure (Figure S1) showing the spatial
#   distribution of haul locations from the EVHOE (Bay of Biscay) and SP-NORTH
#   (Northern Spain) bottom-trawl surveys during 2010–2022. Coordinates are taken
#   from ICES DATRAS haul header (HH) data. The map provides spatial coverage
#   context for the empirical case studies (European hake and Iberian assemblages).
#
# Output:
#   - S1_map_hauls.png (PNG figure saved in the working directory)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- 1. Load basemap (coastline) and data ----------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
extent_sf <- st_as_sfc(st_bbox(c(xmin = -12, xmax =  3,
                                 ymin =  35, ymax = 47),
                               crs = st_crs(world)))

hh_all <- readRDS("data/hh_all.rds")    # all hauls, all species
ca_all <- readRDS("data/ca_all.rds")    # all catch records

# --- 2. Extract unique hauls with coordinates ---------------------------------
hauls_pts <- hh_all %>%
  distinct(Year, Survey, StNo, HaulNo, ShootLong, ShootLat) %>%
  filter(!is.na(ShootLong), !is.na(ShootLat)) %>%
  mutate(Survey = factor(Survey, levels = c("EVHOE","SP-NORTH")))

hauls_sf <- st_as_sf(hauls_pts,
                     coords = c("ShootLong","ShootLat"),
                     crs = 4326,
                     remove = FALSE)

# --- 3. Build the map ---------------------------------------------------------
p_map <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.3) +
  geom_sf(data = hauls_sf, aes(color = Survey), alpha = 0.4, size = 0.7) +
  coord_sf(xlim = c(-12, 3), ylim = c(35, 47), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(title = " ",
       x = "Longitude", y = "Latitude", color = "Survey") +
  theme_minimal(base_size = 12)

# --- 4. Save the figure -------------------------------------------------------
ggsave("S1_map_hauls.png", p_map, width = 7.5, height = 6, dpi = 300)
