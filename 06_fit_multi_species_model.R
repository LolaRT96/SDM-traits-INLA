
# ==============================================================================
# Title: Fit Spatial SDMs for Multispecies Assemblage With and Without Trait (INLA + SPDE)
# Author: M. Grazia Pennino & M.D. Riesgo

# Date:   2025-07-14

# Description:
#   - Load cleaned multispecies SDM dataset from ICES DATRAS (presence–absence, mean length)
#   - Filter extreme values, temperature outliers, and rare species (< 30 presences)
#   - Build SPDE mesh and define priors
#   - Fit two models:
#       1) Spatial model without trait-modulated response
#       2) Spatial model with trait-modulated temperature slope
#   - Evaluate models using DIC, WAIC, and ROC/AUC
#   - Visualize spatial fields and slope–trait relationships
# ==============================================================================


# --- 1. Load libraries --------------------------------------------------------

library(dplyr)
library(ggrepel)
library(tidyr)
library(ggplot2)
library(INLA)
library(viridis)
library(pROC)
library(sp)
library(mapdata)
world <- map_data("world")
class(world)
library(showtext)
library(sysfonts)

font_add("Helvetica", 
         regular = "~/Downloads/helvetica-255/Helvetica.ttf")
showtext_auto()

#dir.create("plots/empirical_multi", recursive = TRUE, showWarnings = FALSE)

# --- 2. Load cleaned SDM dataset ----------------------------------------------

sdm_data_multi <-  readRDS("~/SDMs_Traits/data/sdm_multispecies_clean.rds")

# --- 3. Apply data quality filters --------------------------------------------
# To ensure model robustness and remove artifacts or extreme values, we apply the following filters:
# (1) Remove extremely large length values that may be rare or miscoded.
# (2) Filter out impossible bottom temperature values (e.g., -9°C).
# (3) Remove species with fewer than 30 presence observations to ensure reliable estimation.

# 3.1. Remove implausible or extreme individual sizes (trait) 
#     - Values above 600 mm are rare and may distort trait scaling.

sdm_data_multi <- sdm_data_multi %>%
  filter(mean_length < 600)

# 3.2. Remove implausible bottom temperatures
#     - Temperatures below -2°C or above 20°C are likely errors or outliers.
sdm_data_multi <- sdm_data_multi %>%
  filter(BotTemp > -2 & BotTemp < 20)

# 3.3. Remove species with too few presence records (< 30)
#     - Prevents instability in model fitting for very rare species.
species_with_enough_data <- sdm_data_multi %>%
  group_by(Species) %>%
  summarise(n_presences = sum(presence)) %>%
  filter(n_presences >= 30) %>%
  pull(Species)

sdm_data_multi <- sdm_data_multi %>%
  filter(Species %in% species_with_enough_data)

# 3.4. (Optional) Rescale mean_length by species if modeling within-species responses
#     - For now we keep global scaling as we model across species.
#     - Uncomment below if within-species standardization is preferred:
# sdm_data_multi <- sdm_data_multi %>%
#   group_by(Species) %>%
#   mutate(mean_length_scaled = scale(mean_length)) %>%
#   ungroup()

# Summary of cleaned data
summary(sdm_data_multi)

## 3.5. We do not want every species, select the important

levels(sdm_data_multi$Species)
levels(sdm_data_multi$Survey)

sdm_data_multi_df <- as.data.frame(sdm_data_multi)

sdm_data_multi_filter <- sdm_data_multi_df %>%
 filter(Species %in% c("Lepidorhombus whiffiagonis", 
            "Trachurus trachurus", "Scomber scombrus",  
            "Lepidorhombus boscii", "Merluccius merluccius", 
            "Micromesistius poutassou")) %>%
  droplevels()

levels(sdm_data_multi_filter$Species)

# plot histogram of sizes (cm)

ggplot(data = sdm_data_multi, aes(x=mean_length)) +
  geom_histogram( binwidth=15, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_classic() +
  facet_wrap(~Species)+
  theme(
    plot.title = element_text(size=15)
  )

ggplot(sdm_data_multi, aes(x=Depth, y = mean_length))+
  geom_point(aes(color = BotTemp)) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_classic() +
  facet_wrap(~Species) 


ggplot(sdm_data_multi, aes(x = ShootLong, y = ShootLat, color = BotTemp)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(name = "Bottom temperature (°C)", na.value = "grey80") +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Spatial distribution of bottom temperature"
  )

ggplot(sdm_data_multi, aes(x = ShootLong, y = ShootLat, color = mean_length)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(name = "Mean Length cm", na.value = "grey80") +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Spatial distribution of bottom temperature"
  )


sdm_data_multi$mean_length_cm <- sdm_data_multi$mean_length / 10

ggplot(sdm_data_multi, aes(x = ShootLong, y = ShootLat, color = mean_length_cm)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(name = "Mean Length cm", na.value = "grey80") +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Spatial distribution of bottom temperature"
  )


p_surveys_multi <- ggplot() +
  geom_point(data = sdm_data_multi_filter, 
             aes(x = ShootLong, y = ShootLat),
             alpha = 0.4, size = 2,  color = "steelblue") +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "gray95") + 
  coord_fixed(xlim = c(-10, -1), ylim = c(41, 46)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


# --- 3. Create SPDE mesh and model --------------------------------------------
# Convert coordinates to spatial format
coordinates(sdm_data_multi_filter) <- ~ShootLong + ShootLat

# Build mesh using point locations
mesh <- inla.mesh.2d(
  loc = coordinates(sdm_data_multi_filter),
  max.edge = c(0.5, 2),  # Finer mesh near observations, coarser in outer
  cutoff   = 0.1         # Minimum distance between points
)

plot(mesh, main = "SPDE Mesh for Multispecies SDM")
points(sdm_data_multi_filter@coords, col = "red", pch = 16, cex = 0.4)


# Define SPDE model with PC priors
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(1, 0.01),  # P(range < 1 deg) = 0.01
  prior.sigma = c(1, 0.01)   # P(sigma > 1) = 0.01
)

# Create index for the spatial field
s_index <- inla.spde.make.index("spatial", spde$n.spde)

# Projection matrix linking observations to mesh nodes
A_matrix <- inla.spde.make.A(mesh = mesh, loc = coordinates(sdm_data_multi_filter))

# --- 4. Scale covariates ---------------------------------------
# Extract attribute table from spatial object
sdm_df <- sdm_data_multi_filter@data

# Scale covariates and rename trait variable for clarity
sdm_df <- sdm_df %>%
  mutate(
    Depth_s        = scale(Depth),
    BotTemp_s      = scale(BotTemp),
    mean_length_s  = scale(mean_length)  # (optional: already scaled in previous step)
  )

glimpse(sdm_df)


# --- 5. Define INLA stack and formulas ----------------------------------------

# Convert species to factor and create slope ID for trait-modulated response
sdm_df$Species <- factor(sdm_df$Species)
sdm_df$temp_slope_id <- as.integer(sdm_df$Species)

# Build list of covariates and random effects
effects_list <- list(
  spatial = s_index,
  data.frame(
    intercept       = 1,
    Depth_s         = sdm_df$Depth_s,
    BotTemp_s       = sdm_df$BotTemp_s,
    mean_length_s   = sdm_df$mean_length_s,
    Species         = sdm_df$Species,
    Year            = sdm_df$Year,
    temp_slope_id   = sdm_df$temp_slope_id
  )
)

# Build the INLA stack
stack <- inla.stack(
  data = list(presence = sdm_df$presence),
  A = list(A_matrix, 1),
  effects = effects_list,
  tag = "est"
)


formula_trait <- presence ~ BotTemp_s * mean_length_s + 
  Depth_s +
  f(Species, model = "iid")  +
  f(temp_slope_id, #indice de los niveles del efecto aleatorio (identifica cada especie)
    BotTemp_s, #covariable con la que se va a relacionar, la pendiente de la temperatura
    model = "iid") + #cada especie tiene un efecto indedependiente de pendiente (slopes diferentes)
  f(Year, model = "iid") +
  f(spatial, model = spde)

# Define formula WITHOUT trait
formula_no_trait <- presence ~ BotTemp_s + Depth_s +
  f(Species, model = "iid") +
  f(Year, model = "iid") +
  f(spatial, model = spde)

# Define formula trait WITHOSPATIAL
formula <- presence ~ BotTemp_s * mean_length_s + 
  Depth_s +
  f(Species, model = "iid")  +
  f(temp_slope_id, #indice de los niveles del efecto aleatorio (identifica cada especie)
    BotTemp_s, #covariable con la que se va a relacionar, la pendiente de la temperatura
    model = "iid") + #cada especie tiene un efecto indedependiente de pendiente (slopes diferentes)
  f(Year, model = "iid") 


# --- 6. fit models---------------------------------------


model <- inla(
  formula_trait,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE)
)

model_nt <- inla(
  formula_no_trait,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE)
)

model_nsp <- inla(
  formula,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE)
)

model_nsp$dic$dic
model_nsp$waic$waic

# --- 7. Compare model fits ----------------------------------------------------

comparison <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial", "Model_nsp"),
  DIC   = c(model_nt$dic$dic,  model$dic$dic, model_nsp$dic$dic),
  WAIC  = c(model_nt$waic$waic,  model$waic$waic, model_nsp$waic$waic)
)

comparison <- comparison %>%
  mutate(
    DIC  = formatC(DIC, format = "f", digits = 2),
    WAIC = formatC(WAIC, format = "f", digits = 2)
  )

print(comparison)

# --- 8. Calculate and compare ROC/AUC -----------------------------------------

# Identify indices corresponding to the "est" tag in the stack
idx_est <- inla.stack.index(stack, tag = "est")$data

# Add predicted values for both models to the original dataframe
preds <- sdm_df %>%
  mutate(
    pred_nt    = model_nt$summary.fitted.values$mean[idx_est],
    pred_trait = model_lola2$summary.fitted.values$mean[idx_est],
    pred_trait_np = model_nsp$summary.fitted.values$mean[idx_est]
  )

# Compute ROC and AUC
roc_trait <- roc(preds$presence, preds$pred_trait)
roc_nt    <- roc(preds$presence, preds$pred_nt)
roc_trait_np    <- roc(preds$presence, preds$pred_trait_np)

roc_vals <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial_WithTrait", "NonSpatial_WithTrait"),
  AUC   = c(auc(roc_nt), auc(roc_trait), auc(roc_trait_np))
)

print(roc_vals)


# --- 11. Plot spatial fields --------------------------------------------------
# Create projection grid over mesh
projr <- inla.mesh.projector(mesh, dims = c(200, 200))
field <- inla.mesh.project(projr, model$summary.random$spatial$mean)
field_nt <- inla.mesh.project(projr, model_nt$summary.random$spatial$mean)

df_nt_multiReal <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_nt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs


df_trait_multiReal <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field)) %>%
  filter(!is.na(value))      

#No traits
spatialfield_nt <- ggplot(df_nt_multiReal, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    name = "Spatial effect"
  ) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-12, -1), ylim = c(40, 46)) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic()+
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))


spatialfield <- ggplot(df_trait_multiReal, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    name = "Spatial effect"
  ) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-12, -1), ylim = c(40, 46)) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

(spatialfield_nt | spatialfield)

combination <- (spatialfield_nt | spatialfield)


# --- 13.Plot: Trait-Modulated Thermal Response across Species with Labels-------------------------------------------------

# Extract posterior means and 95% CI of species-specific temperature slopes

slopes <- model$summary.random$temp_slope_id %>%
  as_tibble() %>%
  rename(
    SpeciesID  = ID,
    mean_slope = mean,
    lower      = `0.025quant`,
    upper      = `0.975quant`
  )

# Compute mean observed trait (length) per species
trait_table <- sdm_df %>%
  group_by(Species) %>%
  summarise(
    mean_length = mean(mean_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SpeciesID = as.integer(Species))

# Merge trait values with estimated slopes
slopes_with_trait <- left_join(slopes, trait_table, by = "SpeciesID")

glimpse(slopes_with_trait)

# Create the plot object
plot_slope <- ggplot(slopes_with_trait, aes(x = mean_length, y = mean_slope)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  labs(
    x = "Mean observed length (cm)",
    y = "Slope of temperature effect",
    title = "Trait-Modulated Thermal Response across Species"
  ) +
  theme_minimal()

plot_slope


trait_table <- sdm_df %>%
  group_by(Species) %>%
  summarise(
    mean_length = mean(mean_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SpeciesID = as.integer(Species))


slopes_with_trait <- left_join(slopes_species, trait_table, by = "SpeciesID")


slopes_with_trait <- slopes_with_trait %>%
  arrange(desc(mean_length)) %>%
  mutate(Species = factor(Species, levels = Species))


ggplot(slopes_with_trait, 
       aes(x = mean_slope, 
           y = Species,
           color = (mean_length/10))) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y",
                width = 0,
                linewidth = 1) +
  geom_point(size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4, color = "red") +
  scale_color_paletteer_c("viridis::magma", name = "Mean length (cm)") +
  labs(
    x = "Temperature effect (slope)",
    y = "Species"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14))


library(RColorBrewer)

colors_div <- brewer.pal(n = 8, name = "Spectral")

thermal_slope <- ggplot(slopes_with_trait, 
       aes(x = mean_slope, 
           y = Species,
           color = (mean_length/10))) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y",
                width = 0,
                linewidth = 1) +
  geom_point(size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4, color = "red") +
  #scale_color_paletteer_c("viridis::magma", name = "Mean length (cm)") +
  scale_color_distiller(
    palette = "Spectral",
    name = "Mean length (cm)"
  ) +
  labs(
    x = "Temperature effect (slope)",
    y = "Species"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14))

thermal_slope






