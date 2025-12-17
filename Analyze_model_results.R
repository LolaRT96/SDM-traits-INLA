# ==============================================================================
# Title: Analyze and Visualize Model Results for Traits in Space Pipeline
# Author: M.Grazia Pennino
# Date: 2025-07-08
# Description: Loads fitted spatial and multispecies SDM models,
#              extracts model fit metrics,
#              visualizes fixed effects and species-specific slopes,
#              plots spatial predictions (example),
#              performs trait sensitivity analysis (example),
#              compares model predictive performance with ROC curves.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(pROC)
library(INLA)
library(raster)

# --- 1. Load fitted models -----------------------------------------------------

model_spatial <- readRDS("output/model_sdm_merluza_spatial.rds")

model_multispecies <- readRDS("output/model_sdm_multispecies_hierarchical.rds")

# --- 2. Extract and compare model fit metrics ----------------------------------

metrics <- tibble(
  Model = c("Single-species Spatial SDM", "Multispecies Hierarchical SDM"),
  DIC = c(model_spatial$dic$dic, model_multispecies$dic$dic),
  WAIC = c(model_spatial$waic$waic, model_multispecies$waic$waic)
)

print(metrics)

# --- 3. Extract and plot fixed effects -----------------------------------------

fixef_spatial <- model_spatial$summary.fixed %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Model = "Single-species")

fixef_multi <- model_multispecies$summary.fixed %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Model = "Multispecies")

fixef_all <- bind_rows(fixef_spatial, fixef_multi)

# Plot fixed effects (excluding intercept)
ggplot(fixef_all %>% filter(Parameter != "(Intercept)"), 
       aes(x = Parameter, y = mean, color = Model)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.2, position = position_dodge(width = 0.6)) +
  labs(title = "Fixed Effects Estimates",
       y = "Estimate",
       x = NULL) +
  theme_minimal()

ggsave("plots/fixed_effects_comparison.png", width = 8, height = 5)

# --- 4. Extract and plot species-specific slopes (random effects) ---------------

# Assuming multispecies model has a random effect named "slope_species"
# Adjust name if your model structure differs
species_slopes <- model_multispecies$summary.random$slope_species %>%
  as_tibble() %>%
  mutate(SpeciesID = row_number())

ggplot(species_slopes, aes(x = reorder(SpeciesID, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.3) +
  labs(title = "Species-specific Temperature Slopes",
       x = "Species ID",
       y = "Slope Estimate") +
  coord_flip() +
  theme_minimal()

ggsave("plots/species_slopes.png", width = 7, height = 6)

# --- 5. Plot spatial predictions (example with raster) --------------------------

# Example: create a raster with spatial prediction mean from the single-species model
# (In practice, replace with your spatial prediction extraction)
if ("spatial.field" %in% names(model_spatial$summary.random)) {
  spatial_mean <- model_spatial$summary.random$spatial.field$mean
  mesh <- model_spatial$.args$data$mesh  # Assumed stored mesh, adapt if necessary
  proj <- inla.mesh.projector(mesh, dims = c(100, 100))
  field_projected <- inla.mesh.project(proj, spatial_mean)
  
  raster_pred <- raster::raster(t(field_projected), 
                                xmn = min(proj$x), xmx = max(proj$x), 
                                ymn = min(proj$y), ymx = max(proj$y),
                                crs = NA)
  
  png("plots/spatial_prediction_single_species.png", width = 800, height = 600)
  plot(raster_pred, main = "Spatial Prediction - Single Species Model")
  dev.off()
} else {
  message("No spatial.field component found in model_spatial.")
}

# --- 6. Trait sensitivity analysis (example) -----------------------------------

# Simulate trait effect curves using fixed intercept and varying slope coefficients
temp_seq <- seq(5, 25, length.out = 100)
beta_0 <- model_spatial$summary.fixed["(Intercept)", "mean"]
beta_depth <- model_spatial$summary.fixed["depth", "mean"]

trait_values <- c(0.5, 1.0, 1.5, 2.0)
curves <- map_df(trait_values, function(trait_val) {
  beta_temp <- 0.3 * trait_val
  eta <- beta_0 + beta_temp * temp_seq + beta_depth * mean(model_spatial$.args$data$depth)
  tibble(temp = temp_seq, prob = plogis(eta), trait = as.factor(trait_val))
})

ggplot(curves, aes(x = temp, y = prob, color = trait)) +
  geom_line(size = 1.2) +
  labs(title = "Trait Sensitivity Analysis",
       x = "Temperature (°C)",
       y = "Predicted Probability of Presence",
       color = "Trait Value") +
  theme_minimal()

ggsave("plots/trait_sensitivity_analysis.png", width = 8, height = 5)

# --- 7. ROC curve and AUC comparison --------------------------------------------

# Load prediction data saved after model fitting
preds <- readRDS("output/predictions_merluza.rds")  # Dataframe with presence, pred_spatial, pred_multispecies

roc_spatial <- roc(preds$presence, preds$pred_spatial)
roc_multi <- roc(preds$presence, preds$pred_multispecies)

cat("AUC Single-species Spatial Model: ", auc(roc_spatial), "\n")
cat("AUC Multispecies Hierarchical Model: ", auc(roc_multi), "\n")

png("plots/roc_curve_comparison.png", width = 800, height = 600)
plot(roc_spatial, col = "blue", main = "ROC Curve Comparison")
plot(roc_multi, col = "red", add = TRUE)
legend("bottomright", legend = c("Single-species", "Multispecies"),
       col = c("blue", "red"), lwd = 2)
dev.off()

# --- End of script ---
