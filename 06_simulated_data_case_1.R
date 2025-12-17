# ==============================================================================
# Title: Simulated Single-Species Spatial SDM with Trait-Modulated Response
# Author: M.Grazia Pennino
# Date: 2025-07-08
# Description: Simulates a single species with a trait-modulated temperature response,
#              fits spatial and non-spatial SDMs using INLA+SPDE,
#              generates diagnostic plots and ROC validation,
#              performs trait sensitivity analysis.
# ==============================================================================

# --- Load required packages ----------------------------------------------------

library(INLA)
library(sp)
library(tidyverse)
library(fields)     # for image.plot()
library(viridis)    # color palette
library(pROC)       # ROC curves and AUC

# --- Create output directory ---------------------------------------------------

if(!dir.exists("plots")) dir.create("plots")

# --- 1. Simulate environmental data and species occurrence ----------------------

set.seed(42)
n_sites <- 200

sim_data <- tibble(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100),
  temp = rnorm(n_sites, mean = 15, sd = 3),
  depth = runif(n_sites, 10, 200)
)

# Trait modulating temperature effect (e.g. thermal tolerance)
trait_temp_tolerance <- 1.5

# Model coefficients
beta_0 <- -2.0
beta_temp <- 0.3 * trait_temp_tolerance
beta_depth <- -0.01

# Calculate presence probability and simulate binary response
sim_data <- sim_data %>%
  mutate(
    logit_p = beta_0 + beta_temp * temp + beta_depth * depth,
    p_presence = plogis(logit_p),
    presence = rbinom(n_sites, 1, p_presence)
  )

# --- 2. Build spatial mesh for SPDE --------------------------------------------

coordinates(sim_data) <- ~x + y

mesh <- inla.mesh.2d(
  loc = coordinates(sim_data),
  max.edge = c(10, 20),
  cutoff = 1
)

png("plots/mesh_plot.png", width = 800, height = 600)
plot(mesh)
points(sim_data, col = "red", pch = 16, cex = 0.6)
dev.off()

# --- 3. Define SPDE model and create INLA stack ---------------------------------

spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
s_index <- inla.spde.make.index("spatial.field", spde$n.spde)

A_matrix <- inla.spde.make.A(mesh, loc = coordinates(sim_data))

stack <- inla.stack(
  data = list(presence = sim_data$presence),
  A = list(A_matrix, 1),
  effects = list(s_index, data.frame(temp = sim_data$temp, depth = sim_data$depth)),
  tag = "est"
)

# --- 4. Fit spatial SDM with INLA ----------------------------------------------

formula <- presence ~ temp + depth + f(spatial.field, model = spde)

model_spatial <- inla(
  formula,
  family = "binomial",
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

print(summary(model_spatial))

# --- 5. Plot spatial field (posterior mean) ------------------------------------

spatial_mean <- model_spatial$summary.random$spatial.field$mean
proj <- inla.mesh.projector(mesh, dims = c(200, 200))
field_projected <- inla.mesh.project(proj, spatial_mean)

png("plots/spatial_field.png", width = 800, height = 600)
image.plot(
  proj$x, proj$y, field_projected,
  col = viridis(100),
  xlab = "X", ylab = "Y",
  asp = 1,
  main = "Posterior Mean Spatial Field"
)
contour(proj$x, proj$y, field_projected, add = TRUE, col = "black", lwd = 0.4)
dev.off()

# --- 6. Trait-modulated temperature response curve with credible intervals ------

beta_fixed <- model_spatial$summary.fixed

temp_seq <- seq(min(sim_data$temp), max(sim_data$temp), length.out = 100)
depth_mean <- mean(sim_data$depth)

predict_response <- function(beta0, beta_temp, beta_depth, temps, depth) {
  eta <- beta0 + beta_temp * temps + beta_depth * depth
  plogis(eta)
}

plot_df <- tibble(
  temp = temp_seq,
  p_mean = predict_response(beta_fixed["(Intercept)", "mean"], beta_fixed["temp", "mean"], beta_fixed["depth", "mean"], temp_seq, depth_mean),
  p_lower = predict_response(beta_fixed["(Intercept)", "0.025quant"], beta_fixed["temp", "0.025quant"], beta_fixed["depth", "0.025quant"], temp_seq, depth_mean),
  p_upper = predict_response(beta_fixed["(Intercept)", "0.975quant"], beta_fixed["temp", "0.975quant"], beta_fixed["depth", "0.975quant"], temp_seq, depth_mean)
)

ggplot(plot_df) +
  geom_line(aes(temp, p_mean), color = "blue", linewidth = 1.2) +
  geom_ribbon(aes(x = temp, ymin = p_lower, ymax = p_upper), alpha = 0.3, fill = "white") +
  labs(x = "Temperature (°C)", y = "Probability of Presence",
       title = "Trait-modulated Temperature Response with 95% Credible Interval") +
  theme_minimal()

ggsave("plots/temp_response_inla_spatial.png", width = 7, height = 5)

# --- 7. Fit non-spatial model for comparison -----------------------------------

model_nospatial <- inla(
  presence ~ temp + depth,
  family = "binomial",
  data = as.data.frame(sim_data),
  control.compute = list(dic = TRUE, waic = TRUE)
)

print(summary(model_nospatial))

# --- 8. Compare model metrics: DIC and WAIC ------------------------------------

cat("Spatial model DIC:", model_spatial$dic$dic, "\n")
cat("Non-spatial model DIC:", model_nospatial$dic$dic, "\n")
cat("Spatial model WAIC:", model_spatial$waic$waic, "\n")
cat("Non-spatial model WAIC:", model_nospatial$waic$waic, "\n")

# --- 9. ROC validation for both models -----------------------------------------

df_preds <- as.data.frame(sim_data) %>%
  mutate(
    pred_spatial = model_spatial$summary.fitted.values[inla.stack.index(stack, "est")$data, "mean"],
    pred_nospatial = model_nospatial$summary.fitted.values$mean
  )

roc_spatial <- roc(df_preds$presence, df_preds$pred_spatial)
roc_nospatial <- roc(df_preds$presence, df_preds$pred_nospatial)

cat("AUC Spatial Model:", auc(roc_spatial), "\n")
cat("AUC Non-Spatial Model:", auc(roc_nospatial), "\n")

png("plots/roc_comparison.png", width = 800, height = 600)
plot(roc_spatial, col = "blue", main = "ROC Curve Comparison")
plot(roc_nospatial, col = "red", add = TRUE)
legend("bottomright", legend = c("Spatial", "Non-Spatial"), col = c("blue", "red"), lwd = 2)
dev.off()

# --- 10. Trait sensitivity analysis ---------------------------------------------

trait_values <- c(0.5, 1.0, 1.5, 2.0)
curves <- map_df(trait_values, function(tr_val) {
  beta_temp_val <- 0.3 * tr_val
  p <- predict_response(beta_0, beta_temp_val, beta_depth, temp_seq, depth_mean)
  tibble(temp = temp_seq, prob = p, trait = factor(tr_val))
})

ggplot(curves, aes(temp, prob, color = trait)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Temperature (°C)",
    y = "Probability of Presence",
    color = "Trait Value",
    title = "Trait Sensitivity Analysis"
  ) +
  theme_minimal()

ggsave("plots/sensitivity_trait_curves.png", width = 7, height = 5)
