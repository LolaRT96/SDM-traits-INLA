# ==============================================================================
# Title: Fit Spatial SDMs for Simulated Data With and Without Trait (INLA + SPDE)
# Author: M. Grazia Pennino (adapted)
# Date:   2025-07-14
# Description:
#   - Simulate SDM-like data with trait and spatial coordinates
#   - Scale covariates and check for collinearity
#   - Define two PC priors for the SPDE, select best by DIC
#   - Fit four models:
#       1) Spatial model without trait
#       2) Spatial model with trait
#       3) Non-spatial model without trait
#       4) Non-spatial model with trait
#   - Evaluate models using DIC, WAIC, and ROC/AUC
#   - Plot mesh, spatial field, and ROC curves
# ==============================================================================

# --- 0. Load libraries ---------------------------------------------------------

library(INLA)    # Bayesian spatial modeling
library(sp)      # SpatialPointsDataFrame
library(dplyr)   # Data wrangling
library(fields)  # image.plot()
library(viridis) # Color scales
library(ggplot2) # Plotting
library(pROC)    # ROC/AUC calculations
library(tibble)
library(patchwork)
library(metR)

dir.create("plots/simulated_single", recursive = TRUE, showWarnings = FALSE)

# --- 1. Simulate SDM-like data ------------------------------------------------

set.seed(123)

n <- 1000

df <- tibble::tibble(
  presence = rbinom(n, 1, 0.5),
  ShootLong = runif(n, -10, 10),
  ShootLat  = runif(n, 35, 45),
  Depth     = rnorm(n, 200, 50),
  BotTemp   = rnorm(n, 12, 2),
  BotSal    = rnorm(n, 35, 1),
  mean_length_cm     = rnorm(n, 25, 5),
  FB_max_length_cm   = rnorm(n, 50, 10),
  Year      = sample(2015:2022, n, replace = TRUE)
)

# --- 2. Scale covariates and define variables --------------------------------

df <- df %>%
  mutate(
    mean_len_s = as.numeric(scale(mean_length_cm)),
    Depth_s    = as.numeric(scale(Depth)),
    BotTemp_s  = as.numeric(scale(BotTemp)),
    BotSal_s   = as.numeric(scale(BotSal)),
    FB_len_s   = as.numeric(scale(FB_max_length_cm)),
    year_f     = as.factor(Year)
  )

# --- 3. Check collinearity and select environmental covariates ---------------

cov_env <- df %>% select(Depth_s, BotTemp_s, BotSal_s)

cor_env <- cor(cov_env, use = "complete.obs")

print(cor_env)

high_corr <- which(abs(cor_env) > 0.7 & abs(cor_env) < 1, arr.ind = TRUE)

if (nrow(high_corr) > 0) {
  drop_var <- names(which.max(colMeans(abs(cor_env))))
  message("Dropping environmental covariate due to high collinearity: ", drop_var)
  env_vars <- setdiff(names(cov_env), drop_var)
} else {
  env_vars <- names(cov_env)
}
message("Using environmental covariates: ", paste(env_vars, collapse = ", "))

fixed_no_trait  <- env_vars #Solo las variables ambientales
fixed_with_trait <- c(env_vars, "mean_len_s") #variables ambientales + el trait 


# --- 4. Prepare spatial data ---------------------------------------------------



coordinates(df) <- ~ShootLong + ShootLat
proj4string(df) <- CRS("+proj=longlat +datum=WGS84")
coords <- coordinates(df)



# --- 5. Build triangulation mesh ------------------------------------------------

mesh <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(0.5, 2),
  cutoff   = 0.1
)

plot(mesh); points(df, col = "red", pch = 16, cex = 0.5)

# --- 6. Define PC priors for SPDE and select best ------------------------------

spde_options <- list(
  loose = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(1, 0.01), prior.sigma = c(1, 0.01)),
  
  tight = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05))
)

fit_spatial <- function(spde_model, covariates) {
  
  idx <- inla.spde.make.index("spatial.field", spde_model$n.spde)
  A <- inla.spde.make.A(mesh, loc = coords)
  
  df_cov <- df@data %>% 
    mutate(year_f = as.factor(Year)) %>%
    select(all_of(covariates), year_f)
  
  stk <- inla.stack(
    data = list(presence = df$presence),
    A = list(A, 1),
    effects = list(
      spatial.field = idx,
      data = df_cov
    ),
    tag = "est"
  )
  
  formula <- as.formula(paste(
    "presence ~", paste(c(covariates, "f(year_f, model = 'iid')"), collapse = " + "),
    "+ f(spatial.field, model = spde_model)"
  ))
  
  res <- inla(
    formula,
    family = "binomial",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE)
  )
  
  list(model = res, stack = stk)
}

results <- lapply(spde_options, fit_spatial, covariates = fixed_with_trait)


for (nm in names(results)) {
  m <- results[[nm]]$model
  cat(sprintf("%-6s DIC = %8.2f, WAIC = %8.2f\n", nm, m$dic$dic, m$waic$waic))
}

best_prior <- names(results)[which.min(sapply(results, function(x) x$model$dic$dic))]

message("Selected SPDE prior: ", best_prior)

spatial_with_trait <- results[[best_prior]]$model #extrae la fomula del mejor modelo 

stk_with_trait     <- results[[best_prior]]$stack #extrae la fomula del mejor stack 

##Priors seleccionados: LOOSE
##prior.range = c(1, 0.01), prior.sigma = c(1, 0.01)

# --- 7a. Fit spatial model WITHOUT trait --------------------------------------

#La funcion predefinida anteriormente es: fit_espacial(spde_model, covariates)

f_spatial_nt <- fit_spatial(spde_options[[best_prior]], fixed_no_trait)

spatial_no_trait <- f_spatial_nt$model

stk_no_trait     <- f_spatial_nt$stack

# --- 7b. Fit non-spatial models -----------------------------------------------

fmla_ns_base  <- as.formula(paste("presence ~", paste(fixed_no_trait, collapse = " + ")))

fmla_ns_trait <- as.formula(paste("presence ~", paste(fixed_with_trait, collapse = " + ")))

model_ns_base  <- inla(fmla_ns_base,  family = "binomial", data = df@data,
                       control.compute = list(dic = TRUE, waic = TRUE))

model_ns_trait <- inla(fmla_ns_trait, family = "binomial", data = df@data,
                       control.compute = list(dic = TRUE, waic = TRUE))

# --- 8. Compare model fits ----------------------------------------------------

comparison <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial_WithTrait", "NonSpatial_NoTrait", "NonSpatial_WithTrait"),
  DIC   = c(spatial_no_trait$dic$dic, spatial_with_trait$dic$dic,
            model_ns_base$dic$dic,      model_ns_trait$dic$dic),
  WAIC  = c(spatial_no_trait$waic$waic, spatial_with_trait$waic$waic,
            model_ns_base$waic$waic,      model_ns_trait$waic$waic)
)
print(comparison)

# --- 9. Calculate and compare ROC/AUC -----------------------------------------

preds <- df@data %>%
  mutate(
    pred_sp_nt    = spatial_no_trait$summary.fitted.values[ inla.stack.index(stk_no_trait,     "est")$data, "mean"],
    pred_sp_trait = spatial_with_trait$summary.fitted.values[ inla.stack.index(stk_with_trait, "est")$data, "mean"],
    pred_ns_nt    = model_ns_base$summary.fitted.values$mean,
    pred_ns_trait = model_ns_trait$summary.fitted.values$mean
  )
roc_vals <- tibble::tibble(
  Model = comparison$Model,
  AUC   = c(
    auc(roc(preds$presence, preds$pred_sp_nt)),
    auc(roc(preds$presence, preds$pred_sp_trait)),
    auc(roc(preds$presence, preds$pred_ns_nt)),
    auc(roc(preds$presence, preds$pred_ns_trait))
  )
)

print(roc_vals)

# --- 10. Plot combined spatial field with and without trait --------------------

sp_means_wt <- spatial_with_trait$summary.random$spatial.field$mean
sp_means_nt <- spatial_no_trait$summary.random$spatial.field$mean

projr       <- inla.mesh.projector(mesh, dims = c(200, 200))
field_wt    <- inla.mesh.project(projr, sp_means_wt)
field_nt    <- inla.mesh.project(projr, sp_means_nt)

# Save combined figure with two panels
png("plots/simulated_single/spatial_fields_combined.png", width = 1200, height = 600)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 5))  # Adjust margins

# Panel (a): Without trait
image.plot(projr$x, projr$y, field_nt, col = viridis(100),
           xlab = "Longitude", ylab = "Latitude", asp = 1,
           main = "(a) ")
contour(projr$x, projr$y, field_nt, add = TRUE, col = "black", lwd = 0.4)

# Panel (b): With trait
image.plot(projr$x, projr$y, field_wt, col = viridis(100),
           xlab = "Longitude", ylab = "Latitude", asp = 1,
           main = "(b) ")
contour(projr$x, projr$y, field_wt, add = TRUE, col = "black", lwd = 0.4)


###MEJOR USA GGPLOT 

df_nt <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_nt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs

df_wt <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_wt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs

df_nt_crop <- df_nt |>
  dplyr::filter(
    between(x, -15, 15),
    between(y, 30, 50)
  )

df_wt_crop <- df_wt |>
  dplyr::filter(
    between(x, -15, 15),
    between(y, 30, 50)
  )

p_nt <-ggplot(df_nt_crop, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(aes(z = value), size = 3, stroke = 0.15) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  coord_equal(expand = FALSE) +
  theme_classic()

p_wt <-ggplot(df_wt_crop, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(aes(z = value), size = 3, stroke = 0.15) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  coord_equal(expand = FALSE) +
  theme_classic()

# p_wt <- ggplot(df_wt, aes(x, y, fill = value)) +
#   geom_raster() +
#   geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
#   geom_text_contour(
#     aes(z = value),
#     stroke = 0.15,
#     size = 3,
#     skip = 0      # <- etiqueta TODAS las líneas
#   ) +
#   scale_fill_distiller(
#     palette = "RdBu",
#     direction = -1,
#     name = "Mean"
#   ) +
# coord_cartesian(
#   xlim = c(-14.64, 14.64),
#   ylim = c(30.37, 49.66)
# ) +
#   labs(
#     title = "(b) With trait",
#     x = "Longitude",
#     y = "Latitude"
#   ) +
#   theme_classic()

# #No traits
# p_nt <- ggplot(df_nt, aes(x, y, fill = value)) +
#   geom_raster() +
#   geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
#   geom_text_contour(
#     aes(z = value),
#     stroke = 0.15,
#     size = 3,
#     skip = 0      # <- etiqueta TODAS las líneas
#   ) +
#   scale_fill_distiller(
#     palette = "RdBu",
#     direction = -1,
#     name = "Mean"
#   ) +
#   coord_cartesian(
#     xlim = c(-15, 15),
#     ylim = c(30, 50)
#   ) +
#   labs(
#     title = "(a) Without trait",
#     x = "Longitude",
#     y = "Latitude"
#   ) +
#   theme_classic()

windows();(p_nt | p_wt)

#OBSERVANDO LOS PATRONES DE VARIACIÓN:

# El campo espacial positivo indica zonas donde la probabilidad predicha es mayor 
# de lo que explican las covariables; negativo indica zonas donde es menor. 
# Son efectos residuales suavizados del proceso espacial.
# 



# # ROC comparison plot
# # png("plots/simulated_single/roc_comparison_simulated.png", width = 800, height = 600)
# plot(roc(preds$presence, preds$pred_sp_trait), col = "blue", lwd = 2, main = "ROC Comparison")
# lines(roc(preds$presence, preds$pred_sp_nt),    col = "green",  lwd = 2)
# lines(roc(preds$presence, preds$pred_ns_trait), col = "red",    lwd = 2)
# legend("bottomright",
#        legend = c("Spatial + Trait", "Spatial only", "Non-spatial + Trait"),
#        col    = c("blue", "green", "red"),
#        lwd    = 2)
# 

# 11. Print model summaries

models_SIMULATED <- list(
  Spatial_NoTrait_SIMULATED      = spatial_no_trait,
  Spatial_WithTrait_SIMULATED    = spatial_with_trait,
  NonSpatial_NoTrait_SIMULATED   = model_ns_base,
  NonSpatial_WithTrait_SIMULATED = model_ns_trait
)

saveRDS(
  models_SIMULATED,
  file = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/output/models_SIMULATED.rds"
)

for (nm in names(models)) {
  cat("========================================\n")
  cat("Model:", nm, "\n")
  cat("========================================\n\n")
  print(summary(models[[nm]]))
  cat("\n\n")
}

message("✅ All simulated models fitted, evaluated, and plotted successfully.")
