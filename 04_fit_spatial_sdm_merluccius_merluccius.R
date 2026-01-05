# ==============================================================================
# Title: Fit Spatial SDMs for Merluccius merluccius With and Without Trait (INLA + SPDE)
# Author: M. Grazia Pennino
# Date:   2025-07-08
# Description:
#   - Load prepared SDM data (haul-level mean_length_cm and FishBase max length)
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
library(mapdata)
world <- map_data("world")
class(world)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires) 
library(INLAspacetime)
library(inlabru)
library(patchwork)


# --- 1. Load prepared SDM data ------------------------------------------------

sdm_data <- readRDS("C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/data/sdm_data_merluza.rds")

# --- 2. Scale covariates and define variables --------------------------------

df <- sdm_data %>%
  mutate(
    mean_len_s = as.numeric(scale(mean_length_cm)),        # observed mean haul length
    Depth_s    = as.numeric(scale(Depth)),                  # bottom depth
    BotTemp_s  = as.numeric(scale(BotTemp)),                # bottom temperature
    BotSal_s   = as.numeric(scale(BotSal)),                 # bottom salinity
    FB_len_s   = as.numeric(scale(FB_max_length_cm)),       # FishBase maximum length
    year_f     = as.factor(Year)  # Ensure year is factor
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

# Define fixed effect sets
fixed_no_trait  <- env_vars
fixed_with_trait <- c(env_vars, "mean_len_s")

# --- 4. Prepare spatial data (FOR A BARRIER MODEL, NOT USED)---------------------------------------------------

#We are going to consider a barrier model as we have several coast lines and 
#we have also islands in the area sampled 

land <- ne_countries(scale = "large", returnclass = "sf")
land_utm <- st_transform(land, 32630)

# bbox
bbox <- st_as_sfc(st_bbox(c(xmin=-13.88, xmax=1.35, ymin=40, ymax=55), crs = 4326))
bbox_utm <- st_transform(bbox, 32630)

# recortar tierra al bbox
land_crop <- st_intersection(land_utm, bbox_utm)

# unir todo en un solo polígono
land_union <- st_union(land_crop)

# océano = bbox - tierra
ocean_poly <- st_difference(bbox_utm, land_union)

# Visualizar
plot(st_geometry(ocean_poly), col = "lightblue")
plot(st_geometry(land_union), add = TRUE, col = "grey40")

df_sf <- st_as_sf(df, coords = c("ShootLong", "ShootLat"), crs = 4326)
df_sf <- st_transform(df_sf, crs = st_crs(ocean_poly))
plot(st_geometry(ocean_poly))
plot(st_geometry(land), add = TRUE, col = "grey")
plot(st_geometry(df_sf), add = TRUE, col = "blue", cex = 0.5)

# --- 5.1 Build triangulation mesh (FOR A BARRIER MODEL, NOT USED) ------------------------------------------------

loc_sf <- st_as_sf(
  df,
  coords = c("ShootLong", "ShootLat"),
  crs = 4326
)

loc_utm <- st_transform(loc_sf, crs = st_crs(ocean_poly))
loc_xy <- st_coordinates(loc_utm)


ocean_sp <- as(ocean_poly, "Spatial")
boundary <- fmesher::fm_as_segm(ocean_sp)

mesh <- fmesher::fm_mesh_2d_inla(
  loc = loc_xy,
  boundary = boundary,
  max.edge = c(70000, 250000),
  cutoff = 20000,
  offset = c(70000, 150000), 
  crs = 32630
)

mesh$n

plot(mesh)
plot(ocean_poly, add = TRUE, border = "blue", lwd = 2)
points(loc_xy, col = "red", pch = 16, cex = 0.4)


#Flujo de trabajo con un modelo de barrera (según el tutorial de Krainski)

#Identificar los triangulos dentro del oceáno 

water.tri <- fmesher::fm_contains(
  ocean_poly, 
  y = mesh, 
  type = "centroid"
)
water.tri.idx <- water.tri[[1]]

all.tri <- seq_len(nrow(mesh$graph$tv))
barrier.tri <- setdiff(all.tri, water.tri.idx)

#Construir la matriz de precisión
# We consider the range for the barrier as a fraction of the range over the domain. 
# We just use half of the average rectangle edges as the range in the domain and 10% of it in the barrier.

sigma <- 1
x_range <- range(mesh$loc[, 1])
y_range <- range(mesh$loc[, 2])

dx <- diff(x_range)
dy <- diff(y_range)

r <- mean(c(dx, dy))
r
range_domain  <- 0.5  * r
range_barrier <- 0.05 * r

# We now have to compute the Finite Element matrices needed for the model discretization, as detailed in Bakka et al. (2019).
#Hay que modificar el operador diferencial para que la correlación no cruce as barreras, y que lo haga
#de manera muy débil
#Se construyen con esta funcion dos operadores diferenciales diferentes uno para el oceano y otro para la tierra, generando dos matrices 
#FEM: de masa y de rigidez (laplacina), operando en: κ2−Δ

bfem <- mesh2fem.barrier(mesh, barrier.tri)

#En el segundo paso creamos la matriz de precision del campo latenta

Q <- inla.barrier.q(
  bfem,
  ranges = c(range_domain, range_barrier),
  sigma  = 1
)


#Model fitting 

bmodel <- barrierModel.define(
  mesh = mesh, 
  barrier.triangles = barrier.tri,
  prior.range = c(range_domain, 0.05),  # mejor que 1 fijo
  prior.sigma = c(1, 0.01),
  range.fraction = 0.1
)

model <- presence ~ Intercept(1) +
  Depth_s + BotTemp_s + mean_len_s +
  I(BotTemp_s * mean_len_s) +
  f(year_f, model = "iid") +
  field(
    cbind(ShootLong, ShootLat),
    model = bmodel
  )

result <- bru(
  model,
  data   = df,
  family = "binomial"
)

spatial_field <- result$summary.random$field

# --- 5. Build triangulation mesh 
# --- 5.2 Build triangulation mesh 

coordinates(df) <- ~ShootLong + ShootLat
proj4string(df) <- CRS("+proj=longlat +datum=WGS84")
coords <- coordinates(df)

mesh <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(0.5, 2),
  cutoff   = 0.1
)

plot(mesh); points(df, col = "red", pch = 16, cex = 0.5)

mesh$n

# --- 6. Define PC priors for SPDE and select best ------------------------------

spde_options <- list(
  loose = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(1, 0.01), prior.sigma = c(1, 0.01)),
  tight = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05))
)


fit_spatial <- function(spde_model, covariates) {
  idx <- inla.spde.make.index("spatial.field", spde_model$n.spde)
  A <- inla.spde.make.A(mesh, loc = loc_xy)
  df_cov <- df %>% 
    mutate(time_f = as.factor(year_f)) %>%
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

# Compare DIC for priors
for (nm in names(results)) {
  m <- results[[nm]]$model
  cat(sprintf("%-6s DIC = %8.2f, WAIC = %8.2f\n", nm, m$dic$dic, m$waic$waic))
}

best_prior <- names(results)[which.min(sapply(results, function(x) x$model$dic$dic))]
message("Selected SPDE prior: ", best_prior)
spatial_with_trait <- results[[best_prior]]$model
stk_with_trait     <- results[[best_prior]]$stack

##BEST PRIOR IS TIGHT 

# inla.spde2.pcmatern(mesh, alpha = 2,
#                     prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05)

# --- 7a. Fit spatial model WITH trait --------------------------------------

spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                    prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05))

s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

#Define A
A <- inla.spde.make.A(mesh, coords)

#Define stack
stack.1 <- inla.stack(
  data   = list(y = df$pres),
  A      = list(A, 1),
  effects = list(
    s.index,  
    df %>% transmute(
      intercept = 1,
      Depth_s, BotTemp_s, mean_len_s, year_f
    )
  ),
  tag = "fit"
)

##Formula 

f.1 <- y ~ -1 +intercept + Depth_s + BotTemp_s + mean_len_s +
  BotTemp_s:mean_len_s + 
  f(year_f, model = 'iid') +
  f(spatial.field, model = spde)

spatial_with_trait <- inla(
  f.1,
  data              = inla.stack.data(stack.1),
  family            = "binomial",
  control.predictor = list(
    A       = inla.stack.A(stack.1),
    compute = TRUE
  ),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose           = TRUE
)

spatial_with_trait$dic$dic
spatial_with_trait$waic$waic

# --- 7a. Fit spatial model WITHOUT trait --------------------------------------

f.2 <- y ~ -1 +intercept + Depth_s + BotTemp_s + 
  f(year_f, model = 'iid') +
  f(spatial.field, model = spde)

spatial_without_trait <- inla(
  f.2,
  data              = inla.stack.data(stack.1),
  family            = "binomial",
  control.predictor = list(
    A       = inla.stack.A(stack.1),
    compute = TRUE
  ),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose           = TRUE
)

spatial_without_trait$dic$dic
spatial_without_trait$waic$waic

# --- 7b. Fit non-spatial models -----------------------------------------------

f.3 <- y ~ -1 +intercept + Depth_s + BotTemp_s + 
  f(year_f, model = 'iid')

model_ns_bas <- inla(
  f.3,
  data              = inla.stack.data(stack.1),
  family            = "binomial",
  control.predictor = list(
    A       = inla.stack.A(stack.1),
    compute = TRUE
  ),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose           = TRUE
)

f.4 <- y ~ -1 +intercept + Depth_s + BotTemp_s + mean_len_s +
  BotTemp_s:mean_len_s + 
  f(year_f, model = 'iid')  

model_ns_trait <- inla(
  f.4,
  data              = inla.stack.data(stack.1),
  family            = "binomial",
  control.predictor = list(
    A       = inla.stack.A(stack.1),
    compute = TRUE
  ),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose           = TRUE
)

# --- 8. Compare model fits ----------------------------------------------------

comparison <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial_WithTrait", "NonSpatial_NoTrait", "NonSpatial_WithTrait"),
  DIC   = c(spatial_without_trait$dic$dic, spatial_with_trait$dic$dic,
            model_ns_bas$dic$dic,      model_ns_trait$dic$dic),
  WAIC  = c(spatial_without_trait$waic$waic, spatial_with_trait$waic$waic,
            model_ns_bas$waic$waic,      model_ns_trait$waic$waic)
)
print(comparison)

# --- 9. Calculate and compare ROC/AUC -----------------------------------------

idx <- inla.stack.index(stack.1, "fit")$data

pred_ns_nt <- model_ns_bas$summary.fitted.values[idx, "mean"]
pred_ns_tr <- model_ns_trait$summary.fitted.values[idx, "mean"]
pred_sp_nt <- spatial_without_trait$summary.fitted.values[idx, "mean"]
pred_sp_tr <- spatial_with_trait$summary.fitted.values[idx, "mean"]

preds <- df %>%
  mutate(
    pred_ns_nt = pred_ns_nt,
    pred_ns_tr = pred_ns_tr,
    pred_sp_nt = pred_sp_nt,
    pred_sp_tr = pred_sp_tr
  )

roc_vals <- tibble::tibble(
  Model = comparison$Model,
  AUC   = c(
    auc(roc(preds$pres, preds$pred_sp_nt)),
    auc(roc(preds$pres, preds$pred_sp_tr)),
    auc(roc(preds$pres, preds$pred_ns_nt)),
    auc(roc(preds$pres, preds$pred_ns_tr))
  )
)

print(roc_vals)

# --- Plot combined spatial fields: F ------------------------------------

sp_means_wt <- spatial_with_trait$summary.random$spatial.field$mean
sp_means_nt <- spatial_without_trait$summary.random$spatial.field$mean

projr    <- inla.mesh.projector(mesh, dims = c(200, 200))
field_wt <- inla.mesh.project(projr, sp_means_wt)
field_nt <- inla.mesh.project(projr, sp_means_nt)


df_nt_merluccius <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_nt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs

df_wt_merluccius <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_wt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs


df_nt_crop_merl <- df_nt_merluccius |>
  dplyr::filter(
    between(x, -13.88, 1.35),
    between(y, 39.36, 54.45)
  )

df_wt_crop_merl <- df_wt_merluccius |>
  dplyr::filter(
    between(x, -13.88, 1.35),
    between(y, 39.36, 54.45)
  )

p_nt_merluccius <- ggplot(df_nt_crop_merl, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(aes(z = value), size = 3, stroke = 0.15) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-13.88, 1.35), ylim = c(40, 55)) +
  theme_classic()

p_wt_merluccius <-ggplot(df_wt_crop_merl, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(aes(z = value), size = 3, stroke = 0.15) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-13.88, 1.35), ylim = c(40, 55)) +
  theme_classic()


windows();(p_nt_merluccius | p_wt_merluccius)


#No traits
p_nt_merluccius <- ggplot(df_nt_merluccius, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    name = "Mean"
  ) +
  labs(
    title = "(a) Without trait",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic()

#With trait
p_wt_merluccius <- ggplot(df_wt_merluccius, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    name = "Mean"
  ) +
  labs(
    title = "(b) With trait",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic()

windows();(p_nt_merluccius | p_wt_merluccius)

# 1. Collect your models in a named list
models_merluccius <- list(
  Spatial_NoTrait_merluccius      = spatial_no_trait,
  Spatial_WithTrait_merluccius    = spatial_with_trait,
  NonSpatial_NoTrait_merluccius   = model_ns_base,
  NonSpatial_WithTrait_merluccius = model_ns_trait
)

saveRDS(
  models_merluccius,
  file = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/output/models_SIMULATED.rds"
)

