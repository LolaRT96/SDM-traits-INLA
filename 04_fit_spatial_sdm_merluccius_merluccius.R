# ==============================================================================
# Title: Fit Spatial SDMs for Merluccius merluccius With and Without Trait (INLA + SPDE)
# Author: M. Grazia Pennino (MODIFIED BY LOLA RIESGO)
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
#   - Plot spatial field
#   - Model evaluation of the best model (spatial with trait)
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
glimpse(sdm_data)
sdm_data$Survey <- factor(sdm_data$Survey)
levels(sdm_data$Survey)

# sdm_dataNORTH <- sdm_data %>% filter(Survey %in% c("SP-NORTH")) %>%
#   droplevels()

# --- 2. Scale covariates and define variables --------------------------------

df <- sdm_data %>%
  mutate(
    mean_len_s = as.numeric(scale(mean_length_mm)),        # observed mean haul length
    Depth_s    = as.numeric(scale(Depth)),                  # bottom depth
    BotTemp_s  = as.numeric(scale(BotTemp)),                # bottom temperature
    BotSal_s   = as.numeric(scale(BotSal)),                 # bottom salinity
    FB_len_s   = as.numeric(scale(FB_max_length_cm)),       # FishBase maximum length
    year_f     = as.factor(Year)  # Ensure year is factor
  )

summary(sdm_data)

summary(df) 

ggplot(df, aes(x= year_f, y = mean_length_mm))+
  geom_boxplot()+
  theme_minimal()

ggplot(df, aes(x= year_f, y = BotTemp))+
  geom_point()+
  theme_minimal()

#Si la temperatura es -9 es que NO hay registro 

ggplot(df, aes(x = ShootLong, y = ShootLat, color = BotTemp)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(name = "Bottom temperature (°C)", na.value = "grey80") +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Spatial distribution of bottom temperature"
  )

#Limpiar
df_clean <- df %>%
  mutate(
    BotTemp = na_if(BotTemp, -9),
    BotSal  = na_if(BotSal,  -9)
  )

summary(df_clean$BotTemp)
summary(df_clean$mean_length_mm)
df_clean$mean_length_cm <- df_clean$mean_length_mm /10
summary(df_clean$mean_length_cm)

vars_to_scale <- c("mean_len_s", "Depth_s", "BotTemp_s", "BotSal_s")

df_1 <- df_clean %>%
  filter(if_all(all_of(vars_to_scale), ~ . >= -3.3 & . <= 3.3))

summary(df_1)

# --- 3. Check collinearity and select environmental covariates ---------------

cov_env <- df %>% select(Depth_s, BotTemp_s, BotSal_s)
cor_env <- cor(cov_env, use = "complete.obs")
print(cor_env)

# Define fixed effect sets
fixed_no_trait  <- c("BotTemp_s", "Depth_s")
fixed_with_trait <- c("BotTemp_s", "Depth_s", "mean_len_s")

# --- 5.2 Build triangulation mesh 

loc <- cbind(df_1$ShootLong, df_1$ShootLat)
loc <- as.data.frame(loc)
colnames(loc) <- c("ShootLong", "ShootLat")
coordinates(loc) <- ~ShootLong + ShootLat
proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")
coords <- coordinates(loc)

mesh <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(0.5, 2),
  cutoff   = 0.1
)

plot(mesh);points(df_1, col = "red", pch = 16, cex = 0.5)

mesh$n

# --- 6. Define PC priors for SPDE and select best ------------------------------

df_1 <- as.data.frame(df_1)

spde_options <- list(
  loose = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(1, 0.01), prior.sigma = c(1, 0.01)),
  tight = inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05))
)


fit_spatial <- function(spde_model, covariates, df) {
  
  idx <- inla.spde.make.index("spatial.field", spde_model$n.spde)
  
  coords <- as.matrix(df[, c("ShootLong", "ShootLat")])
  
  A <- inla.spde.make.A(mesh, loc = coords)
  
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
    "presence ~", 
    paste(c(covariates, "f(year_f, model = 'iid')"), collapse = " + "),
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


results <- lapply(
  spde_options, 
  fit_spatial, 
  covariates = fixed_with_trait,
  df = df_1
)


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

coordinates(df_1) <- ~ShootLong + ShootLat
proj4string(df_1) <- CRS("+proj=longlat +datum=WGS84")
coords <- coordinates(df_1)

mesh <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(0.5, 2),
  cutoff   = 0.1
)
plot(mesh);points(df_1, col = "red", pch = 16, cex = 0.5)

mesh$n

df_1 <- as.data.frame(df_1)

spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                    prior.range = c(0.5, 0.05), prior.sigma = c(0.5, 0.05))

s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

#Define A
A <- inla.spde.make.A(mesh, coords)

#Define stack
stack.1 <- inla.stack(
  data   = list(y = df_1$pres),
  A      = list(A, 1),
  effects = list(
    s.index,  
    df_1 %>% transmute(
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

model_ns_bas$dic$dic
model_ns_bas$waic$waic

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
   geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "Spatial effect") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-13.88, 1.35), ylim = c(40, 55)) +
  theme_classic()+
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

p_wt_merluccius <-ggplot(df_wt_crop_merl, aes(x, y, fill = value)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "Spatial effect") +
  # geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "black") + 
  coord_fixed(xlim = c(-13.88, 1.35), ylim = c(40, 55)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))


(p_nt_merluccius | p_wt_merluccius)

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



# EXPLORACIÓN DEL MEJOR MODELO  -------------------------------------------


#Relación de la temperatura y la longitud media 

spatial_with_trait$summary.fixed

temp_seq <- seq(min(df_1$BotTemp_s), max(df_1$BotTemp_s), length.out = 50)
length_seq <- seq(
  min(df_1$mean_len_s, na.rm = TRUE),
  max(df_1$mean_len_s, na.rm = TRUE),
  length.out = 50
)

grid <- expand.grid(BotTemp_s = temp_seq, mean_len_s = length_seq)
grid$intercept <- 1
grid$Depth_s <- mean(df_1$Depth_s)   # fijar otras variables en su media

X <- model.matrix(~ -1 + intercept +Depth_s + BotTemp_s + mean_len_s + BotTemp_s:mean_len_s, data = grid)
beta <- spatial_with_trait$summary.fixed$mean
grid$eta <- as.vector(X %*% beta)
grid$prob <- 1 / (1 + exp(-grid$eta))  # probabilidad binomial






ggplot(grid, aes(x = BotTemp_s, y = mean_len_s, fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "Bottom Temperature (°C)", y = "Mean body size (cm)", fill = "Probability of presence") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

##Valores desescalados

df$temp_s <- (df$BotTemp_s - mean(df$BotTemp_s)) / sd(df$BotTemp_s)

df$length_cm_s <- (df$mean_length_cm - mean(df$mean_length_cm)) / sd(df$mean_length_cm)

mean_temp <- mean(df$BotTemp_s)
sd_temp   <- sd(df$BotTemp_s)

mean_length <- mean(df$mean_length_cm, na.rm = TRUE)
sd_length   <- sd(df$mean_length_cm, na.rm = TRUE)

temp_orig_seq <- seq(min(df$BotTemp_s), max(df$BotTemp_s), length.out = 50)
length_orig_seq <- seq(min(df$mean_length_cm, na.rm = TRUE), max(df$mean_length_cm, na.rm = TRUE), length.out = 50)

grid_orig <- expand.grid(
  temp = temp_orig_seq,
  length_cm = length_orig_seq)


grid_orig$temp_s <- (grid_orig$temp - mean_temp) / sd_temp

grid_orig$length_cm_s <- (grid_orig$length_cm - mean_length) / sd_length

grid_orig$intercept <- 1

grid_orig$bathy_s <- mean(df$Depth_s)  

X <- model.matrix(~ -1 + intercept + bathy_s + temp_s + length_cm_s + temp_s:length_cm_s, data = grid_orig)
beta <- spatial_with_trait$summary.fixed$mean
grid_orig$eta <- as.vector(X %*% beta)

grid_orig$prob <- 1 / (1 + exp(-grid_orig$eta))

ggplot(grid_orig, aes(x = temp, y = length_cm, fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "Bottom Temperature (°C)", y = "Mean body size (cm)", fill = "Probability of presence") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )






#Distribución de las marginales 

marginals_fixed <- spatial_with_trait$marginals.fixed

marginals_df <- lapply(names(marginals_fixed), function(param) {
  df <- as.data.frame(marginals_fixed[[param]])
  colnames(df) <- c("x", "y")
  df$parameter <- param
  df
})

marginals_combined <- bind_rows(marginals_df)

# Paso 2: graficar con facets y escalas libres
ggplot(marginals_combined, aes(x = x, y = y)) +
  geom_line(color = "steelblue", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred") +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Distribuciones Marginales de los Efectos Fijos",
    x = "Valor del coeficiente",
    y = "Densidad"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
