# ==============================================================================
# Title: Fit Spatial SDMs for Simulated Data With and Without Trait (INLA + SPDE)
# Author: M. Grazia Pennino & M.D. Riesgo

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
library(tibble)
library(patchwork)
library(metR)
library(withr)
library(fmesher)
library(showtext)
library(sysfonts)
library(scico)

font_add("Helvetica", 
         regular = "C:/Users/mdolores.riesgo/Downloads/helvetica-255/Helvetica.ttf")
showtext_auto()


# --- 1. Simulate SDM-like data ------------------------------------------------


set.seed(123)

n_time <- 7 #years 2015:2022 

#Create the grid (latitude and longitude) from 0 to 100 
#From the Trust-Your-Model script 

campo <- function(c1, c2, c3, c4) {
  xy <- expand.grid(
    seq(c1, c2, length.out = 100),
    seq(c3, c4, length.out = 100)
  )
  cbind(xy[, 1], xy[, 2])
}

c1 <- 0; c2 <- 100; c3 <- 0; c4 <- 100 

loc_xy <- campo(c1, c2, c3, c4) 

#Mesh

with_seed(123, {mesh2d_sim <- 
  fm_mesh_2d_inla(loc.domain = loc_xy, 
                  max.edge = c(8, 15))  
})

windows();plot(mesh2d_sim, main = "Malla", asp = 1, lwd = 0.5)

mesh2d_sim$n


# Depth  ------------------------------------------------------------------

grid_list <- vector("list", n_time)

for (t in seq_len(n_time)) {
  
  loc_df <- data.frame(
    x <- loc_xy[,1],
    y <- loc_xy[,2]
  )
  
  depth_range <- c(50, 200)
  
  loc_df$bathy <- depth_range[1] +
    ((loc_df$y - c3) / (c4 - c3)) * diff(depth_range)
  
  grid_list[[t]] <- data_frame (
    x =loc_df$x,
    y= loc_df$y,
    bathy = loc_df$bathy,
    time = t )
}

grid_bathy_df <- bind_rows(grid_list)
head(grid_bathy_df)

p_bathy <- ggplot(grid_bathy_df, aes(x =x, y = y, fill = bathy)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "G", direction =-1)+
  labs(
    x = "x",
    y = "y",
    fill = "Bathymetry (m)"  # ← aquí cambias el nombre de la leyenda
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )

p_bathy

# Bottom temperature ---------------------------------------------------------

# Anomalías temporales suaves
anomalies <- cumsum(rnorm(n_time, mean = 0, sd = 0.2))

# Función que convierte profundidad en temperatura

temp_from_depth <- function(depth) {
  17 - (depth - 50) * (12 / 150)
}

grid_temp_list <- vector("list", n_time)

for (t in seq_len(n_time)) {
  
  loc_df <- grid_bathy_df %>% filter(time == t)
  
  anomaly_t <- anomalies[t]
  
  # Temperatura final
  loc_df$temp <- temp_from_depth(loc_df$bathy) + 
    anomaly_t + 
    rnorm(nrow(loc_df), mean = 0, sd = 0.3)  # ruido espacial pequeño
  
  grid_temp_list[[t]] <- loc_df
}

grid_temp_df <- bind_rows(grid_temp_list)

p_temp <- ggplot(grid_temp_df, aes(x = x, y = y, fill = temp)) +
  geom_tile() +
  # facet_wrap(~time, ncol = 4) +
  coord_equal() +
  scale_fill_viridis_c(option = "H") +
  labs(x = "x", y = "y", fill = "Temp (°C)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )

p_temp 

temp_time_df <- grid_temp_df %>%
  group_by(time) %>%
  summarise(temp_mean = mean(temp), 
            temp_sd = sd(temp), .groups = "drop") %>%
  mutate(temp_upper = temp_mean + temp_sd,
         temp_lower = temp_mean - temp_sd)

p_temp_time <- ggplot(temp_time_df, aes(x = time, y = temp_mean)) +
  geom_line(color = "steelblue", size = 1) +
  labs(x = "Time (unit)", y = "Mean Temp (°C)")+
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

p_temp_time

# Bottom salinity ---------------------------------------------------------

sal_from_temp <- function(temp) {
  
  35 - (temp - 5) * (2 / 12)
}

grid_sal_list <- vector("list", n_time)

for (t in seq_len(n_time)) {
  
  loc_df <- grid_temp_df %>% filter(time == t)
  
  # 
  loc_df$sal <- sal_from_temp(loc_df$temp) +
    rnorm(nrow(loc_df), mean = 0, sd = 0.1)  
  
  grid_sal_list[[t]] <- loc_df
}

grid_sal_df <- bind_rows(grid_sal_list)

ggplot(grid_sal_df, aes(x = x, y = y, fill = sal)) +
  geom_tile() +
  facet_wrap(~time, ncol = 3) +
  coord_equal() +
  scale_fill_viridis_c(option = "C") +
  labs(x = "x", y = "y", fill = "Salinity (PSU)")


# Spatio temporal structure (Gaussian Random Field)  ----------------------------------------------

with_seed(123, {
  
  #Precision
  prec <- 1 / 4
  rho <- 25
  sigma <- sqrt(7)
  phi <- 0.8
  
  sigma_eps <- sigma * sqrt(1 - phi^2) #Noise deviation: each time unit we add noise of a specific amplitude. 
  
  u0_nodes <- fmesher::fm_matern_sample(mesh2d_sim, n = 1, rho = rho, sigma = sigma) #Guassian Random Field
  
  u0_nodes <- u0_nodes - mean(u0_nodes) #center at zero to eliminate unwanted deviations from the field. 
  
  n_time <- 7 
  A_grid <- fm_basis(mesh2d_sim, loc_xy) #matrix that associates each node of the mesh with each point of the grid
  
  latent_list <-  vector("list", n_time)
  
  u_prev <- u0_nodes #node starting point
  
  for(t in seq(n_time)) {
    
    eps_t <-  fm_matern_sample(mesh2d_sim, n=1, rho = rho, sigma = sigma) #noise
    eps_t <- eps_t - mean(eps_t) #mean to zero
    
    u_t <- phi * u_prev + eps_t
    u_prev <- u_t #result at the nodes for time t
    
    latent_t <- drop(A_grid %*% u_t) #we interpolate those values u_t at each grid point 
    
    latent_list[[t]] <- data_frame(
      x =loc_xy[,1],
      y=loc_xy[,2],
      latent = latent_t,
      time = t
    )
    
  }
  
  latent_time_df <- bind_rows(latent_list)
  
} )

ggplot(latent_time_df, aes(x=x, y=y, fill=latent)) +
  geom_tile()+
  facet_wrap(~time,  nrow = 3) +
  coord_equal()+
  scale_fill_viridis_c() +
  labs(
    x = "x",
    y = "y",
    fill = "Spatial structure"  
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "bottom",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        
         strip.text        = element_text(face = "bold", size = 12)  
  )


# Area occupied -----------------------------------------------
#Body size is going to regulate the environmental answer to depth to temperature 
#Small body: warm waters 
#Large body: depth water/cold waters 

df_grid <- grid_temp_df %>%
  select(x, y, time, temp) %>%           # temp SOLO aquí
  inner_join(
    latent_time_df %>% select(x, y, time, latent),
    by = c("x","y","time")
  ) %>%
  inner_join(
    grid_bathy_df %>% select(x, y, time, bathy),
    by = c("x","y","time")
  ) %>%
  inner_join(
    grid_sal_df %>% select(x, y, time, sal),
    by = c("x","y","time")
  )

glimpse(df_grid)



n_time <- 7

# pesos ambientales
beta_latent <- 4.5
beta_temp   <- 2.2
beta_sal    <- 1.0
beta_bathy  <- 0.05

# thermal niche 
beta_size        <- 0.6
beta_temp_size   <- -0.08
beta_bathy_size <-  0.03    

# ocupación
frac_ocup_media <- 0.45
sd_frac_ocup    <- 0.03

df_list <- vector("list", n_time)

with_seed(123, {
  
  for (t in seq_len(n_time)) {
    
    grid_t <- df_grid %>% filter(time == t)
    
    frac_ocup_t <- min(
      max(rnorm(1, frac_ocup_media, sd_frac_ocup), 0.15),
      0.40
    )
    
    # gradiente ontogénico espacial (latente)
    grid_t <- grid_t %>%
      mutate(
        size_latent = scale(bathy)[,1] - scale(temp)[,1]
      )
    
    # score ecológico
    grid_t <- grid_t %>%
      mutate(
        score = beta_latent * latent +
          beta_temp   * temp +
          beta_sal    * sal +
          beta_bathy  * bathy +
          beta_size        * size_latent +
          beta_temp_size   * size_latent * temp +
          beta_bathy_size  * size_latent * bathy
      )
    
    # umbral para fijar ocupación
    umbral <- quantile(grid_t$score, probs = 1 - frac_ocup_t)
    
    grid_t <- grid_t %>%
      mutate(pres = ifelse(score > umbral, 1, 0))
    
    df_list[[t]] <- grid_t
  }
})

df_ocup <- bind_rows(df_list)


# Sampling (n=1000) -------------------------------------------------------

n_pts <- 1000
target_ratio_bounds <- c(0.4, 0.5)

ratio_vec <- runif(n_time,
                   target_ratio_bounds[1],
                   target_ratio_bounds[2])

df_muestreo_list <- vector("list", n_time)

with_seed(456, {
  
  for (t in seq_len(n_time)) {
    
    grid_t <- df_ocup %>% filter(time == t)
    
    pres_t <- grid_t %>% filter(pres == 1)
    abs_t  <- grid_t %>% filter(pres == 0)
    
    ratio_t <- ratio_vec[t]
    
    n_abs  <- floor(n_pts / (1 + ratio_t))
    n_pres <- n_pts - n_abs
    
    n_abs  <- min(n_abs,  nrow(abs_t))
    n_pres <- min(n_pres, nrow(pres_t))
    
    muestra_t <- bind_rows(
      sample_n(abs_t,  n_abs),
      sample_n(pres_t, n_pres)
    )
    
    df_muestreo_list[[t]] <- muestra_t
  }
})

df_muestra <- bind_rows(df_muestreo_list)

#Body size vary with temp and depth

df_muestra <- df_muestra %>%
  mutate(
    mu_length =
      35 +
      0.05 * bathy -      
      0.6  * temp,        
    length_cm = rnorm(n(), mean = mu_length, sd = 4)
  )

#Body size VS depth 

ggplot(df_muestra, aes(bathy, length_cm)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(x = "Depth (m)", y = "Body length (cm)")

#Body size VS temp 

ggplot(df_muestra, aes(temp, length_cm)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(x = "Bottom temperature (°C)", y = "Body length (cm)")

#area ocupada
ggplot(df_ocup, aes(x = x, y = y, fill = factor(pres))) +
  geom_tile(color = NA) +
  coord_equal() +
  scale_fill_manual(values = c("0" = "grey90", "1" = "deepskyblue"),
                    labels = c("Absence", "Presence")) +
  facet_wrap(~ time, ncol = 3) +
  labs(
    x = "x", 
    y = "y") +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "bottom",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                       
         strip.text        = element_text(face = "bold", size = 12)  
  )


##TWO DATAFRAMES READY

#df_muestra = df_fit #for fitting 
#df_ocup = df_predict #for predict 

# --- 2. Scale covariates and define variables --------------------------------

df <- df_muestra  %>%
  mutate(
    length_cm_s = as.numeric(scale(length_cm)),
    bathy_s    = as.numeric(scale(bathy)),
    temp_s  = as.numeric(scale(temp)),
    sal_s   = as.numeric(scale(sal)),
    time_f     = as.factor(time)
  )

# --- 3. Check collinearity and select environmental covariates ---------------

cov_env <- df %>% select(bathy_s, temp_s, sal_s)

cor_env <- cor(cov_env, use = "complete.obs")

print(cor_env)

high_corr <- which(abs(cor_env) > 0.7 & abs(cor_env) < 1, arr.ind = TRUE)
#temperatura y salinidad estan correlacionadas 
#nos quedamos con la salinidad 

env_vars <- df %>% select(bathy_s, temp_s)
fixed_no_trait  <- c("bathy_s", "temp_s")
fixed_with_trait <- c("bathy_s", "temp_s", "length_cm_s")


# --- 4. Prepare spatial data ---------------------------------------------------

loc <- data.frame(x = df$x, y = df$y)
coordinates(loc) <- ~ x + y
coords <- coordinates(loc)

# --- 5. Build triangulation mesh ------------------------------------------------

mesh <- fm_mesh_2d_inla(loc.domain = loc, 
                max.edge = c(7, 13),
                cutoff = 0.1)  
mesh$n

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
  A <- inla.spde.make.A(mesh, loc)
  
  df_cov <- df %>% 
    mutate(time_f = as.factor(time)) %>%
    select(all_of(covariates), time_f)
  
  stk <- inla.stack(
    data = list(presence = df$pres),
    A = list(A, 1),
    effects = list(
      spatial.field = idx,
      data = df_cov
    ),
    tag = "est"
  )
  
  formula <- as.formula(paste(
    "presence ~", paste(c(covariates, "f(time_f, model = 'iid')"), collapse = " + "),
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

# --- 7a. SPATIAL WITHOUT TRAITS --------------------------------------

#To have a better control of predictors its better to write all functions 
#Define spde
spde <- inla.spde2.pcmatern(mesh, prior.range = c(1, 0.01), prior.sigma = c(1, 0.01))
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

#Define A
A <- inla.spde.make.A(mesh, loc)


#Define stack
stack.1 <- inla.stack(
  data   = list(y = df$pres),
  A      = list(A, 1),
  effects = list(
    s.index,  
    df %>% transmute(
      intercept = 1,
      temp_s, time, bathy_s, length_cm_s
    )
  ),
  tag = "fit"
)

#Define formula: individuals are goint to response different regarding their length
#effect of temperature depends on the length
#size deendend on temp 

f.1 <- y ~ -1 +intercept + bathy_s + temp_s + 
       f(spatial.field, model = spde) +
       f(time, model = 'iid') 

#Model spatial NO traits 

spatial_no_trait <- inla(
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

# --- 7b. SPATIAL WITH TRAITS --------------------------------------

f.2 <- y ~ -1 + intercept + bathy_s + temp_s + length_cm_s +
  temp_s:length_cm_s + f(spatial.field, model = spde) +
  f(time, model = 'iid') 

#Model spatial WITH traits

spatial_with_trait <- inla(
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

# --- 7c. NO SPATIAL WITH TRAITS --------------------------------------

f.3 <- y ~ -1 + intercept + bathy_s + temp_s + length_cm_s +
  temp_s:length_cm_s + 
  f(time, model = 'iid') 

model_ns_trait <- inla(
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


# --- 7c. NO SPATIAL NO TRAITS --------------------------------------

f.4 <- y ~ -1 + intercept + bathy_s + temp_s +
  f(time, model = 'iid')

model_ns_base <- inla(
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
  DIC   = c(spatial_no_trait$dic$dic, spatial_with_trait$dic$dic,
            model_ns_base$dic$dic,      model_ns_trait$dic$dic),
  WAIC  = c(spatial_no_trait$waic$waic, spatial_with_trait$waic$waic,
            model_ns_base$waic$waic,      model_ns_trait$waic$waic), 
  LPML_CPO = c(
    sum(log(spatial_with_trait$cpo$cpo), na.rm = TRUE),
    sum(log(spatial_no_trait$cpo$cpo), na.rm = TRUE),
    sum(log(model_ns_trait$cpo$cpo), na.rm = TRUE),
    sum(log(model_ns_base$cpo$cpo), na.rm = TRUE)
  )
)

comparison <- comparison %>%
  mutate(
    DIC  = formatC(DIC, format = "f", digits = 2),
    WAIC = formatC(WAIC, format = "f", digits = 2),
    LPML_CPO = formatC(LPML_CPO, format = "f", digits = 2),
  )

print(comparison)

# Relative contribution of traits -----------------------------------------

# Baseline model

dic_base  <- model_ns_base$dic$dic
waic_base <- model_ns_base$waic$waic

# Other models

dic_ns_trait <- model_ns_trait$dic$dic
dic_s_base   <- spatial_no_trait$dic$dic
dic_s_trait  <- spatial_with_trait$dic$dic

waic_ns_trait <- model_ns_trait$waic$waic
waic_s_base   <- spatial_no_trait$waic$waic
waic_s_trait  <- spatial_with_trait$waic$waic

# ΔDIC (improvement vs baseline)

delta_dic_ns_trait <- dic_base - dic_ns_trait
delta_dic_s_base   <- dic_base - dic_s_base
delta_dic_s_trait  <- dic_base - dic_s_trait
delta_dic_spatial_trait  <- dic_s_base - dic_s_trait

# ΔWAIC

delta_waic_ns_trait <- waic_base - waic_ns_trait
delta_waic_s_base   <- waic_base - waic_s_base
delta_waic_s_trait  <- waic_base - waic_s_trait
delta_waic_spatial_trait  <- waic_s_base - waic_s_trait

model_comparison <- data.frame(
  Model = c("Baseline vs No Spatial with trait", 
            "Baseline vs Spatial without trait", 
            "Baseline vs Spatial with trait",
            "Spatial without trait vs Spatial with trait"),
  
  dDIC  = c(delta_dic_ns_trait,
            delta_dic_s_base,
            delta_dic_s_trait,
            delta_dic_spatial_trait),
  
  dWAIC = c(delta_waic_ns_trait,
            delta_waic_s_base,
            delta_waic_s_trait,
            delta_waic_spatial_trait)
)

model_comparison

# --- 9. Calculate and compare ROC/AUC -----------------------------------------

idx <- inla.stack.index(stack.1, "fit")$data

pred_ns_nt <- model_ns_base$summary.fitted.values[idx, "mean"]
pred_ns_tr <- model_ns_trait$summary.fitted.values[idx, "mean"]
pred_sp_nt <- spatial_no_trait$summary.fitted.values[idx, "mean"]
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

# Plot calibration checks  -------------------------------------------------

#Index of real observation samples
idx <- inla.stack.index(stack.1, "fit")$data

# Funtion for calibration
calibration_data <- function(model, model_name, y) {
  
  p <- model$summary.fitted.values[idx, "mean"]
  
  tibble(
    obs = y,
    pred = p
  ) %>%
    mutate(bin = ntile(pred, 10)) %>%
    group_by(bin) %>%
    summarise(
      pred_mean = mean(pred),
      obs_mean  = mean(obs),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(model = model_name)
}

cal_df <- bind_rows(
  calibration_data(
    spatial_with_trait,
    "Spatial, with trait",
    df$pres
  ),
  calibration_data(
    spatial_no_trait,
    "Spatial, no trait",
    df$pres
  ),
  calibration_data(
    model_ns_trait,
    "No spatial, with trait",
    df$pres
  ),
  calibration_data(
    model_ns_base,
    "No spatial, no trait",
    df$pres
  )
)

# Plot
cal.plots <-  ggplot(cal_df,
       aes(x = pred_mean,
           y = obs_mean)) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.8, color = "gray"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ model) +
  coord_equal() +
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  ) +
  labs(
    x = "Mean predicted probability",
    y = "Observed prevalence")

cal.plots

ggsave(
  filename = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/plots/calplots.jpg",
  plot = cal.plots,
  width = 1484,    # ancho en píxeles
  height = 758,   # alto en píxeles
  units = "px",
  dpi = 72        # dpi estándar para píxeles (72 dpi)
)

# --- 10. Plot combined spatial field with and without trait --------------------

sp_means_wt <- spatial_with_trait$summary.random$spatial.field$mean
sp_means_nt <- spatial_no_trait$summary.random$spatial.field$mean

projr <- fmesher::fm_evaluator(
  mesh = mesh,
  dims = c(200, 200)
)

field_wt <- fmesher::fm_evaluate(
  projector = projr,
  field     = sp_means_wt
)

field_nt    <- fmesher::fm_evaluate(
  projector = projr,
  field     = sp_means_nt
)

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

p_nt <-ggplot(df_nt, aes(x, y, fill = value)) +
  geom_raster() +
  #scale_fill_distiller(palette = "RdBu", direction = -1,  name = "Spatial effect") +
  scale_fill_scico(palette = "vik", ,  name = "Spatial effect") +
  xlim(0,100) + ylim(0,100)+
  coord_equal(expand = FALSE) +
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )

p_wt <-ggplot(df_wt, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_scico(palette = "vik",  name = "Spatial effect") +
  coord_equal(expand = FALSE) +
  xlim(0,100) + ylim(0,100)+
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )

combination <- (p_nt | p_wt)

combination

#Differences

df_both <- df_nt %>%
inner_join(
df_wt,
by = c("x", "y"),
suffix = c("_nt", "_wt")
)

df_both <- df_both %>% 
  rename( value_nt = value_nt,
          value_wt = value_wt)

head(df_both)

df_both$diff <- df_both$value_nt - df_both$value_wt
head(df_both)
df_both$sim <- -abs(df_both$diff)
head(df_both)

df_sim <- df_both %>% 
  mutate(
  diff = abs(value_nt - value_wt),
  similarity = ifelse(diff > 0.1, "Not similar", "Similar"))
head(df_sim)


p_sim <-ggplot(df_sim, aes(x, y, fill = sim)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  xlim(0,100) + ylim(0,100) +
    scale_fill_viridis(name = "Similarity")+
  theme_classic() +
  theme( 
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )

p_sim


combination2 <- (p_wt | p_sim)

combination2

ggsave(
  filename = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/plots/simulationSingleCase_change.jpg",
  plot = combination2,
  width = 972,    # ancho en píxeles
  height = 380,   # alto en píxeles
  units = "px",
  dpi = 72        # dpi estándar para píxeles (72 dpi)
)


# 11. Save the models

models_SIMULATED <- list(
  Spatial_NoTrait_SIMULATED      = spatial_no_trait,
  Spatial_WithTrait_SIMULATED    = spatial_with_trait,
  NonSpatial_NoTrait_SIMULATED   = model_ns_base,
  NonSpatial_WithTrait_SIMULATED = model_ns_trait
)

saveRDS(
  models_SIMULATED,
  file = "~/SDMs_Traits/output/models_SIMULATED.rds"
)

##Exploration interaction temp and body size 

spatial_with_trait$summary.fixed
spatial_no_trait$summary.fixed

temp_seq <- seq(min(df$temp_s), max(df$temp_s), length.out = 50)
length_seq <- seq(min(df$length_cm_s), max(df$length_cm_s), length.out = 50)

grid <- expand.grid(temp_s = temp_seq, length_cm_s = length_seq)
grid$intercept <- 1
grid$bathy_s <- mean(df$bathy_s)   # fijar otras variables en su media
grid$time <- mean(df$time)

X <- model.matrix(~ -1 + intercept + bathy_s + temp_s + length_cm_s + temp_s:length_cm_s, data = grid)
beta <- spatial_with_trait$summary.fixed$mean
grid$eta <- as.vector(X %*% beta)
grid$prob <- 1 / (1 + exp(-grid$eta))  # probabilidad binomial

ggplot(grid, aes(x = temp_s, y = length_cm_s, fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "Temperature (scaled)", y = "Length (scaled)", fill = "Prob of presence") +
  theme_minimal()


df$temp_s <- (df$temp - mean(df$temp)) / sd(df$temp)
df$length_cm_s <- (df$length_cm - mean(df$length_cm)) / sd(df$length_cm)
mean_temp <- mean(df$temp)
sd_temp   <- sd(df$temp)
mean_length <- mean(df$length_cm)
sd_length   <- sd(df$length_cm)
temp_orig_seq <- seq(min(df$temp), max(df$temp), length.out = 50)
length_orig_seq <- seq(min(df$length_cm), max(df$length_cm), length.out = 50)
grid_orig <- expand.grid(
   temp = temp_orig_seq,
   length_cm = length_orig_seq)
grid_orig$temp_s <- (grid_orig$temp - mean_temp) / sd_temp
grid_orig$length_cm_s <- (grid_orig$length_cm - mean_length) / sd_length
grid_orig$intercept <- 1
grid_orig$bathy_s <- mean(df$bathy_s)  
grid_orig$time <- mean(df$time)
X <- model.matrix(~ -1 + intercept + bathy_s + temp_s + length_cm_s + temp_s:length_cm_s, data = grid_orig)
beta <- spatial_with_trait$summary.fixed$mean
grid_orig$eta <- as.vector(X %*% beta)
grid_orig$prob <- 1 / (1 + exp(-grid_orig$eta))

prob_len <- ggplot(grid_orig, aes(x = temp, y = length_cm, fill = prob)) +
 geom_tile() +
 scale_fill_viridis_c(option = "magma") +
 labs(x = "Bottom Temperature (°C)", y = "Mean body size (cm)", fill = "Probability of presence") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )
