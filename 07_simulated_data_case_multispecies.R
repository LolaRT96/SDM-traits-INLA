# ==============================================================================
# Title: Simulated Multispecies Spatial SDM with Trait-Modulated Response (INLA)
# Author: M. Grazia Pennino & M.D. Riesgo
# Date: 2025-07-15
# Description:
#   - Simulate multispecies presence–absence data with trait-modulated temperature effects
#   - Fit spatial INLA models with/without trait interaction (slope variation)
#   - Evaluate model fit and plot spatial fields and trait-response relationships

# ==============================================================================

# --- 1. Load libraries --------------------------------------------------------

library(INLA)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggrepel)
library(pROC)
library(fields)
library(sp)
library(metR)
library(patchwork)
library(paletteer)
library(showtext)
library(sysfonts)

font_add("Helvetica", 
         regular = "~/Downloads/helvetica-255/Helvetica.ttf")
showtext_auto()

# --- 2. Simulate spatial domain and mesh -------------------------------------

set.seed(123)

n_points <- 1000

coords <- data.frame(
  x = runif(n_points, 0, 100),
  y = runif(n_points, 0, 100)
)

coordinates(coords) <- ~x + y

mesh <- inla.mesh.2d(loc = coords, max.edge = c(10, 15), cutoff = 0.8)
mesh$n 

plot(mesh, main = "SPDE Mesh")
points(coords, col = "red", pch = 16, cex = 0.5)


# --- 3. Simulate species, traits, and environment ----------------------------

n_species <- 6

species_id <- sample(1:n_species, n_points, replace = TRUE)
                     
years <- sample(2015:2022, n_points, replace = TRUE)

# Traits (mean body length) per species

#OLD CODE 
# trait_species <- rnorm(n_species, mean = 40, sd = 15) 

#More differences between traits 

trait_species <- seq(25, 85, length.out = n_species) + rnorm(n_species, 0, 3)


#There is no intravariation by species 
#in the same group (species 1) you will have the same trait 

#IT IS IMPORTANT TO UNDERSTAD THAT
#Thermal sensitivity is the same between individuals of the same species 
#but different between species 

trait <- trait_species[species_id]

# Environmental variables (scaled)
temp <- scale(rnorm(n_points, mean = 12, sd = 2))
depth <- scale(rnorm(n_points, mean = 200, sd = 50))

# Species-specific temperature sensitivity (slope ~ trait)
#OLD CODE 
#temp_slope <- 0.2 + 0.02 * trait_species + rnorm(n_species, 0, 0.05)
#NEW CODE
temp_slope <- 0.1 + 0.05 * trait_species + rnorm(n_species, 0, 0.03)
#Amplificacion de las variables ecologicas 

#Species with larger body size will have more thermal sensitivity 

linpred <- -1 + 0.5 * depth + temp_slope[species_id] * temp
#depth: fixed effect
#temp_slope[species_id] * temp: variable effect between species

# Add spatial field

spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                            prior.range = c(1, 0.01),
                            prior.sigma = c(1, 0.01))

s_index <- inla.spde.make.index("spatial", spde$n.spde)
A_matrix <- inla.spde.make.A(mesh = mesh, loc = coordinates(coords))

w <- rnorm(spde$n.spde, mean = 0, sd = 1)
spatial_field <- as.vector(A_matrix %*% w)

spatial_strength <- 2.5
noise_sd <- 0.2
 
linpred_spatial <- linpred + spatial_strength * spatial_field + rnorm(n_points, 0, noise_sd)

# # # Add spatial and yearly noise
# linpred_spatial <- linpred + spatial_field + rnorm(n_points, 0, 0.3)

# Simulate presence-absence

prob <- 1 / (1 + exp(-linpred_spatial)) #p=logit−1(η)
             
presence <- rbinom(n_points, 1, prob) #binary data bernuilli

# --- 4. Prepare dataframe -----------------------------------------------------

sim_df <- data.frame(
  presence = presence,
  temp = as.numeric(temp),
  depth = as.numeric(depth),
  x = coordinates(coords)[,1],
  y = coordinates(coords)[,2],
  species = as.factor(species_id),
  year = as.factor(years),
  trait = trait
)

# Scale all numerical variables (tem, depht, trait)
sim_df$trait_s <- scale(sim_df$trait)
sim_df$temp_s <- scale(sim_df$temp)
sim_df$depth_s <- scale(sim_df$depth)

# Species-specific ID for slope
sim_df$species_slope_id <- as.integer(sim_df$species)

#Visualize each species in the field 
sim_pres <- sim_df %>%
  filter(presence == 1)

ggplot(sim_pres, aes(x = x, y = y,
                     color = species,
                     shape = species)) +
  geom_point(size = 3, alpha = 0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "X",
    y = "Y",
    color = "Species",
    shape = "Species"
  )

summary_length_spp <- sim_df %>%
  group_by(species) %>%
  summarise(mean_size = mean(trait, na.rm = TRUE),
            mean_sd =sd(trait))

summary_length_spp

summary_length_spp <- as.data.frame(summary_length_spp)
# write.xlsx(summary_length_spp,
#            file = "~/SDMs_Traits/summary_length_spp.xlsx",
#            rowNames = FALSE)
# 

ggplot(sim_pres, aes(x = x, y = y,
                     color = trait,
                     shape = species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  coord_equal() +
  theme_classic() +
  labs(
    x = "X",
    y = "Y",
    color = "Size",
    shape = "Species"
  )

sim_pres$trait_f <- as.factor(sim_pres$trait)


ggplot(sim_pres, aes(x = x, y = y,
                     color = trait_f,
                     shape = species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_d() +
  coord_equal() +
  theme_classic() +
  labs(
    x = "X",
    y = "Y",
    color = "Size",
    shape = "Species"
  )



# --- 5. INLA stack ------------------------------------------------------------

coordinates(sim_df) <- ~x + y
A <- inla.spde.make.A(mesh = mesh, loc = coordinates(sim_df))

effects <- list(
  spatial = s_index,
  data.frame(
    intercept = 1,
    temp_s = sim_df$temp_s,
    depth_s = sim_df$depth_s,
    species = sim_df$species,
    year = sim_df$year,
    trait_s = sim_df$trait_s,
    species_slope_id = sim_df$species_slope_id
  )
)

stack <- inla.stack(
  data = list(presence = sim_df$presence),
  A = list(A, 1),
  effects = effects,
  tag = "est"
)

# --- 6. Fit models ------------------------------------------------------------

# Model with shared slope (no trait modulation)
# What mean?
# All species behave in the same way in temperature changes, no diffferences in the thermal niche
# Differences are only in the intercept

formula_nt <- presence ~ -1 + temp_s + depth_s +
  f(species, model = "iid") +
  f(year, model = "iid") +
  f(spatial, model = spde)

# Model with thermal-slope 
# What mean?
# We have a interaction between species and temperature 
# and we also have a residual slope by species 
# the thermal niche depends on the trait with variability between species 

formula_trait <- presence ~ -1 + temp_s * trait_s + depth_s + 
  f(species, model = "iid") + 
  f(species_slope_id, temp_s, model = "iid") + 
  f(year, model = "iid") + 
  f(spatial, model = spde)

# Model with thermal-slope but without spatial effect
formula_trait_Nospatial <- presence ~ -1 + temp_s * trait_s + depth_s + 
  f(species, model = "iid") + 
  f(species_slope_id, temp_s, model = "iid") + 
  f(year, model = "iid") 

#En la tercera parte del modelo 

model_nt <- inla(formula_nt, family = "binomial",
                 data = inla.stack.data(stack),
                 control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                 verbose = TRUE)

model_trait <- inla(formula_trait, family = "binomial",
                     data = inla.stack.data(stack),
                     control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     verbose = TRUE)

model_trait_Nospatial <- inla(formula_trait_Nospatial, family = "binomial",
                     data = inla.stack.data(stack),
                     control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     verbose = TRUE)

#Compare model fits

comparison <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial_Trait", "NonSpatial_Trait"),
  DIC   = c(model_nt$dic$dic,
            model_trait$dic$dic, 
            model_trait_Nospatial$dic$dic),
  WAIC  = c(model_nt$waic$waic, 
            model_trait$waic$waic, 
            model_trait_Nospatial$waic$waic),
  LPML_CPO = c(
    sum(log(model_nt$cpo$cpo), na.rm = TRUE),
    sum(log(model_trait$cpo$cpo), na.rm = TRUE),
    sum(log(model_trait_Nospatial$cpo$cpo), na.rm = TRUE)
  )
)

comparison <- comparison %>%
  mutate(
    DIC  = formatC(DIC, format = "f", digits = 2),
    WAIC = formatC(WAIC, format = "f", digits = 2),
    LPML_CPO = formatC(LPML_CPO, format = "f", digits = 2)
  )

print(comparison)


# --- 7. Evaluate performance --------------------------------------------------

idx <- inla.stack.index(stack, "est")$data
pred_df <- sim_df@data
pred_df$pred_nt <- model_nt$summary.fitted.values$mean[idx]
pred_df$pred_trait <- model_trait$summary.fitted.values$mean[idx]
pred_df$pred_trait_nonspatial<- model_trait_Nospatial$summary.fitted.values$mean[idx]

roc_nt <- roc(pred_df$presence, pred_df$pred_nt)
roc_trait <- roc(pred_df$presence, pred_df$pred_trait)
roc_trait_nospatial <- roc(pred_df$presence, pred_df$pred_trait_nonspatial)

roc_vals <- tibble::tibble(
  Model = comparison$Model,
  AUC   = c(
    auc(roc(pred_df$presence, pred_df$pred_nt)),
    auc(roc(pred_df$presence, pred_df$pred_trait)),
    auc(roc(pred_df$presence, pred_df$pred_trait_nonspatial))
  )
)

print(roc_vals)


# Plot calibration checks  -------------------------------------------------

idx <- inla.stack.index(stack, tag = "est")$data

calibration_data <- function(model, model_name, y, idx) {
  
  p <- model$summary.fitted.values$mean[idx]
  
  tibble(
    obs = y,
    pred = p
  ) %>%
    mutate(bin = dplyr::ntile(pred, 10)) %>%
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
  calibration_data(model_nt,               "NoTrait_Spatial",      sim_df$presence, idx),
  calibration_data(model_trait,            "Trait_Spatial",        sim_df$presence, idx),
  calibration_data(model_trait_Nospatial,  "Trait_NoSpatial",      sim_df$presence, idx)
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
  filename = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/plots/calplots_multiespeciesReal.jpg",
  plot = cal.plots,
  width = 1484,    # ancho en píxeles
  height = 758,   # alto en píxeles
  units = "px",
  dpi = 72        # dpi estándar para píxeles (72 dpi)
)


# --- 8. Spatial field plots ---------------------------------------------------

projr <- inla.mesh.projector(mesh, dims = c(300, 300))
field_trait <- inla.mesh.project(projr, model_trait$summary.random$spatial$mean)
field_nt <- inla.mesh.project(projr, model_nt$summary.random$spatial$mean)

df_nt_multiSimu <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_nt)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs

df_wt_multiSimu <- expand.grid(
  x = projr$x,
  y = projr$y
) %>%
  mutate(value = as.vector(field_trait)) %>%
  filter(!is.na(value))      # <- ELIMINA NAs

#No traits
p_nt_multiSimu <- ggplot(df_nt_multiSimu, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_scico(palette = "vik", ,  name = "Spatial effect") +
  labs(
    x = "x",
    y = "y"
  ) +
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


#With trait
p_wt_multiSimu <- ggplot(df_wt_multiSimu, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_scico(palette = "vik", ,  name = "Spatial effect") +
  labs(
    x = "x",
    y = "y"
  ) +
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

windows();(p_nt_multiSimu | p_wt_multiSimu)

combination <- (p_nt_multiSimu | p_wt_multiSimu)


ggsave(
  filename = "C:/Users/mdolores.riesgo/Documents/LolaR/PhD_MB/PhD_SideProjects/SDMs_Traits/plots/simulationMultiCase_change.jpg",
  plot = combination,
  width = 972,    # ancho en píxeles
  height = 380,   # alto en píxeles
  units = "px",
  dpi = 72        # dpi estándar para píxeles (72 dpi)
)


# --- 9. Trait vs slope relationship -------------------------------------------

slopes <- model_trait$summary.random$species_slope_id %>%
  as_tibble() %>%
  rename(
    species_id = ID,
    slope_mean = mean,
    lower = `0.025quant`,
    upper = `0.975quant`
  )

trait_df <- data.frame(
  species_id = 1:n_species,
  trait = trait_species
)

slope_trait <- left_join(slopes, trait_df, by = "species_id")

ggplot(slope_trait, aes(x = trait, y = slope_mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_smooth(method = "lm", color = "blue", linetype = "dashed") +
  labs(x = "Trait (mean body length)", y = "Temperature slope (θⱼ)",
       title = "Trait-Modulated Thermal Response") +
  theme_minimal()


slopes_species <- model_trait$summary.random$species_slope_id %>%
  as_tibble() %>%
  group_by(ID) %>%                
  summarise(
    SpeciesID  = unique(ID),
    mean_slope = mean(mean),
    lower      = mean(`0.025quant`),
    upper      = mean(`0.975quant`)
  )

trait_df <- data.frame(
  SpeciesID = 1:n_species,
  trait = trait_species
)

slope_trait <- left_join(slopes_species, trait_df, by = "SpeciesID")

slopes_with_trait <- slope_trait %>%
  arrange(trait) %>%  # ordenar de menor a mayor para que las grandes queden abajo
  mutate(SpeciesID = factor(SpeciesID, levels = SpeciesID))  # niveles en este orden

windows()

ggplot(slopes_with_trait, 
       aes(x = mean_slope, 
           y = SpeciesID,
           color = trait)) +  
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y",
                width = 0,
                linewidth = 1) +
  geom_point(size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4, color = "red") +
  scale_color_paletteer_c("viridis::turbo", name = "Body size") +
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
  )

thermal_slope <- ggplot(slopes_with_trait, 
       aes(x = mean_slope, 
           y = SpeciesID,
           color = factor(round(trait, 0)))) +  
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y",
                width = 0,
                linewidth = 1) +
  geom_point(size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             alpha = 0.4, color = "red") +
  scale_color_scico(palette = "imola", name = "Mean length (cm)") +
  labs(
    x = "Temperature effect (slope)",
    y = "Species"
  ) +
  theme_classic()+
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
    legend.text  = element_text(size = 14)
  )

thermal_slope 
