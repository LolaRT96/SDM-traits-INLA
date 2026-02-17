# ==============================================================================
# Title: Simulated Multispecies Spatial SDM with Trait-Modulated Response (INLA)
# Author: M. Grazia Pennino
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

# Add spatial and yearly noise
linpred_spatial <- linpred + spatial_field + rnorm(n_points, 0, 0.3)

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

# Model with trait-modulated slope (random slope on temp)
# What mean?
# Species differ in the thermal responses but is not dependent of the trait 
# Trait is not added as a modular effect just additive 

#Esta modificada con respecto a la de Maria por que a ella le salia por duplicado por que 
#cada observacion tenia su propia pendiente, nosotras queremos un slope por especie 

formula_trait <- presence ~ temp_s + depth_s + trait_s +
  f(species, model = "iid") +
  f(species_slope_id, temp_s, model = "iid") +
  f(year, model = "iid") +
  f(spatial, model = spde)

# Model with trait*temp
# What mean?
# Thermal sensitivity depends on the trait 
# It not capture the noise between species 

formula_trait2 <- presence ~ -1 + temp_s * trait_s + depth_s +
  f(species, model = "iid") +
  f(year, model = "iid") +
  f(spatial, model = spde)

#Mixed formula 2 and 3
# What mean?
# We have a interaction between species and temperature 
# and we also have a residual slope by species 
# the thermal niche depends on the trait with variability between species 

# formula_trait3 <- presence ~  -1 + temp_s * trait_s + depth_s +
#   f(species, model = "iid") +
#   f(species_slope_id, temp_s, model = "iid",
#     group = species_slope_id, control.group = list(model = "iid")) +
#   f(year, model = "iid") +
#   f(spatial, model = spde)
#Usando group= en la tercera parte del modelo, se asocia un efecto aleatorio para cada OBSERVACIÓN 


formula_trait3_fix <- presence ~ -1 + temp_s * trait_s + depth_s + 
  f(species, model = "iid") + 
  f(species_slope_id, temp_s, model = "iid") + 
  f(year, model = "iid") + 
  f(spatial, model = spde)

#En la tercera parte del modelo 

model_nt <- inla(formula_nt, family = "binomial",
                 data = inla.stack.data(stack),
                 control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE),
                 verbose = TRUE)

model_trait <- inla(formula_trait, family = "binomial",
                    data = inla.stack.data(stack),
                    control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE),
                    verbose = TRUE)

model_trait2 <- inla(formula_trait2, family = "binomial",
                    data = inla.stack.data(stack),
                    control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE),
                    verbose = TRUE)

model_trait3 <- inla(formula_trait3_fix, family = "binomial",
                     data = inla.stack.data(stack),
                     control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE),
                     verbose = TRUE)


#Compare model fits

comparison <- tibble::tibble(
  Model = c("Spatial_NoTrait", "Spatial_Trait", "Spatial_Trait2", "Spatial_Trait3"),
  DIC   = c(model_nt$dic$dic, model_trait$dic$dic, model_trait2$dic$dic, model_trait3$dic$dic),
  WAIC  = c(model_nt$waic$waic, model_trait$waic$waic, model_trait2$waic$waic, model_trait3$waic$waic)
)

print(comparison)


# --- 7. Evaluate performance --------------------------------------------------

idx <- inla.stack.index(stack, "est")$data
pred_df <- sim_df@data
pred_df$pred_nt <- model_nt$summary.fitted.values$mean[idx]
pred_df$pred_trait <- model_trait$summary.fitted.values$mean[idx]
pred_df$pred_trait2 <- model_trait2$summary.fitted.values$mean[idx]
pred_df$pred_trait3 <- model_trait3$summary.fitted.values$mean[idx]

roc_nt <- roc(pred_df$presence, pred_df$pred_nt)
roc_trait <- roc(pred_df$presence, pred_df$pred_trait)
roc_trait2 <- roc(pred_df$presence, pred_df$pred_trait2)
roc_trait3 <- roc(pred_df$presence, pred_df$pred_trait3)

roc_vals <- tibble::tibble(
  Model = comparison$Model,
  AUC   = c(
    auc(roc(pred_df$presence, pred_df$pred_nt)),
    auc(roc(pred_df$presence, pred_df$pred_trait)),
    auc(roc(pred_df$presence, pred_df$pred_trait2)),
    auc(roc(pred_df$presence, pred_df$pred_trait3))
  )
)

print(roc_vals)


# --- 8. Spatial field plots ---------------------------------------------------

projr <- inla.mesh.projector(mesh, dims = c(200, 200))
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
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(
    aes(z = value),
    stroke = 0.15,
    size = 3,
    skip = 0      # <- etiqueta TODAS las líneas
  ) +
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
p_wt_multiSimu <- ggplot(df_wt_multiSimu, aes(x, y, fill = value)) +
  geom_raster() +
  geom_contour(aes(z = value), colour = "black", linewidth = 0.3) +
  geom_text_contour(
    aes(z = value),
    stroke = 0.15,
    size = 3,
    skip = 0      # <- etiqueta TODAS las líneas
  ) +
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

windows();(p_nt_multiSimu | p_wt_multiSimu)

# --- 9. Trait vs slope relationship -------------------------------------------

slopes <- model_trait3$summary.random$species_slope_id %>%
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

library(ggplot2)
library(ggrepel)
library(dplyr)
library(paletteer)

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
  scale_color_paletteer_c("viridis::turbo", name = "trait") +
  labs(
    x = "Temperature effect (slope)",
    y = "Species",
    title = "Species-Specific Sensitivity to Bottom Temperature"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )


