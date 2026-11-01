library(tidyverse)
library(greta)
library(FNN)
library(ggforce)

# Prediction was done at single plot level due to computing limitations

# Read model and draws ----------------------------------------------------

draws <- read_rds("outputs/draws.rds")
mod   <- read_rds("outputs/mod.rds")


# Read data ----------------------------------------------------
dat <- read_csv( "data/data.csv") 

# Litter data
litter <- dat %>%
  distinct(Plot, Collection, Plot_Trap, Trap, Species_code, Mass, Interval_days, X_trap, Y_trap, Date, Sieve_size)

n_sp<- length(unique(litter$Species_code))
n_min <- 100
sp_to_use <- names(which(table(litter$Species_code) >= n_min))

# Tree locations, kauri dieback canopy scores and DBH size data
trees <- dat %>%
  distinct(Plot, tree_id, Species_code, X_tree, Y_tree, DBH, Canopy)  %>%
  filter(Species_code %in% sp_to_use) %>%
  rename(
    x = X_tree,
    y = Y_tree,
    tree = tree_id
  ) %>%
  mutate(Species_code = as.factor(Species_code),
         kdb = as.factor(Canopy),
         Species = as.factor(if_else(Species_code == "AGAaus", "Kauri", "Non-kauri")))



# New prediction grids ----------------------------------------------------
# Generate 1 m x 1 m grid in each forest plot, 
# can be thought of as new traps

new_traps <- 
  expand_grid(x_trap = seq(0.5, 39.5, 1),
              y_trap = seq(0.5, 49.5, 1),
              Species_code = sp_to_use,
              Plot = "HU",
              Month = 1:12) %>% 
  left_join(trees) %>% 
  mutate(dist = sqrt((x_trap - x)^2 + (y_trap - y)^2)) %>% 
  filter(dist <= 25) %>% 
  group_by(Plot, x_trap, y_trap, Species_code, Month) %>%
  mutate(group_id = cur_group_id()) %>% 
  ungroup

# aggregate new traps to per trap per species per plot per month

new_traps_agg <- 
  new_traps %>% 
  distinct(group_id, Plot, x_trap, y_trap, Species_code, Month) %>%
  group_by(Plot, x_trap, y_trap) %>% 
  mutate(trap_id = cur_group_id()) %>% 
  arrange(group_id) %>% 
  ungroup




# Predictions ----------------------------------------------

# recover parameters
alpha <- mod$target_greta_arrays$alpha
log_X0 <- mod$target_greta_arrays$log_X0
log_beta <- mod$target_greta_arrays$log_beta
log_epsilon <- mod$target_greta_arrays$log_epsilon
log_Xb <- mod$target_greta_arrays$log_Xb
Xb <- exp(log_Xb)
delta <- mod$target_greta_arrays$delta
theta <- mod$target_greta_arrays$theta

# species index ( _new suffix to denote new trap variables)
# to avoid name conflicts with original trap indexing
species <- as.factor(litter$Species_code)
sp_litter <- as.numeric(species)
sp_new <- as.numeric(factor(new_traps$Species_code, levels = levels(species))) 
sp_new_agg <- as.numeric(factor(new_traps_agg$Species_code, levels = levels(species))) 

# KDB indices
kdb_new <- new_traps$Canopy

# sieve size to 1mm, include sieve effects
sieve <- 1


log_months_new <- log(1)

# cosinor (temporal) component
month_new <- new_traps_agg$Month
period <- 12
C_new <- cos(2 * pi * month_new / period)
S_new <- sin(2 * pi * month_new / period)
X_cosinor_new <- cbind(C_new, S_new)

# distance-weighted neighbourhood effect
D_new <- new_traps$DBH
# add 1 cm (0.01 m) to dist to avoid zero distance, which leads to Inf when
# log-transformed
distance_new <- ((new_traps$dist + 0.01) / sd(dat$dist)) 
log_distance_new <- log(distance_new)

log_d_ratio_new <- log_distance_new - log_X0[sp_new] 

log_D_distance_new <- 
  exp(log_beta[sp_new] + log_epsilon[kdb_new]) * log(D_new) - 
  0.5 * (log_d_ratio_new / Xb[sp_new])^2

group_id <- new_traps$group_id

max_v <- tapply(log_D_distance_new, group_id, "max") 
v_n_minus_max_v <- log_D_distance_new - max_v[group_id]
exp_v_n_minus_max_v <- exp(v_n_minus_max_v)

log_sum_exp <- log(tapply(exp_v_n_minus_max_v, group_id, "sum")) 
log_D_distance_sum_new <- max_v + log_sum_exp

# spatial random and latent effects
# for each grid, it's random effect and latent variable are calculated
# using inverse-distance weighted average from the traps
e_trap <- mod$target_greta_arrays$e_trap
trap <- levels(as.factor(paste(litter$Plot, litter$Trap, sep = "_")))
trap_id <- new_traps_agg$trap_id

# weighted by distance so we need original litter trap and grid coordinates
trap_coord <- 
  litter %>% 
  group_by(Plot, Plot_Trap) %>% 
  summarise(x_trap = mean(X_trap),
            y_trap = mean(Y_trap)) %>% 
  ungroup %>% 
  mutate(Plot = as.factor(Plot),
         x_trap = x_trap + (as.numeric(Plot)-1) * 1000,
         y_trap = y_trap + (as.numeric(Plot)-1) * 1000) %>% 
  arrange(Plot_Trap)

# for the grids (new litter traps)
trap_coord_new <- 
  new_traps_agg %>% 
  distinct(Plot, trap_id, x_trap, y_trap) %>% 
  mutate(Plot = factor(Plot, levels = levels(as.factor(litter$Plot))),
         x_trap = x_trap + (as.numeric(Plot)-1) * 1000,
         y_trap = y_trap + (as.numeric(Plot)-1) * 1000) %>% 
  arrange(trap_id)

# limit each grid to K number of nearest litter traps
# Find k number of nearest original traps to each cell in the grid
k <- 9  # number of nearest neighbours 
nn <- get.knnx(data = trap_coord[, c("x_trap", "y_trap")], # Original traps
               query = trap_coord_new[, c("x_trap", "y_trap")], # To query for neighbours within all cells
               k = k)

# distances and indices of the nearest neighbours 
# 1 row is a single new_trap and column is for each nearest original trap
dist <- nn$nn.dist
idx <- nn$nn.index

# weights = inverse distance
w <- 1 / (dist + 1e-12) 
w <- w / rowSums(w)   


e_trap_new <- tapply(e_trap[idx] * as.vector(w), as.vector(row(idx)), "sum")

# interpolate latent variables
# Same as random effect
lvs <- mod$target_greta_arrays$lvs
lvs1_new <- tapply(lvs[idx, 1] * as.vector(w), as.vector(row(idx)), "sum")
lvs2_new <- tapply(lvs[idx, 2] * as.vector(w), as.vector(row(idx)), "sum")
lvs_new <- cbind(lvs1_new, lvs2_new)

loadings <- mod$target_greta_arrays$loadings
random_lv_new <- lvs_new %*% t(loadings)

# adjust trap areas
grid_area <- 1*1
trap_area <- pi * 0.3^2
area_offset <- grid_area / trap_area
log_area_offset <- log(area_offset)

# fitted value at the species level
log_mu_new <- 
  alpha[sp_new_agg] + 
  delta[sp_new_agg] * sieve +
  log_D_distance_sum_new + 
  rowSums(X_cosinor_new * theta[sp_new_agg, ]) +
  random_lv_new[cbind(trap_id, sp_new_agg)] +
  e_trap_new[trap_id] +
  log_months_new +
  log_area_offset

mu_new <- exp(log_mu_new)

# sum across species and month to get total annual plot litterfall
mu_agg <- tapply(mu_new, trap_id, "sum")

mu_sim <-
  calculate(mu_agg,
            values = draws,
            nsim = 100,
            trace_batch_size = 10)

# combine everything into a prediction dataframe
predictions <- 
  new_traps_agg %>% 
  distinct(trap_id, Plot, x_trap, y_trap) %>% 
  arrange(trap_id) %>% 
  mutate(
    estimate = apply(mu_sim$mu_agg, 2, median)
  )



# Plot for one of the plots, HU
p_HU <- predictions %>% 
  filter(Plot == "HU") %>% 
  ggplot() +
  geom_raster(aes(x_trap, y_trap, fill = estimate)) +
  scale_fill_viridis_c(option = "turbo") +
  geom_circle(
    data = trees %>% filter(Plot == "HU"),
    aes(x0 = x, y0 = y, 
        linetype = kdb, 
        r = DBH / 100, color = Species),
    fill = NA
  ) +
  scale_colour_manual(values = c("black", "darkgrey"), name = "Species") +
  scale_linetype_manual( values = c("solid", "longdash", "dotdash", "dotted"),
                         name = "Canopy Health",
                         labels = c("Healthy", "Canopy thinning", "Branch dieback", "Severe dieback")) +
  coord_fixed() +
  labs(title = "Mid-successional early-symptomatic") +
  theme_minimal(base_size = 10) +
  guides(fill = guide_colorbar(order = 1),   # for continuous fill
         colour = guide_legend(order = 2),
         linetype = guide_legend(order = 3)) +
  theme(
    legend.title = element_text(size = 14),    # Legend title font size
    legend.text = element_text(size = 12),     # Legend item font size
    legend.key.size = unit(1.2, "cm"),         # Size of legend keys (squares/lines)
    legend.spacing = unit(0.3, "cm")           # Space between legend items
  )

