library(tidyverse)
library(lubridate)
library(greta)

# Bayesian litterfall model

# Read data ----------------------------------------------------------------
# one row per litter observation x neighbouring tree (within 25 m) x collection
dat <- read_csv("data/data.csv", show_col_types = FALSE)

# litter data, one row per litter observation
litter <- dat %>%
  distinct(Plot, Collection, Trap, Species_code, Mass, Interval_days, Date, Sieve_size)


# Response and offset ------------------------------------------------------
L <- litter$Mass
log_months <- log(litter$Interval_days / (365 / 12))


# Indices ------------------------------------------------------------------
# species of litter and of trees
species <- as.factor(litter$Species_code)
sp_litter <- as.numeric(species)
sp_tree <- as.numeric(factor(dat$Species_code, levels = levels(species)))
n_sp <- length(unique(species))

# kauri dieback canopy score
kdb <- dat$Canopy
n_kdb <- length(unique(kdb))

# trap
trap <- as.numeric(as.factor(paste(litter$Plot, litter$Trap, sep = "_")))
n_trap <- length(unique(trap))

# cosinor (temporal) component
month_num <- month(litter$Date)
period <- 12
X_cosinor <- cbind(C = cos(2 * pi * month_num / period),
                   S = sin(2 * pi * month_num / period))

# index litter of a species collected in a collection from a trap of a plot
# and match it to possible trees of that species within 25 m
plot_col_trap_sp <- as.factor(paste(litter$Plot, litter$Collection, litter$Trap,
                                    litter$Species_code, sep = "_"))
plot_col_trap_sp_litter <- as.numeric(plot_col_trap_sp)

plot_trap_sp_tree <- as.numeric(factor(
  paste(dat$Plot, dat$Collection, dat$Trap, dat$Species_code, sep = "_"),
  levels = levels(plot_col_trap_sp)))


# Priors -------------------------------------------------------------------
# species-level terms: intercept, DBH scaling, kernel centre and spread,
# sieve effect, seasonality
alpha <- normal(0, 2.5, dim = n_sp)
log_beta <- normal(0, 0.5, dim = n_sp)
log_X0 <- normal(0, 0.5, dim = n_sp)
log_Xb <- normal(0, 0.5, dim = n_sp)
X0 <- exp(log_X0)
Xb <- exp(log_Xb)
delta <- normal(0, 1, dim = n_sp)
theta <- normal(0, 1, dim = c(n_sp, ncol(X_cosinor)))

# kauri dieback disease effect (1 is healthy, fixed at 0)
log_epsilon <- zeros(n_kdb)
log_epsilon[-1] <- normal(0, 1, dim = n_kdb - 1)

# trap random effect
sd_trap <- exponential(1)
z_trap <- normal(0, 1, dim = n_trap)
e_trap <- z_trap * sd_trap


# Distance-weighted neighbourhood effect -----------------------------------
D <- dat$DBH
distance <- dat$dist / sd(dat$dist)
log_distance <- log(distance)

# lognormal distance kernel
log_d_ratio <- log_distance - log_X0[sp_tree]
log_D_distance <- exp(log_beta[sp_tree] + log_epsilon[kdb]) * log(D) -
  0.5 * (log_d_ratio / Xb[sp_tree])^2

# sum over neighbouring trees with log-sum-exp
# https://mc-stan.org/docs/stan-users-guide/floating-point.html#log-sum-of-exponentials
max_v <- tapply(log_D_distance, plot_trap_sp_tree, "max")
v_n_minus_max_v <- log_D_distance - max_v[plot_trap_sp_tree]
exp_v_n_minus_max_v <- exp(v_n_minus_max_v)
log_sum_exp <- log(tapply(exp_v_n_minus_max_v, plot_trap_sp_tree, "sum"))
log_D_distance_sum <- max_v + log_sum_exp


# Sieve size and latent variables ------------------------------------------
# sieve size (1 = 1 mm)
sieve <- as.numeric(litter$Sieve_size == "1mm")

# latent variables for edge effects (missing neighbours)
n_lv <- 2
loadings <- zeros(n_sp, n_lv)
diag(loadings) <- normal(0, 2.5, dim = n_lv, truncation = c(0, Inf))
lower_idx <- lower.tri(loadings)
loadings[lower_idx] <- normal(0, 2.5, dim = sum(lower_idx))

lvs <- normal(0, 1, dim = c(n_trap, n_lv))
random_lv <- lvs %*% t(loadings)
random_lv <- random_lv[cbind(trap, sp_litter)]


# Model --------------------------------------------------------------------
log_mu <-
  alpha[sp_litter] + log_D_distance_sum +
  delta[sp_litter] * sieve +
  rowSums(X_cosinor * theta[sp_litter, ]) +
  random_lv +
  e_trap[trap] +
  log_months

# litter mass is lognormal
sigma <- exponential(1)
distribution(L) <- lognormal(log_mu, sigma)

mod <- model(alpha, log_beta, log_epsilon, delta, theta, sd_trap, e_trap, sigma,
             loadings, lvs, log_X0, log_Xb)

draws <- mcmc(mod,
              warmup = 3000,
              sampler = hmc(20, 25),
              chains = 1)

write_rds(draws, "outputs/draws.rds")
write_rds(mod, "outputs/mod.rds")
