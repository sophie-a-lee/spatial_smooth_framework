##############################################
####                                      ####
####  Case study: Dengue in South Brazil  ####
####                                      ####
##############################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load human movement-based connectivity coordinates ####
## NOTE: this code will only work if 00_human_movement_coords.R has been run
connect_coords <- fread("data/human_movement_coords.csv") %>% 
  dplyr::select(municip_code_ibge, connect_coord1, connect_coord2)


#### Load dengue data for South Brazil ####
df <- fread("data/dengue_south_brazil.csv") %>% 
  group_by(municip_code_ibge) %>% 
  # Find average number of cases, round to integer for model
  summarise(dengue_mean = round(mean(dengue_year), digits = 0),
            pop_mean = mean(pop_year)) %>% 
  ungroup() %>% 
  mutate(DIR_mean = (dengue_mean/pop_mean) * 10^5,
         # Offset to convert cases to incidence per 100,000 residents
         E = pop_mean/10^5) %>% 
  full_join(., shp, by = "municip_code_ibge") %>% 
  full_join(., connect_coords, by = "municip_code_ibge") %>% 
  st_as_sf()


#### Fit data using spatial smooth model ####
## Need reasonable estimate of theta + lambda for negative binomial model
gam_nb <- gam(dengue_mean ~ s(lon, lat, k = 10, bs = "tp") +
                s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
              offset = log(E), family = nb, data = df, method = "REML")


# Theta estimate
gam_nb$family$getTheta(TRUE)

# Lambda estimate
gam_nb$sp


#### Fit 3 phi model ####
# Extract spatial smosoth terms
jd_smooth <- jagam(dengue_mean ~ s(lon, lat, k = 10, bs = "tp") +
                     s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
                   offset = log(E),
                   file = "jagam_smooth.txt", 
                   # Set prior for smoothing function (in this case, the lon/lat spatial smooth)
                   sp.prior = "gamma", 
                   family = poisson,
                   diagonalize = F, 
                   data = df)


# Extract basis functions to use as linear predictors in the model
X_dist <- jd_smooth$jags.data$X[,2:10]
X_human <- jd_smooth$jags.data$X[,11:19]

# Set constants for model (number obs +  number coefficients)
Consts <- list(n = nrow(df), m_dist = ncol(X_dist),
               m_human = ncol(X_human), m = (ncol(X_dist) + ncol(X_human)))


## Write model formula
Model <- nimbleCode({ 
  
  # u_dist = spatial smooth term based on distance
  u_dist[1:n] <- X_dist[1:n, 1:m_dist] %*% b[2:(m_dist + 1)] 
  
  # u_human = spatial smooth term based on human movement
  u_human[1:n] <- X_human[1:n, 1:m_human] %*% b[(m_dist + 2):(m + 1)] 
  
  for (i in 1:n) { 
    # y = number of cases
    y[i] ~ dnegbin(p[i], theta)
    
    log(mu[i]) <- b[1] + u_dist[i] + u_human[i] + v[i] + log(e[i])
    
    p[i] <- theta/(theta + mu[i])
    
    # v = iid random effect
    v[i] ~ dnorm(0, sd = sig_re)
  } 
  
  # Priors
  # Scaling parameter
  theta ~ T(dgamma(.05,.005), 0, 10)
  
  # Random effect SD
  sig_re ~ dexp(.1)
  
  # Intercept
  b[1] ~ dnorm(0, sd = 5) 
  
  ## prior for sd(s(lon,lat))
  K1[1:(m_dist),1:(m_dist)] <- 
    S1[1:(m_dist),1:(m_dist)] * lambda[1] + 
    S1[1:(m_dist), (m_dist + 1):(2*m_dist)] * lambda[2]
  
  ## prior for sd(s(connect_coord1, connect_coord2))
  K2[1:m_human,1:m_human] <- 
    S2[1:m_human, 1:m_human] * lambda[3] + 
    S2[1:m_human, (m_human + 1):(2*m_human)] * lambda[4]
  
  # Prior for smooth coefficient
  b[2:(m_dist + 1)] ~ dmnorm(zero[1:m_dist], 
                             K1[1:m_dist, 1:m_dist]) 
  
  b[(m_dist + 2):(m + 1)] ~ dmnorm(zero[1:m_human], 
                                   K2[1:m_human, 1:m_human]) 
  
  ## smoothing parameter priors 
  for (i in 1:4) {
    lambda[i] ~ T(dgamma(.05,.005), 0, 100)
  }
} )


# Convert jagam data into data suitable for nimble
nimbleData <- list(y = df$dengue_mean, X_dist = X_dist, 
                   X_human = X_human,
                   zero = jd_smooth$jags.data$zero,
                   S1 = jd_smooth$jags.data$S1,
                   S2 = jd_smooth$jags.data$S2, e = df$E)


# Set initial values for MCMC
inits <- list(b = rnorm((ncol(X_dist) + ncol(X_human) + 1), sd=0.1), 
              v =  rnorm(nrow(df), 1), theta = 1,
              lambda = rep(3, 4), sig_re = .5)

# Sets up model in nimble code
nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                           constants = Consts, data = nimbleData, 
                           inits = inits)

# Tell model which parameter to estimate and return
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("b","lambda", "theta",
                                                   "u_dist", "u_human",
                                                   "v", "sig_re", "mu"))


modelMCMC <- buildMCMC(MCMCconfig)
compiled_model <- compileNimble(nimbleModel)
compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
results <- runMCMC(compiled_model_MCMC, niter = 50000, nburnin = 10000, 
                   nchains = 3, thin=10, inits=inits, progressBar = T, 
                   samplesAsCodaMCMC = T)


# Check convergence
GR.diag <- gelman.diag(results, multivariate = F)

sum(GR.diag$psrf[,"Point est."] > 1.1)


## Save model
write_rds(results, file = "output/case_study/model_output.rds")


#### Extract predicted intercept ####
b_est <- b_extract_3phi(results) 



#### Extract estimates for phi/mixing parameters ####
model_vars <- extract_vars_3re(results)

## Plot density of mixing parameter simulations
phi_density <- ggplot(data = model_vars) +
  geom_density(aes(x = propn_var_dist, ..scaled..)) +
  geom_density(aes(x = propn_var_human, ..scaled..), col = "blue") +
  geom_density(aes(x = propn_var_iid, ..scaled..), col = "red") +
  labs(x = "Proportion of variance", y = "Density") +
  theme_bw()


ggsave(phi_density,
       filename = "output/case_study/phi_density.png")


## Estimate mean and 95% CI of mixing parameters and format into a table 
phi_ests <- model_vars %>% 
  summarise(phi_dist = mean(propn_var_dist), 
            phi_dist_lci = quantile(propn_var_dist, .025),
            phi_dist_uci = quantile(propn_var_dist, .975),
            phi_human = mean(propn_var_human), 
            phi_human_lci = quantile(propn_var_human, .025),
            phi_human_uci = quantile(propn_var_human, .975),
            phi_iid = mean(propn_var_iid), 
            phi_iid_lci = quantile(propn_var_iid, .025),
            phi_iid_uci = quantile(propn_var_iid, .975)) %>% 
  transmute(phi_dist_fmt = paste0(round(phi_dist, 3), " (",
                                  round(phi_dist_lci, 3), ", ",
                                  round(phi_dist_uci, 3), ")"),
            phi_human_fmt = paste0(round(phi_human, 3), " (",
                                   round(phi_human_lci, 3), ", ",
                                   round(phi_human_uci, 3), ")"),
            phi_iid_fmt = paste0(round(phi_iid, 3), " (",
                                 round(phi_iid_lci, 3), ", ",
                                 round(phi_iid_uci, 3), ")"))


fwrite(phi_ests, "output/case_study/phi_ests.csv")


#### Extract random  effect estimates ####
df_re <- re_est_3phi(results)


# Ditance-based random effects
dist_re_map <- ggplot(data = df_re) +
  geom_sf(aes(fill = spat_dist_re), lwd = .05) +
  scale_fill_viridis_c(name = expression(u[1])) +
  theme_void()

ggsave(dist_re_map, 
       filename = "output/case_study/dist_re_map.png")


# Human movement-based random effects
human_re_map <- ggplot(data = df_re) +
  geom_sf(aes(fill = spat_human_re), lwd = .05) +
  scale_fill_viridis_c(name = expression(u[2])) +
  theme_void()

ggsave(human_re_map, 
       filename = "output/case_study/human_re_map.png")


# Unstructured random effects
iid_re_map <- ggplot(data = df_re) +
  geom_sf(aes(fill = iid_re), lwd = .05) +
  scale_fill_viridis_c(name = "v") +
  theme_void()

ggsave(iid_re_map,
       filename = "output/case_study/iid_re_map.png")


# Total random effects
total_re_map <- ggplot(data = df_re) +
  geom_sf(aes(fill = tot_re), lwd = .05) +
  scale_fill_viridis_c(name = expression(S[i])) +
  theme_void()

ggsave(total_re_map,
       filename = "output/case_study/total_re_map.png")
