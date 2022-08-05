###############################################################
####                                                       ####
####  Script to load packages, data and create functions   ####
####                                                       ####
###############################################################


#### Load packages ####
# If pacman package is not already installed, run following line
# install.packages("pacman")
# library(pacman)

pacman::p_load(tidyverse, data.table, mgcv, sf, nimble, coda, mgcViz, INLA,
               spdep, spam, igraph, sfnetworks, lwgeom, cowplot, pROC)
sf_use_s2(F)


## Set seed to reproduce analysis 
set.seed(123)


#### Load data ####
## Shapefile for South Brazil
shp <- read_rds("data/south_brazil_shp.rds")


## Population data with coordinates
pop <- fread("data/dengue_south_brazil.csv") %>% 
  # Select 2018 to match with REGIC study
  filter(year == 2018) %>% 
  dplyr::select(municip_code_ibge, pop_year, lon, lat) %>% 
  #  Convert into shapefile using coordinates (to calculate distance)
  st_as_sf(coords = c("lon", "lat"))


#### Function opposite of %in% ####
'%!in%' <- function(x,y)!('%in%'(x,y))


##### Functions to simulate data #####
#### Function to create a smooth surface (see ?te) ####
## x, z = scaled connectivity coordinates
smooth_create <- function(x, z, sx = 0.3, sz = 0.4) { 
  
  x <- x*20
  
  (pi**sx*sz)*(1.2*exp(-(x - 0.2)^2/sx^2-(z - 0.3)^2/sz^2)+
                 0.8*exp(-(x - 0.7)^2/sx^2-(z - 0.8)^2/sz^2))
  
}


#### Function to simulate count data with single source of spatial structure ####
simulate_bym2 <- function(phi, data = df) {
  
  ## Return number of values
  n <- nrow(df)
  
  # Calculate lambda by combining smooth and IID random terms
  lam <- exp(sqrt(phi)*data$smooth +  
               sqrt(1-phi)*data$iid + log(data$E))
  
  # Simulate data from a Poisson distribution
  true_y <- rpois(n, lam)
  
  # Combine into a data frame
  df_sim <- data %>% 
    mutate(y = true_y,
           lambda = lam,
           phi_true = phi,
           true_SMR = y/E) %>% 
    arrange(municip_code_ibge) 
  
  return(df_sim)
}


#### Function to simulate count data with 2 sources of spatial structure ####
simulate_3phi <- function(phi, data = df) {
  
  ## Return number of values
  n <- nrow(df)
  
  # Calculate lambda
  lam <- exp(sqrt(phi$phi_dist)*data$smooth_dist +  
               sqrt(phi$phi_human)*data$smooth_human +
               sqrt(phi$phi_iid)*data$iid + log(data$E))
  
  #  Simulate data from a Poisson distribution
  true_y <- rpois(n, lam)
  
  # Combine into a data 
  df <- data %>% 
    mutate(phi_dist_true = phi$phi_dist,
           phi_human_true = phi$phi_human,
           y = true_y,
           lambda = lam,
           true_SMR = y/E) 
  
  return(df)
}


#### Function to simulate binary data with single source of spatial structure ####
simulate_binomial <- function(phi, data = df) {
  
  # Calculate logit(p)
  logit_p <- exp(sqrt(phi)*data$smooth +  
                   sqrt(1 - phi)*data$iid)
  
  p_sim <- logit_p/(1 + logit_p)
  
  sim <- list(NA)
  
  #  Simulate data from a binomial distribution
  # To compare results to true outcome, perform 20 draws from bernoulli per municipality
  for(i in 1:20) {
    true_y <- rbinom(n, 1, p_sim)
    
    sim[[i]] <- shp %>% 
      mutate(event = i,
             y = true_y,
             S = log(logit_p),
             p_true = p_sim,
             phi_true = phi)
    
  }
  
  df_sim <- reduce(sim, rbind)
  
  return(df_sim)
}


#### Function to simulate binary data with 2 sources of spatial structure ####
simulate_binomial_3phi <- function(phi, data = df) {

  # Calculate logit(p)
  logit_p <- exp(sqrt(phi$phi_dist)*data$smooth_dist +
                   sqrt(phi$phi_human)*data$smooth_human +
                   sqrt(phi$phi_iid)*data$iid)
  
  p_sim <- logit_p/(1 + logit_p)
  
  sim <- list(NA)
  
  #  Simulate data from a binomial distribution
  # To compare results to true outcome, perform 20 draws from bernoulli per municipality
  for(i in 1:20) {
    true_y <- rbinom(n, 1, p_sim)
    
    sim[[i]] <- df %>% 
      mutate(event = i,
             y = true_y,
             p_true = p_sim,
             phi_dist = phi$phi_dist,
             phi_human = phi$phi_human,
             phi_iid = phi$phi_iid)
    
  }
  
  df_sim <-reduce(sim, rbind)
  
  return(df_sim)
  
}


#### Functions to run models ####
#### Function to fit Poisson spatial smooth + iid model ####
smooth_iid_fit <- function(data) {
  
  # Use mgcv to extract basis functions for the smooth term
  jd_smooth <- jagam(y ~  s(lon, lat, k = 10, bs = "tp"),
                     offset = log(E),
                     # Save JAGS code to create model
                     file = "jagam_ten.txt", 
                     # Set prior for smoothing function (in this case, the spatial smooth)
                     sp.prior = "gamma", 
                     family = poisson,
                     # If T, smooths are re-parameterised to have an iid Gaussian prior (not appropriate here)
                     # diagonalize cannot be used with >1 dimension
                     diagonalize = F, 
                     data = data)
  
  
  # Extract basis functions to use as linear predictors in the model
  X <- jd_smooth$jags.data$X
  
  # Set constants for model (number obs & number of coefficients)
  Consts <- list(n = length(data$y), m = ncol(X))
  
  
  ## Write model formula
  Model <- nimbleCode({ 
    
    # u = spatial smooth term (using basis funtions)
    u[1:n] <- X[1:n, 2:m] %*% b[2:m] 
    
    for (i in 1:n) { 
      # y = number of cases
      y[i] ~ dpois(mu[i]) 
      
      log(mu[i]) <- b[1] + u[i] + v[i] + log(e[i])
      
      # v = iid random effect
      v[i] ~ dnorm(0, sd = sig_re)
    } 
    
    # Priors
    # Random effect SD
    sig_re ~ dexp(.1)
    
    # Intercept
    b[1] ~ dnorm(0, sd = 5) 
    
    ## prior for sd(smooth function)
    K1[1:(m-1), 1:(m-1)] <- S1[1:(m-1), 1:(m-1)] * lambda[1] + 
      S1[1:(m-1), m:(2*(m-1))] * lambda[2]
    
    # Prior for smooth coefficients
    b[2:m] ~ dmnorm(zero[2:m], K1[1:(m-1), 1:(m-1)]) 
    
    ## smoothing parameter priors 
    for (i in 1:2) {
      # truncate lambdas to avoid simulations getting stuck
      lambda[i] ~ T(dgamma(.05, .005), 0, 5)
    }
  } )
  
  
  # Convert jagam data into data suitable for nimble
  nimbleData <- list(y = data$y, X = X, zero = jd_smooth$jags.data$zero, 
                     S1 = jd_smooth$jags.data$S1, e = data$E)
  
  
  # Set initial values for MCMC
  inits <- list(b = rnorm(ncol(X),sd = 1), lambda = c(3, 3), 
                sig_re = runif(1), v = rnorm(nrow(data), 1))
  
  # Sets up model in nimble code
  nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                             constants = Consts, data = nimbleData, 
                             inits = inits)
  
  # Tell model which parameter to estimate and return
  MCMCconfig <- configureMCMC(nimbleModel,
                              monitors=c("b","lambda", "u", "v",
                                         "sig_re", "mu"),
                              # Return WAIC to compare models
                              enableWAIC = TRUE)
  
  
  # Build the model
  modelMCMC <- buildMCMC(MCMCconfig)
  
  compiled_model <- compileNimble(nimbleModel)
  
  compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
  
  results <- runMCMC(compiled_model_MCMC, thin = 10, 
                     niter = 100000, nburnin = 50000, 
                     nchains = 3, inits=inits, progressBar = T, 
                     samplesAsCodaMCMC = T, WAIC = TRUE)
  
  
  return(results)
  
}


#### Function to fit Poisson INLA model with BYM2 structure ####
inla_fit <- function(df) {
  
  # Set prior using PC priors
  phi.prior = list(prior = "pc", 
                   param = c(0.5, 2/3)) # P(phi < 0.5) = 2/3 
  prec.prior = list(prior = "pc.prec",
                    param = c(1, 0.01)) # P(st. dev > 1) = 0.01
  
  
  inla_model <- inla(y ~ f(index, model = "bym2", 
                           # Input binary  neighbourhood matrix
                           graph = "output/sim_study1/inla_nb.graph",
                           # Set priors
                           hyper = list(prec = prec.prior,
                                        phi = phi.prior)),
                     data = df, family = "poisson", offset = log(E),
                     # Return WAICfor model comparion
                     control.compute = list(config = T, waic = T),
                     inla.mode = "experimental")
  
  return(inla_model)
}


#### Function to fit Poisson INLA model with BYM2 structure ####
# Requires neighbourhood structure (in WB format) and scaling parameter
BYM2_nimble_fit <- function(data, nb = nb.WB, nb.scale = scale) {
  
  ## Set constants for NIMBLE
  Consts <-list(N = nrow(data),    
                L = length(nb$weights), 
                scale = nb.scale)
  
  ## Enter data for model
  nimbleData <- list(y = data$y,
                     E = data$E,                      
                     # elements of neighboring matrix
                     adj = nb$adj, 
                     weights = nb$weights, 
                     num = nb$num)
  
  
  # Write model formula
  BYM2Code <- nimbleCode(
    {
      for (i in 1:N){
        
        y[i] ~ dpois(mu[i])    
        
        log(mu[i]) <- log(E[i]) + b + S[i]
        
        # Define BYM2 structure 
        S[i] <- (1/sqrt(tau.b))*(sqrt((1-phi))*v[i] + sqrt(phi/scale)*u[i])
        
        v[i] ~ dnorm(0, tau = 1)           # IID RE   
      } 
      
      # ICAR RE
      u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 
                           tau = 1, zero_mean = 1) # its scaled so tau = 1
      
      # Priors
      # Intercept
      b ~ dnorm(0, sd = 5) 
      
      # prior for the precision of S
      tau.b ~ T(dgamma(1, 0.01), 0, 100)          
      
      # prior for the mixing parameter
      phi ~ T(dbeta(.5, .5), 0, 1)                 
    }
    
  )
  
  
  N <- nrow(data)
  
  ## Set intial values
  inits <- list(b = rnorm(1), 
                v = rnorm(N, sd = 1), 
                u = rnorm(N, sd = 1), 
                tau.b = runif(1), 
                phi = .5)            
  
  # Specify parameters to return
  params <- c("y", "mu", "S", "b", "u", "v", "phi", "tau.b")
  
  # Set up model in nimble code
  nimbleModel <- nimbleModel(code = BYM2Code, name = 'nimbleModel', 
                             constants = Consts, data = nimbleData, 
                             inits = inits)
  
  # Tell model which parameter to estimate and return
  MCMCconfig <- configureMCMC(nimbleModel,
                              monitors= params,
                              # Return WAIC to compare models
                              enableWAIC = TRUE)
  
  
  # Build the model
  modelMCMC <- buildMCMC(MCMCconfig)
  
  compiled_model <- compileNimble(nimbleModel)
  
  compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
  
  results <- runMCMC(compiled_model_MCMC, thin = 10, 
                     niter = 100000, nburnin = 50000, 
                     nchains = 3, inits=inits, progressBar = T, 
                     samplesAsCodaMCMC = T, WAIC = TRUE)
  
  
  return(results)
  
  
}


#### Function to fit Poisson 3 phi model ####
model_fit_3phi <- function(data) {
  
  jd_smooth <- jagam(y ~ s(lon, lat, k = 10, bs = "tp") +
                       s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
                     offset = log(E),
                     file = "jagam_smooth.txt", 
                     # Set prior for smoothing function (in this case, the lon/lat spatial smooth)
                     sp.prior = "gamma", 
                     family = poisson,
                     diagonalize = F, 
                     data = data)
  
  
  # Extract basis functions to use as linear predictors in the model
  X_dist <- jd_smooth$jags.data$X[,2:10]
  X_human <- jd_smooth$jags.data$X[,11:19]
  
  # Set constants for model (number obs +  number coefficients)
  Consts <- list(n = length(data$y), m_dist = ncol(X_dist),
                 m_human = ncol(X_human), 
                 m = (ncol(X_dist) + ncol(X_human)))
  
  
  ## Write model formula
  Model <- nimbleCode({ 
    
    # u_dist = spatial smooth term based on distance
    u_dist[1:n] <- X_dist[1:n, 1:m_dist] %*% b[2:(m_dist + 1)] 
    
    # u_human = spatial smooth term based on human movement
    u_human[1:n] <- X_human[1:n, 1:m_human] %*% b[(m_dist + 2):(m + 1)] 
    
    for (i in 1:n) { 
      # y = number of cases
      y[i] ~ dpois(mu[i]) 
      
      log(mu[i]) <- b[1] + u_dist[i] + u_human[i] + v[i] +log(e[i])
      
      # v = iid random effect
      v[i] ~ dnorm(0, sd = sig_re)
    } 
    
    # Priors
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
      lambda[i] ~ T(dgamma(.05,.005), 0, 10)
    }
  } )
  
  
  # Convert jagam data into data suitable for nimble
  nimbleData <- list(y = data$y, X_dist = X_dist, X_human = X_human,
                     e = data$E, zero = jd_smooth$jags.data$zero,
                     S1 = jd_smooth$jags.data$S1,
                     S2 = jd_smooth$jags.data$S2)
  
  
  # Set initial values for MCMC
  inits <- list(b = rnorm((ncol(X_dist) + ncol(X_human) + 1), sd=0.1), 
                v =  rnorm(nrow(data), 1),
                lambda = rep(3, 4), sig_re = .5)
  
  # Sets up model in nimble code
  nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                             constants = Consts, data = nimbleData, 
                             inits = inits)
  
  # Tell model which parameter to estimate and return
  MCMCconfig <- configureMCMC(nimbleModel,monitors=c("b","lambda",
                                                     "u_dist", "u_human",
                                                     "v", "sig_re", "mu"))
  
  
  modelMCMC <- buildMCMC(MCMCconfig)
  compiled_model <- compileNimble(nimbleModel)
  compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
  results <- runMCMC(compiled_model_MCMC, niter = 50000, nburnin = 10000, 
                     nchains = 3, thin=10, inits=inits, progressBar = T, 
                     samplesAsCodaMCMC = T)
  
  
  return(results)
  
}

#### Function to fit binomial spatial smooth + iid model ####
smooth_binom_fit <- function(data) {
  
  # Use mgcv to extract basis functions for the smooth term
  jd_bin <- jagam(y/n_yrs ~  s(lon, lat, k = 10, bs = "tp"),
                  weights = data$n_yrs,
                  # Use same formula as if running a GAM
                  file = "jagam_bin.txt", 
                  # Set prior for smoothing function (in this case, the lon/lat spatial smooth)
                  sp.prior = "gamma", 
                  family = binomial,
                  # If T, smooths are re-parameterised to have an iid Gaussian prior (not appropriate here)
                  # diagonalize cannot be used with >1 dimension (e.g tensor) or when forecasting/extrapolating beyond observed points
                  diagonalize = F, 
                  data = data)
  
  # Extract basis functions to use as linear predictors in the model
  X <- jd_bin$jags.data$X
  
  # Set constants for model (number obs +  number coefficients)
  Consts <- list(n = nrow(data), m = ncol(X))
  
  ## Write model formula
  Model <- nimbleCode({ 
    
    # u = spatial smooth term
    u[1:n] <- X[1:n, 2:m] %*% b[2:m] 
    
    for (i in 1:n) { 
      # y = number of cases
      y[i] ~ dbin(mu[i], 20) 
      
      logit(mu[i]) <- b[1] + u[i] + v[i]
      
      # v = iid random effect
      v[i] ~ dnorm(0, sd = sig_re)
    } 
    
    # Priors
    # Random effect SD
    sig_re ~ dexp(.1)
    
    # Intercept
    b[1] ~ dnorm(0, sd = 5) 
    
    ## prior for sd(s(lon,lat))
    K1[1:(m-1),1:(m-1)] <- S1[1:(m-1),1:(m-1)] * lambda[1] + 
      S1[1:(m-1),m:(2*(m-1))] * lambda[2]
    
    # Prior for smooth coefficient
    b[2:m] ~ dmnorm(zero[2:m], K1[1:(m-1),1:(m-1)]) 
    
    ## smoothing parameter priors 
    for (i in 1:2) {
      lambda[i] ~ T(dgamma(.05,.005), 0, 10)
    }
  } )
  
  
  # Convert jagam data into data suitable for nimble
  nimbleData <- list(y = data$y, X = X, zero = jd_bin$jags.data$zero, 
                     S1 = jd_bin$jags.data$S1)
  
  # Set initial values for MCMC
  inits <- list(b = rnorm(ncol(X), sd = 1), lambda=c(3, 3), 
                sig_re = runif(1), v = rnorm(nrow(data), 1))
  
  # Sets up model in nimble code
  nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                             constants = Consts, data = nimbleData, 
                             inits = inits)
  
  # Tell model which parameter to estimate and return
  MCMCconfig <- configureMCMC(nimbleModel,
                              monitors=c("b","lambda", "u", "v", 
                                         "sig_re", "mu", "y"),
                              enableWAIC = TRUE)
  
  
  modelMCMC <- buildMCMC(MCMCconfig)
  compiled_model <- compileNimble(nimbleModel)
  compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
  results <- runMCMC(compiled_model_MCMC, thin = 10, 
                     niter = 100000, nburnin = 50000, 
                     nchains = 3, inits=inits, progressBar = T, 
                     samplesAsCodaMCMC = T, WAIC = TRUE)
  
}


#### Function to fit binomial INLA model with BYM2 structure ####
inla_binom_fit <- function(df) {
  
  # Set prior using PC priors
  phi.prior = list(prior = "pc", 
                   param = c(0.5, 2/3)) # P(phi < 0.5) = 2/3 
  prec.prior = list(prior = "pc.prec",
                    param = c(1, 0.01)) # P(st. dev > 1) = 0.01
  
  
  inla_model <- inla(y ~ f(index, model = "bym2", 
                           graph = "output/sim_studys2/inla_nb.graph",
                           hyper = list(prec = prec.prior,
                                        phi = phi.prior)),
                     data = df, family = "binomial", Ntrials = 20,
                     control.compute = list(config = T, waic = T),
                     inla.mode = "experimental")
  
  return(inla_model)
}


#### Function to fit binomial spatial smooth model with 3 random effects ####
binom_3phi_fit <- function(data) {
  
  jd_smooth <- jagam(y/n_yrs ~ s(lon, lat, k = 10, bs = "tp") +
                       s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
                     weights = data$n_yrs,
                     file = "jagam_cs2.txt", 
                     # Set prior for smoothing function (in this case, the lon/lat spatial smooth)
                     sp.prior = "gamma", 
                     family = binomial,
                     diagonalize = F, 
                     data = data)
  
  # Extract basis functions to use as linear predictors in the model
  X_dist <- jd_smooth$jags.data$X[,2:10]
  X_human <- jd_smooth$jags.data$X[,11:19]
  
  # Set constants for model (number obs +  number coefficients)
  Consts <- list(n = nrow(data), n_yrs = data$n_yrs[1], 
                 m_dist = ncol(X_dist), 
                 m_human = ncol(X_human), 
                 m_tot = ncol(jd_smooth$jags.data$X))
  
  
  ## Write model formula
  Model <- nimbleCode({ 
    
    # u_dist = spatial smooth term based on distance
    u_dist[1:n] <- X_dist[1:n, ] %*% b[2:(m_tot - m_human)] 
    
    # u_human = spatial smooth term based on human movement
    u_human[1:n] <- X_human[1:n, ] %*% b[(m_tot - m_human + 1):(m_tot)] 
    
    for (i in 1:n) { 
      # y = number of cases
      y[i] ~ dbin(p[i], n_yrs) 
      
      logit(p[i]) <- b[1] + u_dist[i] + u_human[i] + v[i]
      
      # v = iid random effect
      v[i] ~ dnorm(0, sd = sig_re)
      
    } 
    
    # Priors
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
    b[2:(m_tot - m_human)] ~ dmnorm(zero[1:m_dist], 
                                    K1[1:m_dist, 1:m_dist]) 
    
    b[(m_tot - m_human + 1):(m_tot)] ~ dmnorm(zero[1:m_human], 
                                              K2[1:m_human, 1:m_human]) 
    
    ## smoothing parameter priors 
    for (i in 1:4) {
      
      lambda[i] ~ T(dgamma(.05,.005), 0, 10)
      
    }
    
  } )
  
  
  # Convert jagam data into data suitable for nimble
  nimbleData <- list(y = data$y, 
                     X_dist = X_dist, 
                     X_human = X_human,
                     zero = jd_smooth$jags.data$zero,
                     S1 = jd_smooth$jags.data$S1,
                     S2 = jd_smooth$jags.data$S2)
  
  
  
  # Set initial values for MCMC
  inits <- list(b = rnorm(ncol(jd_smooth$jags.data$X), sd=0.1), 
                v =  rnorm(nrow(data), 1), 
                lambda = rep(1, 4), sig_re = .5)
  
  # Sets up model in nimble code
  nimbleModel <- nimbleModel(code = Model, name = 'nimbleModel', 
                             constants = Consts, data = nimbleData, 
                             inits = inits)
  
  # Tell model which parameter to estimate and return
  MCMCconfig <- configureMCMC(nimbleModel,monitors=c("b", "lambda", "u_dist", 
                                                     "u_human", "v", 
                                                     "sig_re", "p", "y"))
  
  
  modelMCMC <- buildMCMC(MCMCconfig)
  compiled_model <- compileNimble(nimbleModel)
  compiled_model_MCMC <- compileNimble(modelMCMC, project = nimbleModel)
  results <- runMCMC(compiled_model_MCMC, niter = 50000, nburnin = 10000, 
                     nchains = 3, thin=10, inits=inits, progressBar = T, 
                     samplesAsCodaMCMC = T)
  
}



#### Functions to extract model output ####
#### Function to extract intercept term from smooth models ####
b_extract_smooth <- function(results) {
  
  ## Extract simulations of the intercept
  b_sim <- do.call(rbind, results$samples)[ , "b[1]"]
  
  # Return mean and 95% credible interval & format 
  b_est <- data.table(b_est = mean(b_sim), 
                      b_lq = quantile(b_sim, .025),
                      b_uq = quantile(b_sim, .975),
                      b_format = paste0(round(mean(b_sim), 3), " (",
                                        round(quantile(b_sim, .025), 3), 
                                        ", ", round(quantile(b_sim, .975), 3), ")"))
  
  return(b_est)
  
}


#### Function to extract intercept term from BYM2 NIMBLE model ####
b_extract_bym2 <- function(results) {
  
  ## Extract simulations of the intercept
  b_sim <- do.call(rbind, results$samples)[ , "b"]
  
  # Return mean and 95% credible interval & format 
  b_est <- data.table(b_est = mean(b_sim), 
                      b_lq = quantile(b_sim, .025),
                      b_uq = quantile(b_sim, .975),
                      b_format = paste0(round(mean(b_sim), 3), " (",
                                        round(quantile(b_sim, .025), 3), 
                                        ", ", round(quantile(b_sim, .975), 3), ")"))
  
  return(b_est)
  
}


#### Function to extract intercept term from 3 phi models ####
b_extract_3phi <- function(results) {
  
  ## Extract simulations of the intercept
  b_sim <- do.call(rbind, results)[ , "b[1]"]
  
  # Return mean and 95% credible interval & format 
  b_est <- data.table(b_est = mean(b_sim), 
                      b_lq = quantile(b_sim, .025),
                      b_uq = quantile(b_sim, .975),
                      b_format = paste0(round(mean(b_sim), 3), " (",
                                        round(quantile(b_sim, .025), 3), 
                                        ", ", round(quantile(b_sim, .975), 3), ")"))
  
  return(b_est)
  
}


#### Function to extract intercept term from INLA models ####
b_extract_inla <- function (results) {
  
  # Create table with mean and 95% credible interval of intercept estimates
  b_est <- data.table(b_est_inla = results$summary.fixed[1],
                      b_lq_inla = results$summary.fixed[3],
                      b_uq_inla = results$summary.fixed[5],
                      b_est_inla_format = 
                        paste0(round(results$summary.fixed[1], 3),
                               " (", round(results$summary.fixed[3], 3),
                               ", ", round(results$summary.fixed[5], 3),
                               ")"))
  
}


#### Function to extract variance of random effects from smooth model (1 spatial structure) ####
extract_vars_2re <- function(results) {
  
  n <- nrow(shp)
  
  # Return column numbers with structured random effect simulations
  umin <- which(colnames(results$samples[[1]]) == "u[1]")
  umax <- which(colnames(results$samples[[1]]) == paste0("u[", n, "]"))
  
  # Return column numbers with unstructured random effect simulations
  vmin <- which(colnames(results$samples[[1]]) == "v[1]")
  vmax <- which(colnames(results$samples[[1]]) == paste0("v[", n, "]"))
  
  
  # Extract simulations of structured random effect
  u_mat <- do.call(rbind, results$samples)[, umin:umax]
  # Estimate the variance of each simulation
  u_var <- apply(u_mat, 1, var)
  
  # Extract simulations of unstructured random effect
  v_mat <- do.call(rbind, results$samples)[,vmin:vmax]
  # Estimate the variance of each simulation
  v_var <- apply(v_mat, 1, var)
  
  # Calculate the proportion of variance explained by the structured term
  propn_spat <- u_var/(u_var + v_var)
  
  re_vars <- data.table(spat_var = u_var, 
                        iid_var = v_var,
                        propn_var = propn_spat)
  
  return(re_vars)
  
} 


#### Function to extract phi parameter from BYM2 model (NIMBLE) ####
phi_extract_bym2 <- function(results) {
  
  ## Extract simulations of the intercept
  phi_sim <- do.call(rbind, results$samples)[ , "phi"]
  
  # Return mean and 95% credible interval & format 
  phi_est <- data.table(phi_est = mean(phi_sim), 
                        phi_lq = quantile(phi_sim, .025),
                        phi_uq = quantile(phi_sim, .975),
                        phi_format = paste0(round(mean(phi_sim), 3), " (",
                                            round(quantile(phi_sim, .025), 3), 
                                            ", ", round(quantile(phi_sim, .975), 3), ")"))
  
  return(phi_est)
  
}


#### Function to extract variance of random effects from smooth model (2 spatial structure) ####
extract_vars_3re <- function(results) {
  
  n <- nrow(shp)
  
  # Return column numbers with distance-based random effect simulations
  u_dist_min <- which(colnames(results[[1]]) == "u_dist[1]")
  u_dist_max <- which(colnames(results[[1]]) == paste0("u_dist[", n, "]"))
  
  # Return column numbers with human movement-based random effect simulations
  u_human_min <- which(colnames(results[[1]]) == "u_human[1]")
  u_human_max <- which(colnames(results[[1]]) == paste0("u_human[", n, "]"))
  
  # Return column numbers with unstructured random effect simulations
  vmin <- which(colnames(results[[1]]) == "v[1]")
  vmax <- which(colnames(results[[1]]) == paste0("v[", n, "]"))
  

  # Extract simulations of distance-based structured random effect
  u_dist_mat <- do.call(rbind, results)[,u_dist_min:u_dist_max]
  # Estimate the variance of each simulation
  u_dist_var <- apply(u_dist_mat, 1, var)
  
  # Extract simulations of human movement-based structured random effect
  u_human_mat <- do.call(rbind, results)[,u_human_min:u_human_max]
  # Estimate the variance of each simulation
  u_human_var <- apply(u_human_mat, 1, var)
  
  # Extract simulations of unstructured random effect
  v_mat <- do.call(rbind, results)[,vmin:vmax]
  # Estimate the variance of each simulation
  v_var <- apply(v_mat, 1, var)
  
  # Calculate the proportion of variance explained by each term
  propn_spat_dist <- u_dist_var/(u_dist_var +  u_human_var + v_var)
  propn_spat_human <- u_human_var/(u_dist_var +  u_human_var + v_var)
  propn_spat_iid <- v_var/(u_dist_var +  u_human_var + v_var)
  
  re_vars <- data.table(spat_dist_var = u_dist_var,
                        spat_human_var = u_human_var,
                        iid_var = v_var,
                        propn_var_dist = propn_spat_dist,
                        propn_var_human = propn_spat_human,
                        propn_var_iid = propn_spat_iid)
  
  return(re_vars)
  
} 


#### Function to extract random effect estimates (2 spatial structures) ####
re_est_3phi <- function(model) {

  n <- nrow(shp)
  
  u_dist_min <- which(colnames(model[[1]]) == "u_dist[1]")
  u_dist_max <- which(colnames(model[[1]]) == paste0("u_dist[", n, "]"))

  u_human_min <- which(colnames(model[[1]]) == "u_human[1]")
  u_human_max <- which(colnames(model[[1]]) == paste0("u_human[", n, "]"))

  vmin <- which(colnames(model[[1]]) == "v[1]")
  vmax <- which(colnames(model[[1]]) == paste0("v[", n, "]"))

  ## Calculate the mean of each random effect term
  u_dist_mat <- do.call(rbind, model)[ , u_dist_min:u_dist_max]
  u_dist_mean <- apply(u_dist_mat, 2, mean)

  u_human_mat <- do.call(rbind, model)[ , u_human_min:u_human_max]
  u_human_mean <- apply(u_human_mat, 2, mean)

  v_mat <- do.call(rbind, model)[ , vmin:vmax]
  v_mean <- apply(v_mat, 2, mean)

  # Calculate the total random effect term
  re_mean <- apply((u_dist_mat + u_human_mat + v_mat), 2, mean)

  df_re <- shp %>% 
    mutate(tot_re = re_mean,
           spat_dist_re = u_dist_mean, 
           spat_human_re = u_human_mean,
           iid_re = v_mean)

  return(df_re)
}


#### Function to extract outcome estimate from smooth model (1 spatial structure) ####
lambda_pred_2re <- function(results) {
  
  n <- nrow(shp)
  
  # Return column numbers with lambda simulations
  lammin <- which(colnames(results$samples[[1]]) == "mu[1]")
  lammax <- which(colnames(results$samples[[1]]) == paste0("mu[", n, "]"))
  
  ## # Extract simulations of lambda
  lam_mat <- do.call(rbind, results$samples)[ , lammin:lammax]
  # Estimate mean from each simulation
  lam_mean <- apply(lam_mat, 2, mean)
  
  return(lam_mean)
  
}


#### Function to draw random samples from multivariate normal distribution ####
# n = number simulations
# mu = matrix of means
# sig = matrix of standard  deviations
rmvn <- function(n, mu, sig) { 
  L <- mroot(sig); m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}


#### Function to predict the probability from a binary GAM model ####
pred_binomial_prob <- function(model, data = df_gam) {
  
  n.sims <- 1000
  
  # Draw simulations from the coefficient (beta) posterior distribution
  betas <- rmvn(n.sims, coef(model), model$Vp) 
  
  # Extract matrix of linear predictors from model
  model_mat <- predict(object = model, newdata = data, type = "lpmatrix")
  
  # Use linear prediction matrix and coefficient simulations to estimate posterior mean
  # Convert response to estimate the probability of an outbreak
  mean_prob <- exp(model_mat %*% t(betas))/(1 + exp(model_mat %*% t(betas)))
  
  ## Obtain probability prediction (number of 1s/total simulations)
  mean_prob_pred <- apply(mean_prob, 1, mean)
  
  # Combine information into a df
  predictions <- data.table(source_code = data$source_code,
                            dest_code = data$dest_code,
                            prob_connect = mean_prob_pred,
                            connected = data$connected)
  
  return(predictions)
  
}


