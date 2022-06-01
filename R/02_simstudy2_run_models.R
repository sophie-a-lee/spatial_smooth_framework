################################################################
####                                                        ####
####  Simulation study 2: two sources of spatial structure  ####
####                   Script to fit models                 ####
####                                                        ####
################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudy2_simulation.R has been run
df_3phi_sim <- read_rds("data/df_sim_simstudy2.rds")


#### Step 1: Fit the spatial smooth model ####
## Fit GAM to obtain sensible priors for lambdas
gam_fit <- function(df) {
  model <- gam(y ~  s(lon, lat, k = 10, bs = "tp") + 
                 s(connect_coord1, connect_coord2, k = 10, bs = "tp"), 
               offset = log(E), family = "poisson", 
               data = df, method = "REML")
  
  return(model) 
}

sp_check <- lapply(df_3phi_sim, gam_fit)

# Extract smoothing parameter, sp, from GAMs
lapply(sp_check, function(x) x$sp) 


## Fit a model with 3 random effects in NIMBLE to data with different (known) phis
model_list <- lapply(df_3phi_sim, model_fit_3phi)


## Check convergence using Gelman-Rubin diagnostics with coda 
GR.diag <- lapply(model_list, 
                  function(x) gelman.diag(x, multivariate = F))


lapply(GR.diag, function(x) sum(x$psrf[, "Point est."] > 1.1) )


# Save the model object
write_rds(model_list, 
          file = "output/sim_study2/model_list.rds")



