#########################################################################
####                                                                 ####
####  Simulation study S1: a human movement-based spatial structure  ####
####                       Script to fit models                      ####
####                                                                 ####
#########################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudys1_simulation.R has been run
df_sim_human <- read_rds("data/df_sim_human.rds")


#### Step 1: Fit the spatial smooth model ####
## Fit GAM to obtain sensible priors for lambdas
gam_fit <- function(df) {
  model <- gam(y ~  s(connect_coord1, connect_coord2, k = 10, bs = "tp"), 
               offset = log(E), family = "poisson", 
               data = df, method = "REML")
  
  return(model) 
}

sp_check <- lapply(df_sim_human, gam_fit)

# Extract smoothing parameter, sp, from GAMs
lapply(sp_check, function(x) x$sp) 


## Fit a model with a smooth spatial random effect and an IID random
## effect in NIMBLE to data with different (known) contributions 
# Rename connectivity coordinates to lon/lat for function
for(i in 1:length(df_sim_human)) {
  
  df_sim_human[[i]] <- df_sim_human[[i]] %>% 
    dplyr::select(-lon, -lat) %>% 
    rename(lon = connect_coord1,
           lat = connect_coord2)

}

model_list <- lapply(df_sim_human, smooth_iid_fit)


## Check convergence using Gelman-Rubin diagnostics with coda 
GR.diag <- lapply(model_list, 
                  function(x) gelman.diag(x$samples, multivariate = F))


lapply(GR.diag, function(x) sum(x$psrf[, "Point est."] > 1.1) )


# Save the model object
write_rds(model_list, 
          file = "output/sim_studys1/model_list.rds")
