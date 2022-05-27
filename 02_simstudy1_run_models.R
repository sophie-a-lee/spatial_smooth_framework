##################################################################
####                                                          ####
####  Simulation study 1: a distance-based spatial structure  ####
####                   Script to fit models                   ####
####                                                          ####
##################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudy1_simulation.R has been run
df_sim_lonlat <- read_rds("data/df_sim_lonlat.rds")


#### Step 1: Fit the spatial smooth model ####
## Fit GAM to obtain sensible priors for lambdas
gam_fit <- function(df) {
  model <- gam(y ~  s(lon, lat, k = 10, bs = "tp"), offset = log(E),
               family = "poisson", data = df, method = "REML")
  
  return(model) 
}

sp_check <- lapply(df_sim_lonlat, gam_fit)

# Extract smoothing parameter, sp, from GAMs
lapply(sp_check, function(x) x$sp) 


## Fit a model with a smooth spatial random effect and an IID random
## effect in NIMBLE to data with different (known) contributions 
smooth_model_list <- lapply(df_sim_lonlat, smooth_iid_fit)


## Check convergence using Gelman-Rubin diagnostics with coda 
GR.diag <- lapply(smooth_model_list, 
                  function(x) gelman.diag(x$samples, multivariate = F))


lapply(GR.diag, function(x) sum(x$psrf[, "Point est."] > 1.1) )


# Save the model object
write_rds(smooth_model_list, 
          file = "output/sim_study1/smooth_model_list.rds")



####  Step 2: Fit an INLA model to the simulated data ####
## Create a binary neighbourhood matrix
nb_br <- poly2nb(shp, row.names = shp$index)
names(nb_br) <- attr(nb_br, "region.id")


# Convert into an INLA.nb object
nb2INLA("output/sim_study1/inla_nb.graph", nb_br)


## Use a BYM2 random effect which has a spatial structure (assuming cities are connected 
# if and only if they share a border) and an unstructured element
inla_model_list <- lapply(df_sim_lonlat, inla_fit)

write_rds(inla_model_list, 
          file = "output/sim_study1/inla_model_list.rds")






