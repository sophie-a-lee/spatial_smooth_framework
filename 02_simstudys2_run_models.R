############################################################################
####                                                                    ####
####  Simulation study S2: a binomial distance-based spatial structure  ####
####                         Script to fit models                       ####
####                                                                    ####
############################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudys2_simulation.R has been run
df_binom <- read_rds("data/df_sim_binomial_2s.rds")


#### Step 1: Fit the spatial smooth model ####
## Fit GAM to obtain sensible priors for lambdas
gam_fit <- function(df) {
  model <- gam(cbind(y, n_yrs - y) ~ s(lon, lat, k = 10, bs = "tp"),
               family = binomial("logit"), data = df, method = "REML")
  
  return(model) 
}

sp_check <- lapply(df_binom, gam_fit)

# Extract smoothing parameter, sp, from GAMs
lapply(sp_check, function(x) x$sp) 


## Fit a model with a smooth spatial random effect and an IID random
## effect in NIMBLE to data with different (known) contributions 
smooth_model_list <- lapply(df_binom, smooth_binom_fit)


# Save the model object
write_rds(smooth_model_list, 
          file = "output/sim_studys2/smooth_model_list.rds")


####  Step 2: Fit a binomial INLA model to the simulated data ####
## Create a binary neighbourhood matrix
nb_br <- poly2nb(shp, row.names = shp$index)
names(nb_br) <- attr(nb_br, "region.id")


# Convert into an INLA.nb object
nb2INLA("output/sim_studys2/inla_nb.graph", nb_br)


## Use a BYM2 random effect which has a spatial structure (assuming cities are connected 
# if and only if they share a border) and an unstructured element
inla_model_list <- lapply(df_binom, inla_binom_fit)

write_rds(inla_model_list, 
          file = "output/sim_studys2/inla_model_list.rds")




