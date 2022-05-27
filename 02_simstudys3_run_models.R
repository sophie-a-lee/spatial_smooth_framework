############################################################################
####                                                                    ####
####  Simulation study S3: a binomial model with 2 spatial structures   ####
####                         Script to fit models                       ####
####                                                                    ####
############################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudys2_simulation.R has been run
df_binom <- read_rds("data/df_sim_binomial_s3.rds")


#### Fit the spatial smooth model ####
## Fit GAM to obtain sensible priors for lambdas
gam_fit <- function(df) {
  model <- gam(cbind(y, n_yrs - y) ~ s(lon, lat, k = 10, bs = "tp") +
                 s(connect_coord1, connect_coord2, k = 10, bs = "tp"),
               family = binomial("logit"), data = df, method = "REML")
  
  return(model) 
}

sp_check <- lapply(df_binom, gam_fit)

# Extract smoothing parameter, sp, from GAMs
lapply(sp_check, function(x) x$sp) 


## Fit a model with 3 random effects in NIMBLE to data with different (known) contributions
binom_model_list  <- lapply(df_binom, binom_3phi_fit)


# Save the model object
write_rds(binom_model_list, 
          file = "output/sim_studys3/model_list.rds")





