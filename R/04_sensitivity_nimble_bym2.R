###################################################################
####                                                           ####
####  Sensitivity study 1: a distance-based spatial structure  ####
####                NIMBLE CAR model comparison                ####
####                                                           ####
###################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data ####
## NOTE: this code will only work if 01_simstudy1_simulation.R has been run
df_sim_lonlat <- read_rds("data/df_sim_lonlat.rds")


#### Load spatial smooth model ####
## NOTE: this code will only work if 01_simstudy1_simulation.R and
## 02_simstudy1_run_models.R has been run
model_list <- read_rds("output/sim_study1/smooth_model_list.rds")

#### Load INLA BYM2 model ####
inla_model_list <- read_rds("output/sim_study1/inla_model_list.rds")


#### Fit BYM2 model using NIMBLE ####
## Create nb object from shapefile 
row.names(shp) = shp$municip_code_ibge

nb.south <- poly2nb(shp)

# Convert to WinBUGS object
nb.WB <- nb2WB(nb.south)

# Calculate scale parameter (needed for BYM2 scaling and model fit)
nb.mat <- nb2mat(nb.south, style = "B")
colnames(nb.mat) <- rownames(nb.mat)

nb.scale <- -nb.mat
diag(nb.scale) <- abs(apply(nb.scale, 1, sum))
# solve(W.scale) # this should not work since by definition the matrix is singular

Q = inla.scale.model(nb.scale, 
                     constr = list(A = matrix(1, nrow = 1, 
                                               ncol = nrow(nb.scale)), e = 0))
scale = exp((1/nrow(nb.scale)) * sum(log(1/diag(Q))))


#### Fit BYM2 model with NIMBLE for comparison ####
bym2_results <- lapply(df_sim_lonlat, BYM2_nimble_fit)

## Save model results
# write_rds(bym2_results, file = "output/simstudy1/bym2_sens_models.rds")


#### Extract intercept estimates from simulations ####
b_est_smooth <- lapply(model_list, b_extract_smooth)
b_est_bym2 <- lapply(bym2_results, b_extract_bym2)
b_est_inla <- lapply(inla_model_list, b_extract_inla)



# Combine estimates and add known phi values
b_table_smooth <- reduce(b_est_smooth, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))

b_table_bym2 <- reduce(b_est_bym2, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))

b_table_inla <- reduce(b_est_inla, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))


## Combine intercept estimates into one dataset
b_est_full <- full_join(b_table_smooth, b_table_bym2, by = "phi",
                        suffix = c(".smooth", ".bym2")) %>% 
  full_join(., b_table_inla, by = "phi")


#### Plot intercept estimates for smooth & BYM2 models ####
b_line_full <- ggplot(data = b_est_full) +
  # Plot smooth estimates + 95% CI
  geom_point(aes(x = phi, y = b_est.smooth)) +
  geom_linerange(aes(ymin = b_lq.smooth, ymax = b_uq.smooth, x = phi), lwd = 1) +
  # Add INLA estimates + 95% CI
  geom_point(aes(x = (phi + .01), y = b_est_inla.mean), col = "blue") +
  geom_linerange(aes(ymin = b_lq_inla.0.025quant, 
                     ymax = b_uq_inla.0.975quant, x = (phi + .01)), lwd = 1,
                 col = "blue") +
  # Plot BYM2 NIMBLE estimates + 95% CI
  geom_point(aes(x = (phi - .01), y = b_est.bym2), col = "magenta") +
  geom_linerange(aes(ymin = b_lq.bym2, 
                     ymax = b_uq.bym2, 
                     # Shift to avoid overlap
                     x = (phi - .01)), lwd = 1, col = "magenta") +
  # Add reference line with true value (b = 0)
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(name = expression(phi), breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(name = "Intercept estimate") +
  theme_bw()

ggsave(b_line_full, filename = "output/sim_study1/bym2_nimble_b_comp_full.png")


#### Extract estimates for phi/mixing parameters ####
n <- nrow(shp)

## Spatial smooth model
smooth_model_vars <- lapply(model_list, extract_vars_2re)

## BYM2 nimble model
bym2_model_vars <- lapply(bym2_results, phi_extract_bym2)

## INLA model 
phi_inla <- NA
phi_inla_lq <- NA
phi_inla_uq <- NA

for(i in 1: length(inla_model_list)) {
  phi_inla[i] <- inla_model_list[[i]]$summary.hyperpar[2,1]
  phi_inla_lq[i] <- inla_model_list[[i]]$summary.hyperpar[2, 3]
  phi_inla_uq[i] <- inla_model_list[[i]]$summary.hyperpar[2, 5]
}

inla_phi_est <- data.table(true_phi = seq(0, 1, by = .1),
                           inla_phi = phi_inla,
                           inla_phi_lq = phi_inla_lq,
                           inla_phi_uq = phi_inla_uq)


## Combine the estimated phi from each model and calculate mean and 95% CI
smooth_phi_est <- reduce(smooth_model_vars, rbind) %>% 
  #  Add true phi value
  mutate(phi = rep(seq(0,  1, by = .1), 
                   each = nrow(smooth_model_vars[[1]]))) %>% 
  group_by(phi) %>% 
  summarise(phi_est = mean(propn_var),
            phi_lq = quantile(propn_var, .025),
            phi_uq = quantile(propn_var, .975)) %>%
  ungroup()


bym2_phi_est <- reduce(bym2_model_vars, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))



#### Plot mean & 95% CI phi estimates for all 3 models ####
phi_comp_plot_full <- ggplot() +
  geom_point(data = smooth_phi_est,
             aes(x =  phi, y = phi_est)) +
  geom_linerange(data = smooth_phi_est,
                 aes(ymin = phi_lq, ymax = phi_uq, x = phi), lwd = 1) +
  geom_point(data = bym2_phi_est,
             # Shift to avoid overlap
             aes(x = (phi - .01), y = phi_est), col = "magenta") +
  geom_linerange(data = bym2_phi_est,
                 aes(ymin = phi_lq, ymax = phi_uq, 
                     x = (phi - .01)), lwd = 1, col = "magenta") +
  geom_point(data = inla_phi_est,
             # Shift to avoid overlap
             aes(x = (true_phi + .01), y = inla_phi), col = "blue") +
  geom_linerange(data = inla_phi_est,
                 aes(ymin = inla_phi_lq, ymax = inla_phi_uq, 
                     x = (true_phi + .01)), lwd = 1, col = "blue") +
  # Add reference line with simulation value
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = expression("True" ~ phi), y = expression("Estimated"~ phi)) +
  expand_limits(x = c(0, 1), y = c(0, 1)) +
  theme_bw()

ggsave(phi_comp_plot_full, filename = "output/sim_study1/bym2_nimble_phi_comp_full.png")



## Mean absolute error
## BYM2 NIMBLE
# Extract predicted values 
n <- nrow(df_sim_lonlat[[1]])

mu_est <- lapply(bym2_results, lambda_pred_2re)

# Calculate the absolute error for each model
for(i in 1:length(df_sim_lonlat)) {
  df_sim_lonlat[[i]] <- df_sim_lonlat[[i]] %>% 
    mutate(lambda_est = (mu_est[[i]]/E), 
           absolute_error = abs(lambda - lambda_est))
}

# Calculate mean absolute error
mae_bym2<- lapply(df_sim_lonlat, function(x) mean(x$absolute_error))

# Extract WAIC 
waic_bym2 <- lapply(bym2_results, function(x) x$WAIC$WAIC)


## Import other comparison statistics
## NOTE: this code will only work if 01_simstudy1_simulation.R, 
## 02_simstudy1_run_models.R and 03_simstudy1_model_output have been run
inla_smooth_comp <- fread(file = "output/sim_study1/model_comp.csv")



# Combine stats for smooth INLA and BYM2 models
model_comp <- data.table(inla_smooth_comp,
                         bym2_mae = round(unlist(mae_bym2), 2),
                         bym2_waic = round(unlist(waic_bym2), 2),
                         bym2_phi = round(bym2_phi_est$phi_est, 3))

fwrite(model_comp, file = "output/sim_study1/bym2_nimble_comp_full.csv")







