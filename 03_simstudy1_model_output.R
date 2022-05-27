##################################################################
####                                                          ####
####  Simulation study 1: a distance-based spatial structure  ####
####    Script to explore model ouput and compare models      ####
####                                                          ####
##################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data and models ####
## NOTE: this code will only work if 01_simstudy1_simulation.R and
## 02_simstudy1_run_models.R has been run
df_sim_lonlat <- read_rds("data/df_sim_lonlat.rds")

smooth_model_list <- read_rds("output/sim_study1/smooth_model_list.rds")

inla_model_list <-  read_rds("output/sim_study1/inla_model_list.rds")


#### Extract predicted intercept from models ####
## Spatial smooth model
b_est_smooth <- lapply(smooth_model_list, b_extract_smooth)

b_table_smooth <- reduce(b_est_smooth, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))


## INLA model
b_est_ina <- lapply(inla_model_list, b_extract_inla)


b_table_inla <- reduce(b_est_ina, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))


## Combine intercept estimates into one dataset
b_est_full <- full_join(b_table_smooth, b_table_inla, by = "phi")


#### Plot intercept estimates for smooth & INLA models ####
b_line <- ggplot(data = b_est_full) +
  # Plot smooth estimates + 95% CI
  geom_point(aes(x = phi, y = b_est)) +
  geom_linerange(aes(ymin = b_lq, ymax = b_uq, x = phi), lwd = 1) +
  # Add INLA estimates + 95% CI
  geom_point(aes(x = (phi - .01), y = b_est_inla.mean), col = "blue") +
  geom_linerange(aes(ymin = b_lq_inla.0.025quant, 
                     ymax = b_uq_inla.0.975quant, 
                     # Shift to avoid overlap
                     x = (phi - .01)), lwd = 1, col = "blue") +
  # Add reference line with true value (b = 0)
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(name = expression(phi), breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(name = "Intercept estimate") +
  theme_bw()

# Save intercept plot
ggsave(b_line, filename = "output/sim_study1/intercept_line.png")


# Save estimates in a table
b_table <- b_est_full %>% 
  dplyr::select(phi, b_format, b_est_inla_format)

fwrite(b_table, "output/sim_study1/intercept_est.csv")


#### Extract estimates for phi/mixing parameters ####
n <- nrow(shp)

## Spatial smooth model
smooth_model_vars <- lapply(smooth_model_list, extract_vars_2re)

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


## INLA model
phi_inla <- NA
phi_inla_lq <- NA
phi_inla_uq <- NA

for(i in 1: length(inla_model_list)) {
  phi_inla[i] <- summary(inla_model_list[[i]])$hyperpar[2, 1]
  phi_inla_lq[i] <- summary(inla_model_list[[i]])$hyperpar[2, 3]
  phi_inla_uq[i] <- summary(inla_model_list[[i]])$hyperpar[2, 5]
}

inla_phi_est <- data.table(true_phi = seq(0, 1, by = .1),
                           inla_phi = phi_inla,
                           inla_phi_lq = phi_inla_lq,
                           inla_phi_uq = phi_inla_uq)


#### Plot mean & 95% CI phi estimates for both models ####
phi_comp_plot <- ggplot() +
  geom_point(data = smooth_phi_est,
             aes(x =  phi, y = phi_est)) +
  geom_linerange(data = smooth_phi_est,
                 aes(ymin = phi_lq, ymax = phi_uq, x = phi), lwd = 1) +
  geom_point(data = inla_phi_est,
             # Shift to avoid overlap
             aes(x = (true_phi - .01), y = inla_phi), col = "blue") +
  geom_linerange(data = inla_phi_est,
                 aes(ymin = inla_phi_lq, ymax = inla_phi_uq, 
                     x = (true_phi - .01)), lwd = 1, col = "blue") +
  # Add reference line with simulation value
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = expression("True" ~ phi), y = expression("Estimated"~ phi)) +
  expand_limits(x = c(0, 1), y = c(0, 1)) +
  theme_bw()

ggsave(phi_comp_plot, 
       filename = "output/sim_study1/phi_comp_plot.png")


#### Calculate model comparison statistics ####
## WAIC
## Smooth model
waic_smooth <- lapply(smooth_model_list, function(x) x$WAIC$WAIC)

## INLA model
waic_inla <- lapply(inla_model_list, function(x) x$waic$waic)

# Combine WAIC for smooth and INLA models
waic_comp <- data.table(phi = seq(0, 1, by = .1), 
                        smooth = unlist(waic_smooth),
                        inla = unlist(waic_inla))


## Mean absolute error
## Smooth model
# Extract predicted values from smooth model
n <- nrow(df_sim_lonlat[[1]])

mu_est <- lapply(smooth_model_list, lambda_pred_2re)

# Calculate the absolute error for each model
for(i in 1:length(df_sim_lonlat)) {
  df_sim_lonlat[[i]] <- df_sim_lonlat[[i]] %>% 
    mutate(lambda_est = (mu_est[[i]]/E), 
           absolute_error = abs(lambda - lambda_est))
}

# Calculate mean absolute error
mae_smooth <- lapply(df_sim_lonlat, function(x) mean(x$absolute_error))


## INLA model
# Calculate the absolute error for each model
for(i in 1:length(df_sim_lonlat)) {
  df_sim_lonlat[[i]] <- df_sim_lonlat[[i]] %>% 
    mutate(lambda_est_inla = 
             exp(inla_model_list[[i]]$summary.fitted.values[,1]),
           absolute_error_inla = abs(lambda - lambda_est_inla))
}

# Calculate mean absolute error
mae_inla <- lapply(df_sim_lonlat, function(x) mean(x$absolute_error_inla))


# Combine MAE for smooth and INLA models
mae_comp <- data.table(phi = seq(0, 1, by = .1), 
                       smooth = unlist(mae_smooth), 
                       inla = unlist(mae_inla))


#### Export table with model comparisons ####
inla_smooth_comp <- full_join(mae_comp, waic_comp,
                              by = "phi",
                              suffix = c("_mae", "_waic")) %>% 
  ## Format estimates
  mutate(across(ends_with("_mae"), ~round(.x, 2)),
         across(ends_with("_waic"), ~round(.x, 2)),
         # Add phi estimates for each table
         smooth_phi = round(smooth_phi_est$phi_est, 3),
         inla_phi = round(inla_phi_est$inla_phi, 3)) %>% 
  dplyr::select(phi, smooth_mae, smooth_waic, smooth_phi,
                inla_mae, inla_waic, inla_phi)


fwrite(inla_smooth_comp, 
       file = "output/sim_study1/model_comp.csv")




