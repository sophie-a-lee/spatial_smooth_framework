############################################################################
####                                                                    ####
####  Simulation study S2: a binomial distance-based spatial structure  ####
####         Script to explore model ouput and compare models           ####
####                                                                    ####
############################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data and models ####
## NOTE: this code will only work if 01_simstudys2_simulation.R and
## 02_simstudys2_run_models.R has been run
df_sim <- read_rds("data/df_sim_binary_s2.rds")
df_binom <- read_rds("data/df_sim_binomial_2s.rds")

smooth_model_list <- read_rds("output/sim_studys2/smooth_model_list.rds")

inla_model_list <-  read_rds("output/sim_studys2/inla_model_list.rds")


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
ggsave(b_line, filename = "output/sim_studys2/intercept_line.png")


# Save estimates in a table
b_table <- b_est_full %>% 
  dplyr::select(phi, b_format, b_est_inla_format)

fwrite(b_table, "output/sim_studys2/intercept_est.csv")


#### Extract estimates for phi/mixing parameters ####
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
       filename = "output/sim_studys2/phi_comp_plot.png")


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


## Brier score
# Smooth model
n <-nrow(shp)

# Extract fitted outcome from spatial smooth model
mu_list <- lapply(smooth_model_list, lambda_pred_2re)

# Estimate the absolute square difference between fitted and observed outcome
for(i in 1:length(df_sim)) {
  df_sim[[i]] <- df_sim[[i]] %>% 
    mutate(p_est_smooth = as.vector(rep(mu_list[[i]], times = 20)),
           sqdiff_smooth = (y - p_est_smooth)^2,
           p_est_inla = rep(exp(inla_model_list[[i]]$summary.fitted.values$`0.5quant`) /
                              (1 + exp(inla_model_list[[i]]$summary.fitted.values$`0.5quant`)), 
                            times = 20),
           sqdiff_inla = abs(y - p_est_inla)^2)
}


brier_smooth <- lapply(df_sim, function(x) mean(x$sqdiff_smooth))

brier_inla <- lapply(df_sim, function(x) mean(x$sqdiff_inla))


## ROC curve
# Calculate AUC + create ROC objects
roc_smooth <- lapply(df_sim, function(x) roc(x$y, x$p_est_smooth,
                                             auc = T, ci = T))

roc_inla <- lapply(df_sim, function(x) roc(x$y, x$p_est_inla,
                                           auc = T, ci = T))


# Extract sensitivities + specificities to plot
roc_obj <- list(NA) 

for(i in 1:length(df_sim)) {
  phi <- gsub("\\.", "", df_sim[[i]]$phi_true[1])
  
  roc_obj[[i]] <- data.table(roc_smooth_sens = roc_smooth[[i]]$sensitivities,
                             roc_smooth_spec = 1 - roc_smooth[[i]]$specificities,
                             roc_inla_sens = roc_inla[[i]]$sensitivities,
                             roc_inla_spec = 1 - roc_inla[[i]]$specificities)
  
  # Plot ROC curve for each model
  roc_curve <- ggplot(data = roc_obj[[i]]) +
    # Add a line per threshold
    geom_line(aes(y = roc_smooth_sens,  x = roc_smooth_spec)) +
    geom_line(aes(y = roc_inla_sens,  x = roc_inla_spec), linetype = "dashed",
              col = "blue") +
    # Add reference line (y = x line represents chance)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "True negative rate", y = "True positive rate") +
    theme_light()
  
  ggsave(roc_curve, 
         filename = paste0("output/sim_studys2/roc_curve_", phi, ".png"))
}


## Save comparisons in table
model_comp <- data.table(phi = seq(0, 1, .1),
                         brier_smooth  = round(unlist(brier_smooth), 3),
                         brier_inla = round(unlist(brier_inla), 3), 
                         waic_smooth = round(unlist(waic_smooth), 1),
                         waic_inla = round(unlist(waic_inla), 1),
                         roc_smooth = round(unlist(lapply(roc_smooth, 
                                                          function(x) x$auc)), 3),
                         roc_inla = round(unlist(lapply(roc_inla, 
                                                        function(x) x$auc)), 3), 
                         phi_smooth = round(smooth_phi_est$phi_est, 3),
                         phi_inla = round(inla_phi_est$inla_phi, 3))

fwrite(model_comp, file = "output/sim_studys2/model_comp.csv")

