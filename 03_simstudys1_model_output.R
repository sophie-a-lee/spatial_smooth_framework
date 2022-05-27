#########################################################################
####                                                                 ####
####  Simulation study S1: a human movement-based spatial structure  ####
####                 Script to explore model output                  ####
####                                                                 ####
#########################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data and models ####
## NOTE: this code will only work if 01_simstudys1_simulation.R and
## 02_simstudys1_run_models.R has been run
df_sim_human <- read_rds("data/df_sim_human.rds")

model_list <- read_rds("output/sim_studys1/model_list.rds")


#### Extract predicted intercept from models ####
## Spatial smooth model
b_est <- lapply(model_list, b_extract_smooth)

b_table <- reduce(b_est, rbind) %>% 
  # Add true phi value
  mutate(phi = seq(0, 1, by = .1))


#### Plot intercept estimates for models ####
b_line <- ggplot(data = b_table) +
  # Plot estimates + 95% CI
  geom_point(aes(x = phi, y = b_est)) +
  geom_linerange(aes(ymin = b_lq, ymax = b_uq, x = phi), lwd = 1) +
  # Add reference line with true value (b = 0)
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(name = expression(phi), breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(name = "Intercept estimate") +
  theme_bw()


# Save intercept plot
ggsave(b_line, filename = "output/sim_studys1/intercept_line.png")


#### Extract estimates for phi/mixing parameters ####
## Spatial smooth model
model_vars <- lapply(model_list, extract_vars_2re)

## Combine the estimated phi from each model and calculate mean and 95% CI
phi_est <- reduce(model_vars, rbind) %>% 
  #  Add true phi value
  mutate(phi = rep(seq(0,  1, by = .1), 
                   each = nrow(model_vars[[1]]))) %>% 
  group_by(phi) %>% 
  summarise(phi_est = mean(propn_var),
            phi_lq = quantile(propn_var, .025),
            phi_uq = quantile(propn_var, .975)) %>%
  ungroup()


#### Plot mean & 95% CI phi estimates for models ####
phi_plot <- ggplot() +
  geom_point(data = phi_est,
             aes(x =  phi, y = phi_est)) +
  geom_linerange(data = phi_est,
                 aes(ymin = phi_lq, ymax = phi_uq, x = phi), lwd = 1) +
  # Add reference line with simulation value
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = expression("True" ~ phi), y = expression("Estimated"~ phi)) +
  expand_limits(x = c(0, 1), y = c(0, 1)) +
  theme_bw()

ggsave(phi_plot, 
       filename = "output/sim_studys1/phi_plot.png")



