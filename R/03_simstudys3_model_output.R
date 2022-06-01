############################################################################
####                                                                    ####
####  Simulation study S3: a binomial model with 2 spatial structures   ####
####                Script to explore model output                      ####
####                                                                    ####
############################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load simulated data and models ####
## NOTE: this code will only work if 01_simstudys3_simulation.R and
## 02_simstudys3_run_models.R has been run
df_sim <- read_rds("data/df_sim_binary_s3.rds")
df_binom <- read_rds("data/df_sim_binomial_s3.rds")

model_list <- read_rds("output/sim_studys3/model_list.rds")


#### Extract predicted intercept from models ####
## Spatial smooth model
b_est <- lapply(model_list, b_extract_3phi)

b_table <- reduce(b_est, rbind) %>% 
  # Add true phi value
  mutate(phi_human = seq(0, .8, by = .1))


#### Plot intercept estimates for models ####
b_line <- ggplot(data = b_table) +
  # Plot estimates + 95% CI
  geom_point(aes(x = phi_human, y = b_est)) +
  geom_linerange(aes(ymin = b_lq, ymax = b_uq, x = phi_human), lwd = 1) +
  # Add reference line with true value (b = 0)
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(name = expression(phi[human]), 
                     breaks = seq(0, .8, by = .1)) +
  scale_y_continuous(name = "Intercept estimate") +
  theme_bw()


# Save intercept plot
ggsave(b_line, filename = "output/sim_studys3/intercept_line.png")


#### Extract estimates for phi/mixing parameters ####
model_vars <- lapply(model_list, extract_vars_3re)


## Combine the estimated phi from each model and calculate mean and 95% CI
phi_est <- reduce(model_vars, rbind) %>% 
  #  Add true phi value
  mutate(phi_human = rep(seq(0, .8, by = .1), 
                         each = nrow(model_vars[[1]])),
         phi_iid = .2, 
         phi_dist = 1 - (phi_human + phi_iid)) %>% 
  group_by(phi_iid, phi_dist) %>% 
  summarise(phi_dist_est = mean(propn_var_dist),
            phi_dist_lq = quantile(propn_var_dist, .025),
            phi_dist_uq = quantile(propn_var_dist, .975),
            phi_human_est = mean(propn_var_human),
            phi_human_lq = quantile(propn_var_human, .025),
            phi_human_uq = quantile(propn_var_human, .975),
            phi_iid_est = mean(propn_var_iid),
            phi_iid_lq = quantile(propn_var_iid, .025),
            phi_iid_uq = quantile(propn_var_iid, .975)) %>% 
  mutate(phi_human = 1 - (phi_dist + phi_iid)) %>% 
  ungroup()


#### Plot mean & 95% CI phi estimates ####
## Distance-based random effect
phi_dist_line <- ggplot(data = phi_est) +
  geom_point(aes(x = phi_dist, y = phi_dist_est)) +
  geom_linerange(aes(ymin = phi_dist_lq, ymax = phi_dist_uq, 
                     x = phi_dist), lwd = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-.05, .85), y = c(-.01, 1), expand = F) +
  labs(x = expression("True" ~ phi["distance"]),
       y = expression("Estimated" ~ phi["distance"])) +
  scale_x_continuous(breaks = seq(0, .8, by = .1)) +
  theme_bw()


ggsave(phi_dist_line, 
       filename = "output/sim_studys3/phi_line_dist.png")


## Human movement-based random effect
phi_human_line <- ggplot(data = phi_est) +
  geom_point(aes(x = phi_human, y = phi_human_est)) +
  geom_linerange(aes(ymin = phi_human_lq, ymax = phi_human_uq, 
                     x = phi_human), lwd = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-.05, .85), y = c(-.01, 1), expand = F) +
  labs(x = expression("True" ~ phi["human"]),
       y = expression("Estimated" ~ phi["human"])) +
  scale_x_continuous(breaks = seq(0, .8, by = .1)) +
  theme_bw()


ggsave(phi_human_line, 
       filename = "output/sim_studys3/phi_line_human.png")


## Unstructured random effects
phi_iid_line <- ggplot(data = phi_est) +
  geom_point(aes(x = phi_human, y = phi_iid_est)) +
  geom_linerange(aes(ymin = phi_iid_lq, ymax = phi_iid_uq, 
                     x = phi_human), lwd = 1) +
  geom_hline(yintercept = .2, linetype = "dashed") +
  coord_cartesian(xlim = c(-.05, .85), y = c(-.01, 1), expand = F) +
  labs(x = expression("True" ~ phi["human"]),
       y = expression("Estimated" ~ phi["iid"])) +
  scale_x_continuous(breaks = seq(0, .8, by = .1)) +
  theme_bw()


ggsave(phi_iid_line, 
       filename = "output/sim_studys3/phi_line_iid.png")



