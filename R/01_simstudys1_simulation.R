#########################################################################
####                                                                 ####
####  Simulation study S1: a human movement-based spatial structure  ####
####                Script to simulate and explore data              ####
####                                                                 ####
#########################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load human movement-based connectivity coordinates ####
## NOTE: this code will only work if 00_human_movement_coords.R has been run
connect_coords <- fread("data/human_movement_coords.csv")


##### Step 1: Smooth surface simulation #####
### Scale coordinates to lie between [0, 1] 
connect_coords <- connect_coords %>% 
  mutate(shp,
         coord1_scale = (connect_coord1 - min(connect_coord1))/
           (max(connect_coord1) - min(connect_coord1))/20,
         coord2_scale = (connect_coord2 - min(connect_coord2))/
           (max(connect_coord2) - min(connect_coord2)))


## Use the smooth_create function to create smooth surface over coordinates
smooth_sim <- smooth_create(connect_coords$coord1_scale, 
                            connect_coords$coord2_scale)


# Centre estimates around 0 and scale
smooth_sim <- (smooth_sim - mean(smooth_sim)) * 20


## Plot smooth function on South Brazil map to check simulation
shp %>% 
  mutate(smooth = smooth_sim) %>% 
  ggplot( ) +
  geom_sf(aes(fill = smooth), lwd =  .05) +
  scale_fill_viridis_c(name = expression(u[i])) +
  theme_void()


#### Step 2: IID random effect simulation ####
## Simulate values from a normal distribution
n <- nrow(shp)

iid_re <- rnorm(n, sd = 1)


## Plot random variates on South Brazil map
shp %>% 
  mutate(iid = iid_re) %>% 
  ggplot() +
  geom_sf(aes(fill = iid), lwd =  .05) +
  scale_fill_viridis_c(name = "IID") + 
  theme_void()



#### Step 3: Combine random effects and simulate disease counts ####
## Create data frame to generate simulations with smooth and IID random effect 
df <- st_drop_geometry(shp) %>% 
  select(index, municip_code_ibge, lon, lat) %>% 
  mutate(connect_coord1 = connect_coords$connect_coord1,
         connect_coord2 = connect_coords$connect_coord2,
         smooth = smooth_sim,
         iid = iid_re,
         E = pop$pop_year/10^5)

## Create a list of phi values to use to combine random effects
# phi = 0 represents a completely random distribution (no spatial structure)
# phi = 1 represents a purely spatial structure
phi_list <- seq(0, 1, by = .1)


## Simulate counts from a poisson distribution using different contributions
# of each random element, determined by phi
df_sim_human <- lapply(phi_list, simulate_bym2)


#### Step 4: Explore simulations #### 
## Plot examples of simulated counts and smooth functions
# Phi = 0 (completely random structure)
phi0_Si_map <- df_sim_human[[1]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = log(lambda/E)), lwd = .05) +
  scale_fill_viridis_c(name = expression(S[i])) + 
  theme_void()

ggsave(phi0_Si_map, filename = "output/sim_studys1/Si_phi0_map.png")


y_phi0_map <- df_sim_human[[1]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = y), lwd = .05) +
  scale_fill_viridis_c(name = "Cases", trans = "log1p",
                       breaks = c(0, 10, 40)) + 
  expand_limits(fill = c(0, 45)) +
  theme_void()

ggsave(y_phi0_map, filename = "output/sim_studys1/y_phi0_map.png")


# Phi = 0.5 (equally structured and unstructured)
Si_phi05_map <- df_sim_human[[6]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = log(lambda/E)), lwd = .05) +
  scale_fill_viridis_c(name = expression(S[i])) + 
  theme_void()

ggsave(Si_phi05_map, filename = "output/sim_studys1/Si_phi05_map.png")


y_phi05_map <- df_sim_human[[6]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = y), lwd = .05) +
  scale_fill_viridis_c(name = "Cases", trans = "log1p",
                       breaks = c(0, 10, 40)) + 
  expand_limits(fill = c(0, 45)) +
  theme_void()

ggsave(y_phi05_map, filename = "output/sim_studys1/y_phi05_map.png")


# Phi = 1 (completely structured)
Si_phi1_map <- df_sim_human[[11]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = log(lambda/E)), lwd = .05) +
  scale_fill_viridis_c(name = expression(S[i])) + 
  theme_void()

ggsave(Si_phi1_map, filename = "output/sim_studys1/Si_phi1_map.png")


y_phi1_map <- df_sim_human[[11]] %>% 
  full_join(., shp,  by = c("index", "municip_code_ibge")) %>% 
  st_as_sf() %>% 
  ggplot( ) +
  geom_sf(aes(fill = y), lwd = .05) +
  scale_fill_viridis_c(name = "Cases", trans = "log1p",
                       breaks = c(0, 10, 40)) + 
  expand_limits(fill = c(0, 45)) +
  theme_void()

ggsave(y_phi1_map, filename = "output/sim_studys1/y_phi1_map.png")


#### Save data ####
write_rds(df_sim_human, "data/df_sim_human.rds")






