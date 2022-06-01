############################################################################
####                                                                    ####
####  Simulation study S2: a binomial distance-based spatial structure  ####
####                 Script to simulate and explore data                ####
####                                                                    ####
############################################################################

#### Load packages, data and functions ####
source("00_load_data_functions.R")


##### Step 1: Smooth surface simulation #####
### Scale coordinates to lie between [0, 1] 
shp <- shp %>% 
  mutate(lon_scale = (lon - min(lon))/(max(lon) - min(lon))/20,
         lat_scale = (lat - min(lat))/(max(lat) - min(lat)))


## Use the smooth_create function to create smooth surface over coordinates
smooth_sim <- smooth_create(shp$lon_scale, shp$lat_scale)


# Centre estimates around 0 and scale
smooth_sim <- (smooth_sim - mean(smooth_sim)) * 10


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


#### Step 3: Combine random effects and simulate binary outocmes ####
## Create data frame to generate simulations with smooth and IID random effect 
df <- st_drop_geometry(shp) %>% 
  select(index, municip_code_ibge, lon, lat) %>% 
  mutate(smooth = smooth_sim,
         iid = iid_re)

phi_list <- seq(0, 1, by = .1)

df_sim <- lapply(phi_list, simulate_binomial)


## Aggregate event outcomes to total number /20
df_binom <- lapply(df_sim, function(x) st_drop_geometry(x) %>% 
                     group_by(municip_code_ibge, municip_name, municip_code,
                              lon, lat, index, p_true, phi_true) %>% 
                     summarise(y = sum(y),
                               n_yrs = n()) %>% 
                     ungroup() )


#### Step 4: Explore simulations #### 
## Plot examples of simulated counts and smooth functions
# Phi = 0 (completely random structure)
y_map_phi0 <- ggplot(data = full_join(shp, df_binom[[1]], 
                                      by = "municip_code_ibge")) +
  geom_sf(aes(fill = y), lwd = .05) + 
  scale_fill_viridis_c(name = "y") +
  expand_limits(fill = c(0, 20))  +
  theme_void()

ggsave(y_map_phi0, 
       filename = "output/sim_studys2/y_map_phi0.png")


prob_map_phi0 <- ggplot(data = full_join(shp, df_binom[[1]], 
                                         by = "municip_code_ibge")) +
  geom_sf(aes(fill = p_true), lwd = .05) + 
  scale_fill_viridis_c(name = "Probability") +
  expand_limits(fill = c(0, 1)) +
  theme_void()

ggsave(prob_map_phi0, 
       filename = "output/sim_studys2/prob_map_phi0.png")


# Phi = 0.5 (equally structured and unstructured)
y_map_phi05 <- ggplot(data = full_join(shp, df_binom[[6]], 
                                       by = "municip_code_ibge")) +
  geom_sf(aes(fill = y), lwd = .05) + 
  scale_fill_viridis_c(name = "y") +
  expand_limits(fill = c(0, 20)) +
  theme_void()

ggsave(y_map_phi05, filename = "output/sim_studys2/y_map_phi05.png")


prob_map_phi05 <- ggplot(data = full_join(shp, df_binom[[6]], 
                                          by = "municip_code_ibge")) +
  geom_sf(aes(fill = p_true), lwd = .05) + 
  scale_fill_viridis_c(name = "Probability") +
  expand_limits(fill = c(0, 1)) +
  theme_void()

ggsave(prob_map_phi05, 
       filename = "output/sim_studys2/prob_map_phi05.png")


# Phi = 1 (completely structured)
y_map_phi1 <- ggplot(data = full_join(shp, df_binom[[11]], 
                                      by = "municip_code_ibge")) +
  geom_sf(aes(fill = y), lwd = .05) + 
  scale_fill_viridis_c(name = "y") +
  expand_limits(fill = c(0, 20)) +
  theme_void()

ggsave(y_map_phi1, filename = "output/sim_studys2/y_map_phi1.png")


prob_map_phi1 <- ggplot(data = full_join(shp, df_binom[[11]], 
                                         by = "municip_code_ibge")) +
  geom_sf(aes(fill = p_true), lwd = .05) + 
  scale_fill_viridis_c(name = "Probability") +
  expand_limits(fill = c(0, 1)) +
  theme_void()

ggsave(prob_map_phi1, 
       filename = "output/sim_studys2/prob_map_phi1.png")


#### Save data ####
write_rds(df_sim, file = "data/df_sim_binary_s2.rds")
write_rds(df_binom, file = "data/df_sim_binomial_2s.rds")


