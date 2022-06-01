############################################################################
####                                                                    ####
####  Simulation study S2: a binomial model with 2 spatial structures   ####
####                 Script to simulate and explore data                ####
####                                                                    ####
############################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Load human movement-based connectivity coordinates ####
## NOTE: this code will only work if 00_human_movement_coords.R has been run
connect_coords <- fread("data/human_movement_coords.csv")


#### Step 1: distance-based smooth surface simulation ####
### Scale coordinates to lie between [0, 1] 
shp <- shp %>% 
  mutate(lon_scale = (lon - min(lon))/(max(lon) - min(lon))/20,
         lat_scale = (lat - min(lat))/(max(lat) - min(lat)))


## Use the smooth_create function to create smooth surface over coordinates
smooth_sim_dist <- smooth_create(shp$lon_scale, shp$lat_scale)


# Centre estimates around 0 and scale
smooth_sim_dist <- (smooth_sim_dist - mean(smooth_sim_dist)) * 10


## Plot smooth function on South Brazil map to check simulation
shp %>% 
  mutate(smooth = smooth_sim_dist) %>% 
  ggplot( ) +
  geom_sf(aes(fill = smooth), lwd =  .05) +
  scale_fill_viridis_c(name = expression(u[i])) +
  theme_void()


#### Step 2: human movement-based smooth surface simulation ####
## Scale the coordinates to use function
connect_coords <- connect_coords %>% 
  mutate(shp,
         coord1_scale = (connect_coord1 - min(connect_coord1))/
           (max(connect_coord1) - min(connect_coord1))/20,
         coord2_scale = (connect_coord2 - min(connect_coord2))/
           (max(connect_coord2) - min(connect_coord2)))


## Use the smooth_create function to create smooth surface over coordinates
# (see ?te for more details)
smooth_sim_human <- smooth_create(connect_coords$coord1_scale, 
                                  connect_coords$coord2_scale)

# Centre estimates around 0 and scale 
smooth_sim_human <- (smooth_sim_human - mean(smooth_sim_human)) * 20


## Plot smooth function on South Brazil map to check simulation
shp %>% 
  mutate(smooth = smooth_sim_human) %>% 
  ggplot( ) +
  geom_sf(aes(fill = smooth), lwd =  .05) +
  scale_fill_viridis_c(name = "Smooth")


#### Step 3: IID random effect simulation ####
## Simulate values from a normal distribution
n <- nrow(shp)
iid_re <- rnorm(n, sd = 1)


## Plot random variates on South Brazil map
shp %>% 
  mutate(iid = iid_re) %>% 
  ggplot() +
  geom_sf(aes(fill = iid), lwd =  .05) +
  scale_fill_viridis_c(name = "IID")


#### Step 4: Simulate data using different phi estimates ####
## Combine simulated random terms into one dataset
df <- st_drop_geometry(shp) %>% 
  mutate(connect_coord1 = connect_coords$connect_coord1,
         connect_coord2 = connect_coords$connect_coord2,
         smooth_dist = smooth_sim_dist,
         smooth_human = smooth_sim_human,
         iid = iid_re)


## List phi values for simulations, fix phi_iid = .1
phi_iid <- .2

phi_human <- seq(0, .8, by = .1)

df_phi <- list(rep(data.table(phi_dist = NA,
                              phi_human = NA,
                              phi_iid = NA), 
                   length(phi_human)))


# Create list of phis to apply
for(i in 1:length(phi_human)) {
  df_phi[[i]] <- data.table(phi_dist = 1 - (phi_human[i] + phi_iid),
                            phi_human = phi_human[i],
                            phi_iid = phi_iid)
}



df_sim <- lapply(df_phi, simulate_binomial_3phi)


## Aggregate event outcomes to total number /20
df_binom <- lapply(df_sim, function(x) x %>% 
                     group_by(municip_code_ibge, municip_name, municip_code,
                              lon, lat, connect_coord1, connect_coord2,
                              index, p_true, phi_dist,
                              phi_human, phi_iid) %>% 
                     summarise(y = sum(y),
                               n_yrs = n()) %>% 
                     ungroup() )


#### Step 5: Explore simulations #### 
# phi_dist = .4, phi_human = .4,  phi_iid = .2
y_map <- ggplot(data = full_join(shp, df_binom[[5]], 
                                 by = "municip_code_ibge")) +
  geom_sf(aes(fill = y), lwd = .05) + 
  scale_fill_viridis_c(name = "y") +
  expand_limits(fill = c(0, 20)) +
  theme_void()

ggsave(y_map, filename = "output/sim_studys3/y_map.png")


prob_map <- ggplot(data = full_join(shp, df_binom[[5]], 
                                    by = "municip_code_ibge")) +
  geom_sf(aes(fill = p_true), lwd = .05) + 
  scale_fill_viridis_c(name = "Probability") +
  expand_limits(fill = c(0, 1)) +
  theme_void()

ggsave(prob_map, 
       filename = "output/sim_studys3/prob_map.png")


#### Save data ####
write_rds(df_sim, file = "data/df_sim_binary_s3.rds")
write_rds(df_binom, file = "data/df_sim_binomial_s3.rds")


