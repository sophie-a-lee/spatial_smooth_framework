################################################################
####                                                        ####
####  Simulation study 2: two sources of spatial structure  ####
####      Script to create human movement coordinates       ####
####                                                        ####
################################################################


#### Load packages, data and functions ####
source("00_load_data_functions.R")


#### Step 1: Extract connectivity matrix from REGIC shapefile ####
## Load data with links between cities in South Brazil
regic_link <- read_rds("data/regic_network_south.rds")


# Convert data into an sfnetwork object
regic_net <- as_sfnetwork(regic_link, directed = F)


## Plot map with connections
regic_south <- ggplot() +
  geom_sf(data = shp, fill = NA, lwd = .25) +
  geom_sf(data = st_as_sf(regic_link, "edges"),
          lwd = .1)+
  theme_void()

ggsave(regic_south, 
       filename = "output/sim_study2/regic_map.png")


# Extract nodes of the network (cities from links)
regic_nodes <- st_as_sf(regic_net, "nodes")


# Extract municipality codes for the origin nodes 
# (returns first if more than one, returns NA if number not in the origin column)
start <- lapply(
  X = st_equals(regic_nodes, st_startpoint(regic_link)),
  FUN = function(x) regic_link$origin_code[x[1]]
)


# Extract municipality codes for the destination nodes 
# (returns first if more than one, returns NA if number not in the origin column)
end <- lapply(
  X = st_equals(regic_nodes, st_endpoint(regic_link)),
  FUN = function(i) regic_link$dest_code[i[1]]
)


# Combine municipality codes from origin and destination nodes
idxs <- mapply(dplyr::coalesce, start, end)

# Add names to the network
regic_net <- regic_net %>%
  mutate(name = idxs)


# Create sparse binary matrix with connections
regic_adj <- igraph::as_adjacency_matrix(regic_net, names = TRUE)

# Order matrix to preserve order when combining with data
regic_adj <- regic_adj[order(rownames(regic_adj)), 
                       order(colnames(regic_adj))]


# Convert REGIC matrix to a nb list compatible with models
nb_net <- mat2listw(regic_adj, style = "B")$neighbours


## Create a binary neighbourhood matrix (1 = connected, 0  =  not)
nb_regic <- nb2mat(nb_net, style = "B")
colnames(nb_regic) <- rownames(nb_regic)

# Order according to municipality codes to fit with data
nb_regic <- nb_regic[order(rownames(nb_regic)), 
                     order(colnames(nb_regic))]


### Step 2: Extract continuous measure of connectivity based on REGIC, distance and population ####
## Create df with connectivity status between each pair of cities. Include
# populations at each origin  source city, and the distance between them 
df_nb <- as.data.frame(nb_regic)
df_nb$source_code <-  rownames(df_nb)


# Convert connectivity matrix into long-format (row per pair of cities)
df_connect <- reshape::melt(df_nb) %>%
  mutate(source_code = as.numeric(source_code),
         dest_code = as.numeric(as.character(variable)),
         connected = value) %>%
  filter(source_code != dest_code) %>%
  dplyr::select(source_code, dest_code, connected)


## Calculate distance between cities
# Set up blank matrix, add source and dest codes as row/col names
dist_mat <- matrix(nrow = nrow(pop), ncol = nrow(pop))
colnames(dist_mat) <- pop$municip_code_ibge
rownames(dist_mat) <- pop$municip_code_ibge


# Add distance between cities into matrix
# WARNING: this command takes a long time to run
for(i in rownames(dist_mat)) {
  for(j in colnames(dist_mat)){
    dist_mat[which(rownames(dist_mat) == i), 
             which(colnames(dist_mat) == j)] <- 
      st_distance(pop[which(rownames(dist_mat) == i),], 
                  pop[which(colnames(dist_mat) == j),])
  }
}


# Convert distance matrix into long df
df_dist <- as.data.frame(dist_mat)
df_dist$source_code <-  rownames(df_dist)

df_dist <- reshape::melt(df_dist) %>%
  mutate(source_code = as.numeric(source_code),
         dest_code = as.numeric(as.character(variable)),
         distance = value) %>%
  filter(source_code != dest_code) %>%
  dplyr::select(source_code, dest_code, distance)


## Combine population, distance + REGIC connectivity data
df_gam <- full_join(df_connect, df_dist, by = c("source_code", "dest_code")) %>%
  full_join(., st_drop_geometry(pop), 
            by = c("source_code" = "municip_code_ibge")) %>%
  rename(pop_source = pop_year) %>%
  full_join(., pop, by = c("dest_code" = "municip_code_ibge")) %>%
  rename(pop_dest = pop_year) %>% 
  dplyr::select(-geometry)


## Fit GAM to connectivity indicator with smooth fn of pop + dist 
# Use a tensor product as variables are on different scales
connect_model <- bam(connected ~ te(pop_source, pop_dest, distance, bs = "tp"),
                     data = df_gam, family = binomial, method = "REML")


## Extract 'probability of connectivity' from GAM output 
df_pred <- pred_binomial_prob(connect_model)


## Convert prediction df into matrix
pred_mat <- df_pred %>%
  dplyr::select(-connected) %>%
  # Convert to wide format (row per destination city, column per source city)
  pivot_wider(names_from = dest_code, values_from = prob_connect, 
              values_fill = 0) %>%
  column_to_rownames(var="source_code")


## Apply multidimensional scaling to convert probability matrix into 2-d coordinates
msd_scale <- cmdscale(pred_mat, eig = T, k = 2)
# Plot to view coordinate system
plot(msd_scale$points[,1], msd_scale$points[,2])

# Convert into a df with 'connectivity coordinates' per city
connect_coords <-  as.data.frame(msd_scale$points)
names(connect_coords) <- c("connect_coord1", "connect_coord2")
connect_coords$municip_code_ibge <- as.numeric(rownames(connect_coords))


#### Save connectivity coordinates for model simulations ####
fwrite(connect_coords, file = "data/human_movement_coords.csv")





