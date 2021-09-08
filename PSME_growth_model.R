library(ForestPlot)
library(dplyr)
library(ggplot2)


#====================================
# Part 1. Selecting neighborhood size
#====================================

# Create vector of neighborhood radii
radii <- seq(2, 20, 2)

# Create data frame to store mse values
nb_rad_comp <- data.frame(
  radius = radii,
  mse = rep(NA, times = length(radii)),
  rsq = rep(NA, times = length(radii))
)

# Loop through radii, calculating mse for each
for(i in 1:length(radii)){
  
  # Extract neighborhood radius
  nb_rad <- radii[i]
  
  # Construct neighborhoods
  nbhds <- neighborhoods(mapping, radius = nb_rad)
  
  # Describe neighborhoods
  nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id", radius = nb_rad,
                                    densities = "angular")
  
  # Combine neighborhoods with their summaries
  nbhds <- nbhds %>%
    left_join(nbhd_summ, by = "tree_id")
  
  # Remove trees whose max neighborhood overlaps a plot boundary
  #nbhds <- nbhds %>%
  #  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad &
  #           y_coord >= nb_rad & y_coord <= 100 - nb_rad)
  nbhds <- nbhds %>%
    filter(x_coord >= max(radii) & x_coord <= 100 - max(radii) &
             y_coord >= max(radii) & y_coord <= 100 - max(radii))
  
  # Calculate annual growth for all trees
  growth <- growth_summary(tree)
  
  # Remove trees that were only measured once and/or had negative growth
  growth <- growth %>%
    filter(first_record != last_record &
             annual_growth >= 0)
  
  # Add growth data to neighborhoods
  nbhds <- nbhds %>%
    inner_join(growth %>% select(tree_id, size_corr_growth),
               by = "tree_id")
  
  # Add plot abiotic data
  nbhds <- nbhds %>%
    left_join(stand_abiotic, by = "stand_id")
  
  # Subset to focal species (PSME)
  nbhds <- nbhds %>%
    filter(species == "PSME")
  
  # Drop columns not needed for the model
  nbhds <- nbhds %>%
    select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))
  
  # Run model
  psme_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                           rare_comps = 100, density_suffix = "_angle_sum")
  
  # Store mean square error and R^2
  nb_rad_comp[i, "mse"] <- psme_mod$mod_coef$mse[1]
  nb_rad_comp[i, "rsq"] <- psme_mod$R_squared
}

# Plot relationship between neighborhood radius and model mse
ggplot(data = nb_rad_comp, aes(x = radius, y = mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  theme_classic()

#============================
# Part 2. Fitting final model
#============================

# Extract optimal neighborhood size
nb_rad <- nb_rad_comp %>%
  filter(mse == min(mse)) %>%
  pull(radius)

# Construct neighborhoods
nbhds <- neighborhoods(mapping, radius = nb_rad)

# Describe neighborhoods
nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id", radius = nb_rad,
                                  densities = "angular")

# Combine neighborhoods with their summaries
nbhds <- nbhds %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove trees whose max neighborhood overlaps a plot boundary
nbhds <- nbhds %>%
  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad &
           y_coord >= nb_rad & y_coord <= 100 - nb_rad)

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees that were only measured once and/or had negative growth
growth <- growth %>%
  filter(first_record != last_record &
           annual_growth >= 0)

# Add growth data to neighborhoods
nbhds <- nbhds %>%
  inner_join(growth %>% select(tree_id, size_corr_growth),
             by = "tree_id")

# Add plot abiotic data
nbhds <- nbhds %>%
  left_join(stand_abiotic, by = "stand_id")

# Subset to focal species (PSME)
nbhds <- nbhds %>%
  filter(species == "PSME")

# Drop columns not needed for the model
nbhds <- nbhds %>%
  select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))

# Run model
final_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                          rare_comps = 100, density_suffix = "_angle_sum")
