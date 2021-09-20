library(ForestPlot)
library(dplyr)
library(ggplot2)

#====================================
# Part 1. Selecting neighborhood size
#====================================

# Define focal species
focal_sps <- c("ABAM", "PSME", "TSHE")

# Create vector of neighborhood radii
radii <- seq(2, 20, 2)

# Create data frame to store mse values
nb_rad_comp <- data.frame(
  radius = radii,
  ABAM_mse = rep(NA, times = length(radii)),
  ABAM_rsq = rep(NA, times = length(radii)),
  PSME_mse = rep(NA, times = length(radii)),
  PSME_rsq = rep(NA, times = length(radii)),
  TSHE_mse = rep(NA, times = length(radii)),
  TSHE_rsq = rep(NA, times = length(radii))
)

# Loop through radii, calculating mse for each
for(i in 1:length(radii)){
  
  # Extract neighborhood radius
  nb_rad <- radii[i]
  
  # Construct neighborhoods
  nbhds <- neighborhoods(mapping, radius = nb_rad)
  
  # Describe neighborhoods
  nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                    radius = nb_rad, densities = "proportional")
  
  # Combine neighborhoods with their summaries
  nbhds <- nbhds %>%
    left_join(nbhd_summ, by = "tree_id")
  
  # Remove trees whose max neighborhood overlaps a plot boundary
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
  
  # Loop through focal species
  for(sps in focal_sps){
   
    # Subset nbhds to focal species
    one_sps <- nbhds %>%
      filter(species == sps)
    
    # Drop columns not needed for the model
    one_sps <- one_sps %>%
      select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))
    
    # Run model
    mod <- growth_model(one_sps, outcome_var = "size_corr_growth",
                        rare_comps = 100, density_suffix = "_density")
    
    # Store mean square error and R^2
    nb_rad_comp[i, paste0(sps, "_mse")] <- mod$mod_coef$mse[1]
    nb_rad_comp[i, paste0(sps, "_rsq")] <- mod$R_squared 
  }
}

# Plot relationship between neighborhood radius and model mse for each species
nb_comp_ABAM <- ggplot(data = nb_rad_comp, aes(x = radius, y = ABAM_mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  theme_classic()
nb_comp_PSME <- ggplot(data = nb_rad_comp, aes(x = radius, y = PSME_mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  theme_classic()
nb_comp_TSHE <- ggplot(data = nb_rad_comp, aes(x = radius, y = TSHE_mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  theme_classic()

# Chosen neighborhood sizes: ABAM = 12m, PSME = 10m, TSHE = 12m

#============================
# Part 2. Fitting final model
#============================

# Select focal species for fitting final model
sps <- "ABAM"

# Define optimal neighborhood size
if(sps == "PSME"){
  nb_rad <- 10
} else {
  nb_rad <- 12
}

# Construct neighborhoods
nbhds <- neighborhoods(mapping, radius = nb_rad)

# Describe neighborhoods
nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id", radius = nb_rad,
                                  densities = "proportional")

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

# Subset to focal species
nbhds <- nbhds %>%
  filter(species == sps)

# Drop columns not needed for the model
nbhds <- nbhds %>%
  select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))

# Run model
final_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                          iterations = 50, rare_comps = 100,
                          density_suffix = "_density")

# Format data for plotting
if(sps == "ABAM"){
  coef_summ <- final_mod$mod_coef %>%
    mutate(ABAM = 0.5 * (sps_compABAM + ABAM_density),
           ABLA = 0.5 * (sps_compABLA + ABLA_density),
           CANO = 0.5 * (sps_compCANO + CANO_density),
           PSME = 0.5 * (sps_compPSME + PSME_density),
           RARE = 0.5 * (sps_compRARE + RARE_density),
           TABR = 0.5 * (sps_compTABR + TABR_density),
           THPL = 0.5 * (sps_compTHPL + THPL_density),
           TSHE = 0.5 * (sps_compTSHE + TSHE_density),
           TSME = 0.5 * (sps_compTSME + TSME_density)) %>%
    select(ABAM, ABLA, CANO, PSME, RARE, TABR, THPL, TSHE, TSME)
  plot_data <- data.frame(
    competitor <- c("ABAM", "ABLA", "CANO", "PSME", "RARE", "TABR", "THPL",
                    "TSHE", "TSME"),
    means <- apply(coef_summ, MARGIN = 2, FUN = mean),
    sds <- apply(coef_summ, MARGIN = 2, FUN = sd)
  )
} else if(sps == "PSME"){
  coef_summ <- final_mod$mod_coef %>%
    mutate(ABAM = 0.5 * (sps_compABAM + ABAM_density),
           ABLA = 0.5 * (sps_compABLA + ABLA_density),
           PICO = 0.5 * (sps_compPICO + PICO_density),
           PIMO = 0.5 * (sps_compPIMO + PIMO_density),
           PSME = 0.5 * (sps_compPSME + PSME_density),
           RARE = 0.5 * (sps_compRARE + RARE_density),
           THPL = 0.5 * (sps_compTHPL + THPL_density),
           TSHE = 0.5 * (sps_compTSHE + TSHE_density)) %>%
    select(ABAM, ABLA, PICO, PIMO, PSME, RARE, THPL, TSHE)
  plot_data <- data.frame(
    competitor <- c("ABAM", "ABLA", "PICO", "PIMO", "PSME", "RARE", "THPL",
                    "TSHE"),
    means <- apply(coef_summ, MARGIN = 2, FUN = mean),
    sds <- apply(coef_summ, MARGIN = 2, FUN = sd)
  )
} else if(sps == "TSHE"){
  coef_summ <- final_mod$mod_coef %>%
    mutate(ABAM = 0.5 * (sps_compABAM + ABAM_density),
           ABLA = 0.5 * (sps_compABLA + ABLA_density),
           CANO = 0.5 * (sps_compCANO + CANO_density),
           PICO = 0.5 * (sps_compPICO + PICO_density),
           PIMO = 0.5 * (sps_compPIMO + PIMO_density),
           PSME = 0.5 * (sps_compPSME + PSME_density),
           RARE = 0.5 * (sps_compRARE + RARE_density),
           TABR = 0.5 * (sps_compTABR + TABR_density),
           THPL = 0.5 * (sps_compTHPL + THPL_density),
           TSHE = 0.5 * (sps_compTSHE + TSHE_density),
           TSME = 0.5 * (sps_compTSME + TSME_density)) %>%
    select(ABAM, ABLA, CANO, PICO, PIMO, PSME, RARE, TABR, THPL, TSHE, TSME)
  plot_data <- data.frame(
    competitor <- c("ABAM", "ABLA", "CANO", "PICO", "PIMO", "PSME", "RARE",
                    "TABR", "THPL", "TSHE", "TSME"),
    means <- apply(coef_summ, MARGIN = 2, FUN = mean),
    sds <- apply(coef_summ, MARGIN = 2, FUN = sd)
  )
}

# Create plot
ggplot(plot_data, aes(x = competitor, y = means)) +
  geom_point() +
  geom_errorbar(aes(ymin = means - sds, ymax = means + sds), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle(sps) +
  theme_classic()

#==================================
# Part 3. Creating multi-panel plot
#==================================
#library(gridExtra)
library(gtable)
library(grid)
g1 <- ggplotGrob(nb_comp_ABAM)
g2 <- ggplotGrob(nb_comp_PSME)
g3 <- ggplotGrob(nb_comp_TSHE)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

# Might want to look into ggarrange()