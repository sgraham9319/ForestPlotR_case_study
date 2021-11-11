library(ForestPlot)
library(dplyr)
library(ggplot2)
library(ggpubr)

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

# Loop through radii, calculating mse for each (this may take 2 minutes to run)
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
  ggtitle("Abies amabilis") +
  theme_classic() +
  theme(plot.title = element_text(face = "italic"))
nb_comp_PSME <- ggplot(data = nb_rad_comp, aes(x = radius, y = PSME_mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  ggtitle("Pseudotsuga menziesii") +
  theme_classic() +
  theme(plot.title = element_text(face = "italic"))
nb_comp_TSHE <- ggplot(data = nb_rad_comp, aes(x = radius, y = TSHE_mse)) +
  geom_line(col = "green") +
  labs(x = "Neighborhood radius (m)", y = "Mean square error") +
  ggtitle("Tsuga heterophylla") +
  theme_classic() +
  theme(plot.title = element_text(face = "italic"))

# Chosen neighborhood sizes: ABAM = 12m, PSME = 10m, TSHE = 12m

#============================
# Part 2. Fitting final model
#============================

#----------------
# 2.1. ABAM model
#----------------

# Construct neighborhoods
nbhds <- neighborhoods(mapping, radius = 12)

# Describe neighborhoods
nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                  radius = 12, densities = "proportional")

# Combine neighborhoods with their summaries
nbhds <- nbhds %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove trees whose max neighborhood overlaps a plot boundary
nbhds <- nbhds %>%
  filter(x_coord >= 12 & x_coord <= 100 - 12 &
           y_coord >= 12 & y_coord <= 100 - 12)

# Add growth data to neighborhoods
nbhds <- nbhds %>%
  inner_join(growth %>% select(tree_id, size_corr_growth), by = "tree_id")

# Add plot abiotic data
nbhds <- nbhds %>%
  left_join(stand_abiotic, by = "stand_id")

# Subset to focal species and drop unneeded columns
nbhds <- nbhds %>%
  filter(species == "ABAM") %>%
  select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))

# Fit model
ABAM_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                         iterations = 10, rare_comps = 100,
                         density_suffix = "_density")

# Format data for plotting
coef_summ <- ABAM_mod$mod_coef %>%
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
plot_data_ABAM <- data.frame(
  competitor = c("ABAM", "ABLA", "CANO", "PSME", "RARE", "TABR", "THPL",
                  "TSHE", "TSME"),
  means = apply(coef_summ, MARGIN = 2, FUN = mean),
  sds = apply(coef_summ, MARGIN = 2, FUN = sd)
)

# Create plot
sps_int_ABAM <- ggplot(plot_data_ABAM, aes(x = competitor, y = means)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = means - sds, ymax = means + sds), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Neighbor species") +
  ylab("Growth response") +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -0.002, ymax = 0.004,
           alpha = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#----------------
# 2.2. PSME model
#----------------

# Construct neighborhoods
nbhds <- neighborhoods(mapping, radius = 10)

# Describe neighborhoods
nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                  radius = 10, densities = "proportional")

# Combine neighborhoods with their summaries
nbhds <- nbhds %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove trees whose max neighborhood overlaps a plot boundary
nbhds <- nbhds %>%
  filter(x_coord >= 10 & x_coord <= 100 - 10 &
           y_coord >= 10 & y_coord <= 100 - 10)

# Add growth data to neighborhoods
nbhds <- nbhds %>%
  inner_join(growth %>% select(tree_id, size_corr_growth), by = "tree_id")

# Add plot abiotic data
nbhds <- nbhds %>%
  left_join(stand_abiotic, by = "stand_id")

# Subset to focal species and drop unneeded columns
nbhds <- nbhds %>%
  filter(species == "PSME") %>%
  select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))

# Fit model
PSME_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                         iterations = 10, rare_comps = 100,
                         density_suffix = "_density")

# Format data for plotting
coef_summ <- PSME_mod$mod_coef %>%
  mutate(ABAM = 0.5 * (sps_compABAM + ABAM_density),
         ABLA = 0.5 * (sps_compABLA + ABLA_density),
         PICO = 0.5 * (sps_compPICO + PICO_density),
         PIMO = 0.5 * (sps_compPIMO + PIMO_density),
         PSME = 0.5 * (sps_compPSME + PSME_density),
         RARE = 0.5 * (sps_compRARE + RARE_density),
         THPL = 0.5 * (sps_compTHPL + THPL_density),
         TSHE = 0.5 * (sps_compTSHE + TSHE_density)) %>%
  select(ABAM, ABLA, PICO, PIMO, PSME, RARE, THPL, TSHE)
plot_data_PSME <- data.frame(
  competitor = c("ABAM", "ABLA", "PICO", "PIMO", "PSME", "RARE", "THPL",
                  "TSHE"),
  means = apply(coef_summ, MARGIN = 2, FUN = mean),
  sds = apply(coef_summ, MARGIN = 2, FUN = sd)
)

# Create plot
sps_int_PSME <- ggplot(plot_data_PSME, aes(x = competitor, y = means)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = means - sds, ymax = means + sds), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Neighbor species") +
  ylab("Growth response") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -0.0015, ymax = 0.0045,
           alpha = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#----------------
# 2.3. TSHE model
#----------------

# Construct neighborhoods
nbhds <- neighborhoods(mapping, radius = 12)

# Describe neighborhoods
nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                  radius = 12, densities = "proportional")

# Combine neighborhoods with their summaries
nbhds <- nbhds %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove trees whose max neighborhood overlaps a plot boundary
nbhds <- nbhds %>%
  filter(x_coord >= 12 & x_coord <= 100 - 12 &
           y_coord >= 12 & y_coord <= 100 - 12)

# Add growth data to neighborhoods
nbhds <- nbhds %>%
  inner_join(growth %>% select(tree_id, size_corr_growth), by = "tree_id")

# Add plot abiotic data
nbhds <- nbhds %>%
  left_join(stand_abiotic, by = "stand_id")

# Subset to focal species and drop unneeded columns
nbhds <- nbhds %>%
  filter(species == "TSHE") %>%
  select(-c(stand_id, species, dbh, abh, x_coord, y_coord, id_comp, abh_comp))

# Fit model
TSHE_mod <- growth_model(nbhds, outcome_var = "size_corr_growth",
                         iterations = 10, rare_comps = 100,
                         density_suffix = "_density")

# Format data for plotting
coef_summ <- TSHE_mod$mod_coef %>%
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
plot_data_TSHE <- data.frame(
  competitor = c("ABAM", "ABLA", "CANO", "PICO", "PIMO", "PSME", "RARE",
                  "TABR", "THPL", "TSHE", "TSME"),
  means = apply(coef_summ, MARGIN = 2, FUN = mean),
  sds = apply(coef_summ, MARGIN = 2, FUN = sd)
)

# Create plot
sps_int_TSHE <- ggplot(plot_data_TSHE, aes(x = competitor, y = means)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = means - sds, ymax = means + sds), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Neighbor species") +
  ylab("Growth response") +
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -0.0027, ymax = 0.01,
           alpha = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#==================================
# Part 3. Creating multi-panel plot
#==================================

# Initiate plot saving
jpeg(filename = "Figure2.jpg",
     type = "cairo",
     units = "px",
     width = 6000,
     height = 4000,
     pointsize = 12,
     res = 500)

# Create plot
ggarrange(nb_comp_ABAM, nb_comp_PSME, nb_comp_TSHE,
          sps_int_ABAM, sps_int_PSME, sps_int_TSHE,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)

# Save plot
dev.off()
