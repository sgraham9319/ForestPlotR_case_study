library(ForestPlotR)
library(dplyr)
library(ggplot2)
library(ggpubr)

#====================================
# Part 1. Selecting neighborhood size
#====================================

# Fit models with different neighborhood sizes
nbhd_size_comp <- select_nbhd_size(radii = seq(2, 20, 2), map_data = mapping,
                                   model_process = "growth", tree_data = tree,
                                   abiotic_data = stand_abiotic,
                                   focal_sps = c("ABAM", "PSME", "TSHE"),
                                   dens_type = "proportional",
                                   max_x = 100, max_y = 100)

# View results
nbhd_size_comp$ABAM_plot
nbhd_size_comp$PSME_plot
nbhd_size_comp$TSHE_plot

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

# Remove evenness variable because it will be NA when species richness = 1
nbhd_summ <- nbhd_summ %>%
  select(-evenness)

# Combine neighborhoods with their summaries
nbhds <- nbhds %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove trees whose max neighborhood overlaps a plot boundary
nbhds <- nbhds %>%
  filter(x_coord >= 12 & x_coord <= 100 - 12 &
           y_coord >= 12 & y_coord <= 100 - 12)

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees that were only measured once and/or had negative growth
growth <- growth %>%
  filter(first_record != last_record & annual_growth >= 0)

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
ABAM_mod <- growth_mort_model(nbhds, outcome_var = "size_corr_growth",
                              rare_comps = 100, density_suffix = "_density")

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
  competitor = names(coef_summ),
  growth_response = t(coef_summ)[,1]
)

# Identify competitor species with coefficient of zero
pt_cols <- if_else(plot_data_ABAM$growth_response == 0, "red", "black")

# Create plot
sps_int_ABAM <- ggplot(plot_data_ABAM, aes(x = competitor,
                                           y = growth_response)) +
  geom_point(size = 3, color = pt_cols) +
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

# Remove evenness variable because it will be NA when species richness = 1
nbhd_summ <- nbhd_summ %>%
  select(-evenness)

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
PSME_mod <- growth_mort_model(nbhds, outcome_var = "size_corr_growth",
                              rare_comps = 100, density_suffix = "_density")

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
  competitor = names(coef_summ),
  growth_response = t(coef_summ)[,1]
)

# Identify competitor species with coefficient of zero
pt_cols <- if_else(plot_data_PSME$growth_response == 0, "red", "black")

# Create plot
sps_int_PSME <- ggplot(plot_data_PSME, aes(x = competitor,
                                           y = growth_response)) +
  geom_point(size = 3, color = pt_cols) +
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

# Remove evenness variable because it will be NA when species richness = 1
nbhd_summ <- nbhd_summ %>%
  select(-evenness)

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
TSHE_mod <- growth_mort_model(nbhds, outcome_var = "size_corr_growth",
                              rare_comps = 100, density_suffix = "_density")

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
  competitor = names(coef_summ),
  growth_response = t(coef_summ)[,1]
)

# Identify competitor species with coefficient of zero
pt_cols <- if_else(plot_data_TSHE$growth_response == 0, "red", "black")

# Create plot
sps_int_TSHE <- ggplot(plot_data_TSHE, aes(x = competitor,
                                           y = growth_response)) +
  geom_point(size = 3, color = pt_cols) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Neighbor species") +
  ylab("Growth response") +
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -0.0027, ymax = 0.01,
           alpha = 0.3) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.0025, 0.0105)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#==================================
# Part 3. Creating multi-panel plot
#==================================

# Clean the neighborhood size comparison figures
nb_comp_ABAM <- nbhd_size_comp$ABAM_plot +
  geom_vline(xintercept = 12, linetype='dashed') +
  ggtitle("Abies amabilis") +
  theme(plot.title = element_text(face = "italic"))
nb_comp_PSME <- nbhd_size_comp$PSME_plot +
  geom_vline(xintercept = 10, linetype='dashed') +
  ggtitle("Pseudotsuga menziesii") +
  theme(plot.title = element_text(face = "italic"))
nb_comp_TSHE <- nbhd_size_comp$TSHE_plot +
  geom_vline(xintercept = 12, linetype='dashed') +
  ggtitle("Tsuga heterophylla") +
  theme(plot.title = element_text(face = "italic"))

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
