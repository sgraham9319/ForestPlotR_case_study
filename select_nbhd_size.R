test <- select_nbhd_size(radii = seq(3, 15, 3), map_data = mapping,
                         growth_data = tree, abiotic_data = stand_abiotic,
                         focal_sps = c("ABAM", "PSME", "TSHE"),
                         dens_type = "proportional", max_x = 100, max_y = 100)
test$ABAM_plot
nb_comp_ABAM <- test$ABAM_plot +
  ggtitle("Abies amabilis") +
  theme(plot.title = element_text(face = "italic"))


select_nbhd_size <- function(radii, map_data, growth_data,
                             abiotic_data = NULL, focal_sps, dens_type,
                             max_x, max_y, rare_sps = 100){
  
  # Create data frame to store results
  nb_rad_comp <- as.data.frame(matrix(NA, nrow = length(radii),
                                      ncol = length(focal_sps) + 1))
  names(nb_rad_comp) <- c("radius", paste0(focal_sps, "_mse"))
  nb_rad_comp$radius <- radii
  
  # Loop through the neighborhood radii
  for(i in 1:length(radii)){
    
    # Select neighborhood radius
    nb_rad <- radii[i]
    
    # Construct neighborhoods
    nbhds <- neighborhoods(map_data, radius = nb_rad)
    
    # Describe neighborhoods
    nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                      radius = nb_rad, densities = dens_type)
    
    # Combine neighborhoods with their summaries
    nbhds <- nbhds %>%
      left_join(nbhd_summ, by = "tree_id")
    
    # Remove trees whose max neighborhood overlaps a plot boundary
    nbhds <- nbhds %>%
      filter(x_coord >= max(radii) & x_coord <= max_x - max(radii) &
               y_coord >= max(radii) & y_coord <= max_y - max(radii))
    
    # Calculate annual growth for all trees
    growth <- growth_summary(growth_data)
    
    # Remove trees that were only measured once and/or had negative growth
    growth <- growth %>%
      filter(first_record != last_record & annual_growth >= 0)
    
    # Add growth data to neighborhoods
    nbhds <- nbhds %>%
      inner_join(growth %>% select(tree_id, size_corr_growth),
                 by = "tree_id")
    
    # Add plot abiotic data, if provided
    if(!is.null(abiotic_data)){
      nbhds <- nbhds %>%
        left_join(abiotic_data, by = "stand_id")
    }
    
    # Loop through focal species
    for(sps in focal_sps){
      
      # Subset nbhds to focal species
      one_sps <- nbhds %>%
        filter(species == sps)
      
      # Drop columns not needed for the model
      one_sps <- one_sps %>%
        select(-c(stand_id, species, dbh, abh, x_coord, y_coord,
                  id_comp, abh_comp))
      
      # Run model
      if(dens_type == "angular"){
        mod <- growth_model(one_sps, outcome_var = "size_corr_growth",
                            rare_comps = rare_sps,
                            density_suffix = "_angle_sum")
      } else{
        mod <- growth_model(one_sps, outcome_var = "size_corr_growth",
                            rare_comps = rare_sps,
                            density_suffix = "_density")
      }
      
      # Store mean square error
      nb_rad_comp[i, paste0(sps, "_mse")] <- mod$mod_coef$mse[1] 
    }
  }
  
  # Create list to store overall output
  output_list <- as.list(rep(NA, times = length(focal_sps) + 1))
  names(output_list) <- c("mse_vals", paste0(focal_sps, "_plot"))
  
  # Make table of mse values the first element of output list
  output_list$mse_vals <- nb_rad_comp
  
  # Make a graph for each species, storing in output list
  for(a in 1:length(focal_sps)){
    
    # the local() function is used so that ggplot code is forced to evaluate
    # locally (i.e. within the loop)
    output_list[[a + 1]] <- local({
      local_a <- a
      sps_plot <- ggplot(data = nb_rad_comp,
                         aes(x = radius,
                             y = get(paste0(focal_sps[local_a], "_mse")))) +
        geom_line(col = "green") +
        labs(x = "Neighborhood radius (m)", y = "Mean square error") +
        theme_classic()
    })
  }
  
  # Return the output list
  return(output_list)
}
