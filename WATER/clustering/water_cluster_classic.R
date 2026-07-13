water_cluster_classic <- function(
    water,
    k_max = NULL,
    method = "centroid",
    r_threshold = 0.8
) {
  
  method_specified <- !missing(method)
  
  # Load required WATER functions
  source("clustering/water_cluster_initial.R")
  source("clustering/water_cut.R")   
  source("clustering/water_label.R")
  
  # ----------------------------
  # Step 1: Silhouette + heatmaps
  # ----------------------------
  message("  [2.1] Initial clustering (silhouette + heatmaps)")
  water <- water_cluster_initial(water, k_max = k_max)
  
  # ----------------------------
  # Step 2: Prompt user for k_final
  # ----------------------------
  message("  [2.2] Select number of clusters (k)")
  k_recommended <- water$cluster_classic$k_recommended
  k_default <- max(k_recommended)
  
  message("\nRecommended k values: ", paste(k_recommended, collapse = ", "))
  message("Press Enter to use the default (k = ", k_default, ")")
  message("Or type a k value and press Enter: ")
  
  input <- readline()
  
  if (nchar(trimws(input)) == 0) {
    k_final <- k_default
    message("Using default k = ", k_final)
  } else {
    k_final <- suppressWarnings(as.integer(trimws(input)))
    if (is.na(k_final)) {
      warning("Invalid input. Defaulting to k = ", k_default)
      k_final <- k_default
    }
  }
  
  # ----------------------------
  # Step 3: Cut clusters at chosen k
  # ----------------------------
  message("  [2.3] Cutting clusters")
  water <- water_cut(water, k_final = k_final)
  
  # ----------------------------
  # Step 4: Prompt user for labeling method
  # ----------------------------
  message("  [2.4] Label trajectories")
  if (!method_specified) {
    
    message("How would you like to collapse clusters into trajectories?")
    message("  1. Centroid clustering (default, r_threshold = ", r_threshold, ")")
    message("  2. Manual label_map")
    message("  3. None - keep all ", k_final, " clusters as separate trajectories (not recommended for large k)")
    message("Press Enter for default (centroid), or type 1, 2, or 3: ")
    
    method_input <- readline()
    
    if (nchar(trimws(method_input)) == 0 || trimws(method_input) == "1") {
      method <- "centroid"
    } else if (trimws(method_input) == "2") {
      method <- "manual"
    } else if (trimws(method_input) == "3") {
      message("Warning: this will keep all ", k_final,
              " clusters as individual trajectories.")
      message("Are you sure? Type 'yes' to confirm, or press Enter to use centroid clustering: ")
      confirm <- trimws(readline())
      if (tolower(confirm) == "yes") {
        method <- "none"
      } else {
        message("Defaulting to centroid clustering.")
        method <- "centroid"
      }
    } else {
      message("Invalid input. Defaulting to centroid clustering.")
      method <- "centroid"
    }
  }
  
  # ----------------------------
  # Step 5: Label trajectories
  # ----------------------------
  water <- water_label(
    water,
    method      = method,
    r_threshold = r_threshold
  )

  message("Clustering complete.")
  
  return(water)
}