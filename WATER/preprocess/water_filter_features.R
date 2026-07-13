#' Filter and prepare features for WATER analysis
#'
#' Averages replicates, optionally filters features, and bins features
#' by differential abundance using log2FC values.
#' 
water_filter_features <- function(
    water,
    feature_list = NULL,
    log2fc_threshold = 1,
    average_replicates = TRUE
) {
  
  # ----------------------------
  # Use replicate-level filtered data
  # ----------------------------
  if (is.null(water$data_filt)) {
    stop("Run water_filter_expression() first")
  }
  
  data <- water$data_filt   
  n_start <- nrow(data)
  
  # ----------------------------
  # 1. Optional feature filtering
  # ----------------------------
  if (!is.null(feature_list)) {
    
    if (!file.exists(feature_list)) {
      message("Feature list file not found: ", feature_list)
      message("Proceeding without feature filtering.")
    } else {
      
      feature_ids <- suppressWarnings(readLines(feature_list))
      feature_ids <- trimws(feature_ids)
      
      if (length(feature_ids) == 0) {
        message("Feature list is empty. Proceeding without filtering.")
      } else {
        
        keep_idx <- rownames(data) %in% feature_ids
        
        n_keep <- sum(keep_idx)
        n_total <- nrow(data)
        
        pct_keep <- round(100 * n_keep / n_total, 1)
        
        message(sprintf(
          "Feature list filter: retained %d of %d features (%.1f%%)",
          n_keep, n_total, pct_keep
        ))
        
        if (n_keep == 0) {
          message("WARNING: No features matched filter list. Proceeding without filtering.")
        } else {
          data <- data[keep_idx, , drop = FALSE]
        }
      }
    }
    
  } else {
    message("No feature list provided. Skipping feature filtering.")
  }
  
  # ----------------------------
  # 2. Align with log2fc
  # ----------------------------
  if (!is.null(water$log2fc)) {
    
    log2fc <- water$log2fc
    
    common <- intersect(rownames(data), rownames(log2fc))
    
    if (length(common) == 0) {
      stop("No overlapping features between filtered data and log2fc")
    }
    
    message(sprintf(
      "Aligned with log2FC data: retained %d of %d features",
      length(common), nrow(data)
    ))
    
    data <- data[common, , drop = FALSE]
    log2fc <- log2fc[common, , drop = FALSE]
    
    # enforce same order
    log2fc <- log2fc[rownames(data), , drop = FALSE]
    
  } else {
    stop("log2fc data required for differential abundance binning")
  }
  
  # ----------------------------
  # 3. Average replicates (AFTER filtering)
  # ----------------------------
  if (average_replicates) {
    
    tp_levels <- water$timepoint_levels
    tp_labels <- water$timepoint_labels
    
    data_avg <- vapply(tp_levels, function(tp) {
      cols <- which(tp_labels == tp)
      rowMeans(data[, cols, drop = FALSE])
    }, numeric(nrow(data)))
    
    data_avg <- as.data.frame(data_avg)
    rownames(data_avg) <- rownames(data)
    colnames(data_avg) <- tp_levels
    
  } else {
    data_avg <- data
  }
  
  # ----------------------------
  # 4. Differential abundance binning
  # ----------------------------
  is_da <- apply(abs(log2fc) >= log2fc_threshold, 1, any)
  is_da <- as.logical(is_da)
  names(is_da) <- rownames(data_avg)
  
  is_non_da <- !is_da
  names(is_non_da) <- rownames(data_avg)
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$data_avg <- data_avg
  water$log2fc <- log2fc
  water$is_da <- is_da
  water$is_non_da <- is_non_da
  
  message("Differential abundance (DA) classification: ")
  message("DA features identified: ", sum(water$is_da))
  message("NonDA features identified: ", sum(water$is_non_da))
  
  n_final <- nrow(data_avg)
  
  water$history[[length(water$history) + 1]] <- list(
    step = "filter",
    n_start = n_start,
    n_final = n_final,
    n_removed = n_start - n_final,
    n_da = sum(is_da),
    log2fc_threshold = log2fc_threshold,
    timestamp = Sys.time()
  )
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (water$save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    # --------------------------
    # Announce output directory
    # --------------------------
    filter_path <- water_get_output_path(water, "filter", "")
    
    message("\nSaving outputs to:")
    message(filter_path)
    message("")
    
    water_save_output(
      water$data_avg,
      "water_data_avg",
      water,
      step = "filter"
    )
    
    water_save_output(
      data.frame(is_da = water$is_da),
      "water_is_da",
      water,
      step = "filter"
    )
    
    water_save_output(
      data.frame(is_non_da = water$is_non_da),
      "water_is_non_da",
      water,
      step = "filter"
    )
  }
  
  return(water)
}