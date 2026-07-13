water_filter_expression <- function(water, tpm_threshold = 0.5) {
  
  if (is.null(water$data)) {
    stop("No raw data found")
  }
  
  data <- water$data
  
  # ----------------------------
  # Timepoint-level expression filtering
  # ----------------------------
  tp_levels <- water$timepoint_levels
  tp_labels <- water$timepoint_labels
  
  expr_matrix <- sapply(tp_levels, function(tp) {
    cols <- which(tp_labels == tp)
    rowMeans(data[, cols, drop = FALSE]) >= tpm_threshold
  })
  
  # Keep genes expressed in at least one timepoint
  keep <- apply(expr_matrix, 1, any)
  
  data_filt <- data[keep, , drop = FALSE]
  
  if (nrow(data_filt) == 0) {
    stop("No features remain after TPM filtering")
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$data_filt <- data_filt
  water$feature_ids <- rownames(data_filt)
  
  water$history[[length(water$history) + 1]] <- list(
    step = "filter_expression",
    tpm_threshold = tpm_threshold,
    method = "timepoint_mean",
    n_retained = nrow(data_filt),
    timestamp = Sys.time()
  )
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (water$save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    expr_path <- water_get_output_path(water, "filter_expression", "")
    
    message("\nSaving outputs to:")
    message(expr_path)
    message("")
    
    water_save_output(
      water$data_filt,
      "water_data_filtered_expression",
      water,
      step = "filter_expression"
    )
  }
  
  pct <- round(100 * sum(keep) / length(keep), 1)
  
  message(sprintf(
    "Expression filter (TPM ≥ %.2f, timepoint mean): retained %d of %d features (%.1f%%)",
    tpm_threshold,
    sum(keep),
    length(keep),
    pct
  ))
  
  return(water)
}