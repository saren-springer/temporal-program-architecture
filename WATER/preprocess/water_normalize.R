#' Quantile normalize WATER data (replicate level)
#'
#' @param water WATER object
#'
#' @return Updated WATER object
#' @export

water_normalize <- function(water) {
  
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' is required")
  }
  
  # ----------------------------
  # 1. Normalize replicate data
  # ----------------------------
  if (is.null(water$data_filt)) {
    stop("Run water_filter_expression() first")
  }
  
  data_rep <- water$data_filt
  
  data_rep_qn <- preprocessCore::normalize.quantiles(
    data.matrix(data_rep)
  )
  
  rownames(data_rep_qn) <- rownames(data_rep)
  colnames(data_rep_qn) <- colnames(data_rep)
  
  data_rep_qn <- as.data.frame(data_rep_qn)
  
  # Remove zero-variance
  row_vars <- apply(data_rep_qn, 1, var)
  data_rep_qn <- data_rep_qn[row_vars > 0, , drop = FALSE]
  
  water$data_norm <- data_rep_qn
  
  # ----------------------------
  # 2. Average the QN replicate data by timepoint
  # ----------------------------
  tp_levels <- water$timepoint_levels
  tp_labels <- water$timepoint_labels
  
  data_avg_qn <- vapply(tp_levels, function(tp) {
    cols <- which(tp_labels == tp)
    rowMeans(data_rep_qn[, cols, drop = FALSE])
  }, numeric(nrow(data_rep_qn)))
  
  data_avg_qn <- as.data.frame(data_avg_qn)
  rownames(data_avg_qn) <- rownames(data_rep_qn)
  colnames(data_avg_qn) <- tp_levels
  
  water$data_avg_norm <- data_avg_qn
  
  # ----------------------------
  # 3. Subset DA / non-DA from averaged QN data
  # ----------------------------
  common <- intersect(rownames(data_avg_qn), names(water$is_da))
  
  water$data_da     <- data_avg_qn[common[water$is_da[common]],  , drop = FALSE]
  water$data_non_da <- data_avg_qn[common[water$is_non_da[common]], , drop = FALSE]
  
  # ----------------------------
  # History
  # ----------------------------
  water$history[[length(water$history) + 1]] <- list(
    step = "normalize",
    n_rep_features = nrow(water$data_norm),
    n_avg_features = nrow(water$data_avg_norm),
    avg_method = "mean_of_QN_replicates",
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
    norm_path <- water_get_output_path(water, "normalize", "")
    
    message("\nSaving outputs to:")
    message(norm_path)
    message("")
    
    # --------------------------
    # Save files
    # --------------------------
    water_save_output(
      water$data_norm,
      "water_data_normalized_replicates",
      water,
      step = "normalize"
    )
    
    water_save_output(
      water$data_avg_norm,
      "water_data_normalized_avg",
      water,
      step = "normalize"
    )
  }
  
  return(water)
}