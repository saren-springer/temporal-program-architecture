#' Initialize a WATER object
#'
#' @param data Data frame or matrix. First column = feature IDs,
#' remaining columns = samples (ordered by timepoint and replicate)
#'
#' @param replicates Named list specifying replicates per timepoint
#' (e.g., list(T1 = 2, T2 = 3, T3 = 3))
#'
#' @param log2fc_data Optional data frame/matrix with log2FC values
#'
#' @param modality Character string describing data type
#' (e.g., "RNA", "CUT&RUN", "scRNA", "ATAC")
#'
#' @return WATER object
#' @export

water_initialize <- function(
    data,
    replicates,
    log2fc_data = NULL,
    modality = "unknown",
    output_dir = NULL,
    save_outputs = TRUE
) {
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  
  # ----------------------------
  # Input formatting
  # ----------------------------
  if (is.data.frame(data)) {
    feature_ids <- data[[1]]
    data_matrix <- as.matrix(data[, -1])
  } else if (is.matrix(data)) {
    feature_ids <- rownames(data)
    data_matrix <- data
  } else {
    stop("data must be a data.frame or matrix")
  }
  
  if (is.null(feature_ids)) {
    stop("Feature IDs must be provided (first column or rownames)")
  }
  
  rownames(data_matrix) <- feature_ids
  
  # ensure column names exist
  if (is.null(colnames(data_matrix))) {
    colnames(data_matrix) <- paste0("Sample_", seq_len(ncol(data_matrix)))
    message("No sample names detected — assigning generic names")
  }
  
  n_samples <- ncol(data_matrix)
  
  # ----------------------------
  # Replicate structure (REQUIRED)
  # ----------------------------
  if (is.null(replicates)) {
    stop(
      "Replicate structure must be provided.\n\n",
      "Example:\n",
      "replicates = list(T1 = 2, T2 = 3, T3 = 3)\n\n",
      "Columns must be ordered accordingly."
    )
  }
  
  if (is.null(names(replicates))) {
    stop("replicates must be a named list (e.g., list(T1=2, T2=3))")
  }
  
  timepoint_labels <- unlist(
    mapply(function(tp, n) rep(tp, n),
           names(replicates),
           replicates,
           SIMPLIFY = FALSE),
    use.names = FALSE
  )
  
  if (length(timepoint_labels) != n_samples) {
    stop("Number of samples does not match total replicates specified")
  }
  
  # ----------------------------
  # Timepoint metadata
  # ----------------------------
  timepoint_levels <- names(replicates)
  
  replicate_index <- as.integer(ave(
    timepoint_labels,
    timepoint_labels,
    FUN = seq_along
  ))
  
  timepoint_map <- data.frame(
    sample = colnames(data_matrix),
    timepoint = timepoint_labels,
    replicate = replicate_index,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------
  # Handle log2FC (optional)
  # ----------------------------
  log2fc_matrix <- NULL
  
  if (!is.null(log2fc_data)) {
    
    if (is.data.frame(log2fc_data)) {
      log2fc_matrix <- as.matrix(log2fc_data[, -1])
      rownames(log2fc_matrix) <- log2fc_data[[1]]
    } else if (is.matrix(log2fc_data)) {
      log2fc_matrix <- log2fc_data
    } else {
      stop("log2fc_data must be a data.frame or matrix")
    }
    
    common_features <- intersect(
      rownames(data_matrix),
      rownames(log2fc_matrix)
    )
    
    if (length(common_features) == 0) {
      stop("No overlapping feature IDs between data and log2fc_data")
    }
    
    data_matrix <- data_matrix[common_features, , drop = FALSE]
    log2fc_matrix <- log2fc_matrix[common_features, , drop = FALSE]
    log2fc_matrix <- log2fc_matrix[rownames(data_matrix), , drop = FALSE]
  }
  
  # ----------------------------
  # Create WATER object
  # ----------------------------
  water <- list(
    data = data_matrix,
    log2fc = log2fc_matrix,
    
    feature_ids = rownames(data_matrix),
    
    timepoint_labels = timepoint_labels,
    timepoint_levels = timepoint_levels,
    replicate_index = replicate_index,
    timepoint_map = timepoint_map,
    
    modality = modality,
    output_dir = output_dir,
    save_outputs = save_outputs,
    
    clusters = NULL,
    trajectories = NULL,
    history = list()
  )
  
  class(water) <- "WATER"
  
  water$history[[1]] <- list(
    step = "initialize",
    modality = modality,
    n_features = nrow(data_matrix),
    n_samples = ncol(data_matrix),
    n_timepoints = length(timepoint_levels),
    timestamp = Sys.time()
  )
  
  tp_table <- table(timepoint_labels)
  
  message("\nTimepoint structure:")
  for (i in seq_along(tp_table)) {
    message(sprintf("  %s: %d replicates", names(tp_table)[i], tp_table[i]))
  }
  message(sprintf("  Total samples: %d", sum(tp_table)))
  
  message("\nSample mapping:")
  for (i in seq_len(nrow(timepoint_map))) {
    message(sprintf("  %s -> %s (rep %d)",
                    timepoint_map$sample[i],
                    timepoint_map$timepoint[i],
                    timepoint_map$replicate[i]))
  }
  
  # Save to file (only if enabled)
  if (save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    table_dir <- water_get_output_path(water, "initialize", "tables")
    
    message("\nSaving outputs:")
    message(sprintf("  Tables: %s", table_dir))
    
    water_save_output(
      water$data,
      "water_raw_data",
      water,
      step = "initialize"
    )
    
    water_save_output(
      timepoint_map,
      "water_sample_mapping",
      water,
      step = "initialize"
    )
    
    note_path <- file.path(table_dir, "water_sample_mapping.txt")
    writeLines(capture.output(print(timepoint_map)), con = note_path)
  }
  
  message("\nInitialization complete.")
  message("Output directory: ", output_dir)
  message(sprintf("Features: %d | Samples: %d | Timepoints: %d",
                  nrow(data_matrix),
                  ncol(data_matrix),
                  length(timepoint_levels)))
  
  return(water)
}