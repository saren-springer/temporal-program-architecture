water_save_output <- function(object, filename, water, step) {
  
  # ----------------------------
  # Validate inputs
  # ----------------------------
  if (is.null(water$output_dir)) {
    stop("Invalid WATER object: missing output_dir")
  }
  
  if (missing(step)) {
    stop("Argument 'step' must be provided")
  }
  
  # sanitize filename
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", filename)
  
  # check Excel availability once
  has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
  
  # ----------------------------
  # Define output directories
  # ----------------------------
  rds_dir   <- water_get_output_path(water, step, "rds")
  table_dir <- water_get_output_path(water, step, "tables")
  
  rds_path  <- normalizePath(file.path(rds_dir,   paste0(safe_name, ".rds")),  winslash = "/", mustWork = FALSE)
  xlsx_path <- normalizePath(file.path(table_dir, paste0(safe_name, ".xlsx")), winslash = "/", mustWork = FALSE)
  
  # ----------------------------
  # Save RDS
  # ----------------------------
  saveRDS(object, file = rds_path)
  
  # ----------------------------
  # Save Excel (if applicable)
  # ----------------------------
  if ((is.data.frame(object) || is.matrix(object)) && has_openxlsx) {
    
    df <- as.data.frame(object)
    
    df <- data.frame(
      Feature_ID = if (!is.null(rownames(df))) rownames(df) else seq_len(nrow(df)),
      df,
      row.names = NULL
    )
    
    if (nrow(df) > 1e6) {
      warning("Object too large for Excel export; saving RDS only")
    } else {
      openxlsx::write.xlsx(
        df,
        file = xlsx_path,
        rowNames = FALSE
      )
    }
    
  } else if (!has_openxlsx) {
    warning("Package 'openxlsx' not installed; skipping Excel export")
  }
  
  # ----------------------------
  # Message
  # ----------------------------
  message(
    "Saved: ", safe_name,
    " → ", step,
    " (.rds", 
    if (has_openxlsx) " + .xlsx" else "",
    ")"
  )
  
  invisible(TRUE)
}