water_get_output_path <- function(water, step, type) {
  
  # ----------------------------
  # Validate inputs
  # ----------------------------
  if (is.null(water$output_dir)) {
    stop("water$output_dir is NULL — initialize WATER object properly")
  }
  
  if (missing(step) || missing(type)) {
    stop("Both 'step' and 'type' must be provided")
  }
  
  # sanitize inputs (avoid weird folder names)
  step <- gsub("[^A-Za-z0-9_\\-]", "_", step)
  type <- gsub("[^A-Za-z0-9_\\-]", "_", type)
  
  # ----------------------------
  # Build path
  # ----------------------------
  path <- file.path(water$output_dir, step, type)
  
  # normalize (important for Windows / pipelines)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  
  # ----------------------------
  # Create if needed
  # ----------------------------
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  return(path)
}