#' Run PCA on normalized WATER data
#'
#' @param water WATER object
#' @param center Logical; center data (default TRUE)
#' @param scale Logical; scale data (default TRUE)
#'
#' @return Updated WATER object with PCA results
#' @export

water_PCA <- function(water, plot = TRUE, center = TRUE, scale = TRUE) {
  
  if (is.null(water$data_norm)) {
    stop("Run water_normalize() before PCA")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for PCA plotting")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for PCA metadata merging")
  }
  
  if (is.null(water$data_avg_norm)) {
    stop("Run water_normalize() before PCA")
  }
  
  common <- intersect(rownames(water$data_norm), rownames(water$data_avg))
  data <- water$data_norm[common, , drop = FALSE]
  
  # ----------------------------
  # Remove zero-variance features
  # ----------------------------
  row_vars <- apply(data, 1, var)
  data_filt <- data[row_vars > 0, , drop = FALSE]
  if (nrow(data_filt) == 0) {
    stop("No features remain after variance filtering for PCA")
  }
  
  # ----------------------------
  # PCA
  # ----------------------------
  pca_data <- prcomp(
    t(data_filt),
    center = center,
    scale. = scale
  )
  
  # ----------------------------
  # Scores
  # ----------------------------
  pca_scores <- as.data.frame(pca_data$x)
  pca_scores$sample <- rownames(pca_scores)
  
  # ----------------------------
  # Attach metadata
  # ----------------------------
  pca_scores <- dplyr::left_join(
    pca_scores,
    water$timepoint_map,
    by = "sample"
  )
  
  pca_scores$timepoint <- factor(
    pca_scores$timepoint,
    levels = water$timepoint_levels
  )
  
  # ----------------------------
  # Variance explained
  # ----------------------------
  pc1_var <- round(100 * summary(pca_data)$importance[2, 1], 1)
  pc2_var <- round(100 * summary(pca_data)$importance[2, 2], 1)
  
  # ----------------------------
  # Plot
  # ----------------------------
  if (plot) {
    
    p <- ggplot2::ggplot(
      pca_scores,
      ggplot2::aes(PC1, PC2, color = timepoint)
    ) +
      ggplot2::geom_point(size = 4) +
      ggplot2::labs(
        title = paste0("WATER PCA (", water$modality, ")"),
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)")
      ) +
      ggplot2::theme_minimal(base_size = 14)
    
    print(p)
    
    if (water$save_outputs) {
      if (!requireNamespace("svglite", quietly = TRUE)) {
        warning("Package 'svglite' not installed; skipping plot save")
      } else {
        plot_dir <- water_get_output_path(water, "PCA", "plots")
        plot_path <- file.path(plot_dir, "water_PCA.svg")
        plot_path <- normalizePath(plot_path, winslash = "/", mustWork = FALSE)
        svglite::svglite(plot_path, width = 8, height = 6)
        print(p)
        dev.off()
      }
    }
    
    water$pca_plot <- p
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$pca <- list(
    model = pca_data,
    scores = pca_scores,
    variance_explained = summary(pca_data)$importance[2, ]
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step = "PCA",
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
    pca_path <- water_get_output_path(water, "PCA", "")
    
    message("\nSaving outputs to:")
    message(pca_path)
    message("")
    
    # --------------------------
    # Save files
    # --------------------------
    water_save_output(
      water$pca$scores,
      "water_pca_scores",
      water,
      step = "PCA"
    )
  }
  
  return(water)
}