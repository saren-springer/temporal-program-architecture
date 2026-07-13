#' Cut hierarchical clustering tree at chosen k
#'
#' Cuts the clustering tree from water_cluster_classic() at a chosen k,
#' generates faceted trajectory plots for inspection, and prepares the
#' water object for trajectory labeling via water_label_classic().
#'
#' @param water WATER object
#' @param k_final Integer; number of clusters to cut the tree at.
#'   Defaults to the largest recommended k from water_cluster_classic().
#'
#' @return Updated WATER object with cluster assignments
#' @export

water_cut <- function(water, k_final = NULL) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$cluster_classic)) {
    stop("Run water_cluster_initial() before water_cut()")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' required")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required")
  }
  
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' required")
  }
  
  `%>%` <- magrittr::`%>%`
  
  # ----------------------------
  # Default k_final
  # ----------------------------
  if (is.null(k_final)) {
    k_final <- max(water$cluster_classic$k_recommended)
    message("No k_final specified. Defaulting to k = ", k_final,
            " (largest recommended value).")
    message("To use a different k, run water_cut_classic(water, k_final = X)")
  }
  
  k_final <- as.integer(k_final)
  
  if (k_final < 2) {
    stop("k_final must be >= 2")
  }
  
  message("Cutting tree at k = ", k_final)
  
  # ----------------------------
  # Cut tree
  # ----------------------------
  color_palette <- colorRampPalette(
    c("#408493", "#94D2BD", "#E9D8A6", "#EE9C21", "#B73A27")
  )
  
  ph <- pheatmap::pheatmap(
    water$cluster_classic$input_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = color_palette(100),
    cutree_rows = k_final,
    show_rownames = FALSE,
    silent = TRUE
  )
  
  clusters <- cutree(ph$tree_row, k = k_final)
  
  cluster_df <- data.frame(
    feature_id = names(clusters),
    Cluster    = factor(clusters),
    stringsAsFactors = FALSE
  )
  
  # Cluster size summary
  cluster_size_df <- cluster_df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop")
  
  message("Cluster sizes:")
  print(as.data.frame(cluster_size_df), row.names = FALSE)
  
  # ----------------------------
  # Row-scale averaged data
  # ----------------------------
  data_avg <- water$cluster_classic$input_matrix
  data_z <- t(scale(t(data_avg)))
  data_z <- as.data.frame(data_z)
  data_z$Cluster  <- factor(clusters[rownames(data_z)])
  data_z$feature_id <- rownames(data_z)
  
  # ----------------------------
  # Faceted trajectory plot
  # ----------------------------
  tp_levels <- water$timepoint_levels
  
  data_long <- data_z %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(tp_levels),
      names_to  = "Timepoint",
      values_to = "Z_Score"
    ) %>%
    dplyr::mutate(
      Timepoint = factor(Timepoint, levels = tp_levels)
    )
  
  n_clusters <- k_final
  n_col      <- min(5, n_clusters)
  
  traj_plot <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(x = Timepoint, y = Z_Score)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = feature_id),
      color = "grey",
      alpha = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = Cluster),
      fun      = mean,
      geom     = "line",
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ Cluster, ncol = n_col) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0("WATER Clusters (k = ", k_final, ", ",
                     nrow(data_avg), " features)"),
      x     = "",
      y     = "Row Z-score"
    ) +
    ggplot2::theme(
      legend.position       = "none",
      panel.grid.major.y    = ggplot2::element_blank(),
      panel.grid.minor      = ggplot2::element_blank(),
      strip.text            = ggplot2::element_text(face = "bold")
    )
  
  print(traj_plot)
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (water$save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    # --------------------------
    # Resolve paths
    # --------------------------
    plot_dir  <- water_get_output_path(water, "cut_classic", "plots")
    table_dir <- water_get_output_path(water, "cut_classic", "tables")
    
    # --------------------------
    # Announce output directories
    # --------------------------
    message("\nSaving outputs to:")
    message("  Plots : ", plot_dir)
    message("  Tables: ", table_dir)
    message("")
    
    # --------------------------
    # Save plots
    # --------------------------
    ggplot2::ggsave(
      file.path(plot_dir, paste0("trajectory_clusters_k", k_final, ".png")),
      traj_plot,
      width  = 12,
      height = 6,
      dpi    = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, paste0("trajectory_clusters_k", k_final, ".pdf")),
      traj_plot,
      width  = 12,
      height = 6
    )
    
    # --------------------------
    # Save tables
    # --------------------------
    water_save_output(
      cluster_df,
      "water_cluster_assignments",
      water,
      step = "cut_classic"
    )
    
    water_save_output(
      cluster_size_df,
      "water_cluster_sizes",
      water,
      step = "cut_classic"
    )
    
    water_save_output(
      data_z,
      "water_data_zscored",
      water,
      step = "cut_classic"
    )
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$cut_classic <- list(
    k_final        = k_final,
    cluster_df     = cluster_df,
    cluster_sizes  = cluster_size_df,
    data_z         = data_z,
    trajectory_plot = traj_plot
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step      = "cut_classic",
    k_final   = k_final,
    n_features = nrow(data_avg),
    timestamp = Sys.time()
  )
  
  return(water)
}