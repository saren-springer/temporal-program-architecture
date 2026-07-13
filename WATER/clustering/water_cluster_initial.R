#' Run classic hierarchical clustering for WATER trajectory analysis
#'
#' Performs silhouette analysis and generates heatmaps for candidate k values.
#' Inspect the silhouette plot and heatmaps before choosing k_final with
#' water_cut_clusters().
#'
#' @param water WATER object
#' @param k_max Integer; maximum k for silhouette analysis. Defaults to
#'   n_timepoints * 7 + 2.
#'
#' @return Updated WATER object with clustering results
#' @export

water_cluster_initial <- function(water, k_max = NULL) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$data_norm)) {
    stop("Run water_normalize() before clustering")
  }
  
  if (is.null(water$is_da)) {
    stop("DA classification not found. Run water_filter_features() first.")
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required")
  }
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' required")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }
  
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("Package 'factoextra' required")
  }
  
  # ----------------------------
  # Step 1: Build matrix
  # ----------------------------
  
  # Subset data_norm to feature-list filtered genes, then to DA genes
  common <- intersect(rownames(water$data_norm), rownames(water$data_avg))
  data_rep <- as.matrix(water$data_norm[common, , drop = FALSE])
  data_rep <- data_rep[names(water$is_da)[water$is_da], , drop = FALSE]
  
  # Average QN replicates by timepoint
  tp_levels <- water$timepoint_levels
  tp_labels <- water$timepoint_labels
  
  data_avg <- vapply(tp_levels, function(tp) {
    cols <- which(tp_labels == tp)
    rowMeans(data_rep[, cols, drop = FALSE])
  }, numeric(nrow(data_rep)))
  
  data_avg <- as.data.frame(data_avg)
  rownames(data_avg) <- rownames(data_rep)
  colnames(data_avg) <- tp_levels
  
  # ----------------------------
  # Step 2: Clean matrix
  # ----------------------------
  keep <-
    apply(data_avg, 1, function(x) all(is.finite(x))) &
    apply(data_avg, 1, sd) > 0
  
  data_avg <- data_avg[keep, , drop = FALSE]
  
  if (nrow(data_avg) < 2) {
    stop("Not enough features after filtering")
  }
  
  message("Clustering ", nrow(data_avg), " DA features")
  
  # ----------------------------
  # Step 3: Scale + distance
  # ----------------------------
  scaled_data <- scale(data_avg)
  diss <- dist(scaled_data)
  
  # Compute and store hclust object for use in water_cut_classic
  hc <- hclust(diss, method = "ward.D2")
  
  n_timepoints <- ncol(data_avg)
  
  if (is.null(k_max)) {
    k_max <- n_timepoints * 7 + 2
  }
  
  # ----------------------------
  # Step 4: Silhouette analysis
  # ----------------------------
  
  message("Running silhouette analysis across k = 2 to ", k_max, "...")
  message("This may take a few minutes for large datasets.")
  
  sil_plot <- factoextra::fviz_nbclust(
    scaled_data,
    FUN = factoextra::hcut,
    method = "silhouette",
    k.max = k_max,
    diss = diss
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "WATER Silhouette Analysis",
      subtitle = paste0("k.max = ", k_max)
    )
  
  print(sil_plot)
  
  plot_dir <- water_get_output_path(water, "cluster", "plots")
  
  if (water$save_outputs) {
    png(file.path(plot_dir, "silhouette_plot.png"),
        width = 1600, height = 1200, res = 150)
    print(sil_plot)
    dev.off()
  }
  
  # ----------------------------
  # Step 5: Extract silhouette values
  # ----------------------------
  sil_df <- sil_plot$data
  
  if (!all(c("clusters", "y") %in% colnames(sil_df))) {
    stop("Unexpected silhouette plot structure from factoextra")
  }
  
  sil_df <- data.frame(
    k = as.integer(as.character(sil_df$clusters)),
    avg_silhouette = sil_df$y
  )
  
  sil_vals <- sil_df$avg_silhouette
  ks <- sil_df$k
  
  # ----------------------------
  # Step 6: Select candidate k values
  # ----------------------------
  
  # All k with silhouette > 0.5
  k_good <- ks[sil_vals > 0.5]
  
  # Tier 1: k values just before a drop >= 0.1
  drops <- diff(sil_vals)
  drop_positions <- which(drops <= -0.075)
  k_tier1 <- ks[drop_positions]
  k_tier1 <- k_tier1[k_tier1 %in% k_good]
  
  # Tier 2: last k before silhouette drops below 0.5
  k_tier2 <- max(ks[sil_vals >= 0.5], na.rm = TRUE)
  
  # Combined recommendation
  k_recommended <- sort(unique(c(k_tier1, k_tier2)))
  
  message("Generating heatmaps for inspection...\n")
  
  # ----------------------------
  # Step 7: Heatmaps for recommended k values
  # ----------------------------
  
  color_palette <- colorRampPalette(
    c("#408493", "#94D2BD", "#E9D8A6", "#EE9C21", "#B73A27")
  )
  
  save_pheatmap_png <- function(x, filename) {
    png(filename, width = 1400, height = 3000, res = 150)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  if (length(k_recommended) > 0) {
    
    # Save individual heatmaps and collect gtables for side-by-side display
    grob_list <- list()
    
    for (k in k_recommended) {
      message("Generating heatmap for k = ", k)
      
      ph <- pheatmap::pheatmap(
        data_avg,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = color_palette(100),
        cutree_rows = k,
        show_rownames = FALSE,
        main = paste0("k = ", k, " (sil = ",
                      round(sil_df$avg_silhouette[sil_df$k == k], 2), ")"),
        silent = TRUE
      )
      
      grob_list[[as.character(k)]] <- ph$gtable
      
      if (water$save_outputs) {
        save_pheatmap_png(
          ph,
          file.path(plot_dir, paste0("heatmap_k_", k, ".png"))
        )
      }
    }
    
    # Display side by side in plots tab with title
    n <- length(grob_list)
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 2, ncol = n,
          heights = grid::unit(c(0.05, 0.95), "npc")
        )
      )
    )
    
    # Title row
    grid::pushViewport(
      grid::viewport(layout.pos.row = 1, layout.pos.col = 1:n)
    )
    grid::grid.text(
      "Recommended k values for clustering",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    grid::upViewport()
    
    # Heatmap row
    for (i in seq_along(grob_list)) {
      grid::pushViewport(
        grid::viewport(layout.pos.row = 2, layout.pos.col = i)
      )
      grid::grid.draw(grob_list[[i]])
      grid::upViewport()
    }
    grid::upViewport()  # pop back to root
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$cluster_classic <- list(
    input_matrix = data_avg,
    scaled_data = scaled_data,
    hclust = hc,
    silhouette_plot = sil_plot,
    silhouette = sil_df,
    k_good = k_good,
    k_recommended = k_recommended
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step = "cluster_classic",
    k_max = k_max,
    n_features = nrow(data_avg),
    drop_threshold = 0.075,
    k_recommended = k_recommended,
    timestamp = Sys.time()
  )
  
  return(water)
}