#' Initial LOO stability analysis for WATER classic pipeline
#'
#' Performs leave-one-replicate-out stability analysis at both the cluster
#' and trajectory level following initial trajectory assignment.
#'
#' @param water WATER object
#' @param seed Integer; random seed for reproducibility. Default 7.
#'
#' @return Updated WATER object with LOO results
#' @export

water_loo_initial <- function(water, seed = 7) {
  
  set.seed(seed)
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$cluster_classic)) {
    stop("Run water_cluster_classic() before water_loo_initial()")
  }
  
  if (is.null(water$cut_classic)) {
    stop("Run water_cut_classic() before water_loo_initial()")
  }
  
  if (is.null(water$label_classic)) {
    stop("Run water_label_classic() before water_loo_initial()")
  }
  
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' required")
  }
  
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' required")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required")
  }
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' required")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }
  
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' required")
  }
  
  `%>%` <- magrittr::`%>%`
  
  # ----------------------------
  # Retrieve baseline objects
  # ----------------------------
  k_final    <- water$cut_classic$k_final
  tp_levels  <- water$timepoint_levels
  tp_labels  <- water$timepoint_labels
  
  # Features that made it into clustering
  features <- rownames(water$cluster_classic$input_matrix)
  
  # Replicate-level QN data subset to clustered features
  data_rep <- water$data_norm[
    intersect(features, rownames(water$data_norm)), , drop = FALSE
  ]
  
  # Baseline cluster assignments (named by feature_id)
  baseline_cluster_df <- water$cut_classic$cluster_df
  baseline_clusters <- setNames(
    as.integer(as.character(baseline_cluster_df$Cluster)),
    baseline_cluster_df$feature_id
  )
  
  # Baseline trajectory assignments (named by feature_id)
  baseline_traj_df <- water$label_classic$cluster_df
  baseline_traj <- setNames(
    baseline_traj_df$Trajectory,
    baseline_traj_df$feature_id
  )
  
  # Trajectory map: cluster number -> trajectory name
  trajectory_map <- water$label_classic$trajectory_map
  
  # Baseline z-scored matrix
  rowz <- function(M) t(scale(t(M)))
  
  Z_base <- rowz(water$cluster_classic$input_matrix)
  Z_base <- Z_base[rownames(Z_base) %in% names(baseline_clusters), , drop = FALSE]
  
  # ----------------------------
  # Helper: compute cluster centroids
  # ----------------------------
  cluster_centroids <- function(Z, labels) {
    labs <- as.factor(labels[rownames(Z)])
    split_rows <- split(rownames(Z), labs)
    cents <- lapply(split_rows, function(ix) {
      colMeans(Z[ix, , drop = FALSE])
    })
    cents <- do.call(rbind, cents)
    rn_num <- suppressWarnings(as.numeric(rownames(cents)))
    if (!all(is.na(rn_num))) {
      cents <- cents[order(rn_num), , drop = FALSE]
    }
    cents
  }
  
  # ----------------------------
  # Helper: match LOO clusters to baseline via Hungarian algorithm
  # ----------------------------
  match_clusters <- function(Z_b, lab_b, Z_l, lab_l) {
    
    g <- intersect(rownames(Z_b), rownames(Z_l))
    Zb <- Z_b[g, , drop = FALSE]
    Zl <- Z_l[g, , drop = FALSE]
    lb <- lab_b[g]
    ll <- lab_l[g]
    
    Cb <- cluster_centroids(Zb, lb)
    Cl <- cluster_centroids(Zl, ll)
    
    R <- cor(t(Cl), t(Cb), method = "pearson", use = "pairwise.complete.obs")
    sol <- clue::solve_LSAP(1 - R)
    
    map <- data.frame(
      LOO_cluster      = rownames(R),
      Baseline_cluster = colnames(R)[sol],
      corr             = R[cbind(seq_len(nrow(R)), sol)],
      stringsAsFactors = FALSE
    )
    
    lut <- setNames(map$Baseline_cluster, map$LOO_cluster)
    ll_mapped <- unname(lut[as.character(ll)])
    
    list(mapping = map, ll_mapped = ll_mapped, genes = g)
  }
  
  # ----------------------------
  # Helper: cluster LOO data
  # ----------------------------
  cluster_k <- function(M_tp, k) {
    Z <- rowz(M_tp)
    keep <- is.finite(rowSums(Z)) & apply(Z, 1, sd) > 0
    Z <- Z[keep, , drop = FALSE]
    hc <- hclust(dist(Z), method = "ward.D2")
    lab <- cutree(hc, k = k)
    list(labels = lab, Z = Z, keep = keep)
  }
  
  # ----------------------------
  # Run LOO
  # ----------------------------
  rep_cols <- colnames(data_rep)
  loo_summaries <- list()
  
  message(sprintf("  Running %d LOO iterations", length(rep_cols)))
  
  plot_dir  <- water_get_output_path(water, "loo_initial", "plots")
  table_dir <- water_get_output_path(water, "loo_initial", "tables")
  
  for (drop_col in rep_cols) {
    
    drop_tp <- tp_labels[which(colnames(data_rep) == drop_col)]
    message(sprintf("    Dropping sample: %s (%s)", drop_col, drop_tp))
    
    # Recompute timepoint means without dropped replicate
    M_loo <- do.call(cbind, lapply(tp_levels, function(tp) {
      cols <- which(tp_labels == tp)
      cols <- cols[colnames(data_rep)[cols] != drop_col]
      if (length(cols) < 1) return(NULL)
      rowMeans(data_rep[, cols, drop = FALSE], na.rm = TRUE)
    }))
    
    if (is.null(M_loo) || ncol(M_loo) < length(tp_levels)) {
      message("  Skipping ", drop_col, " (timepoint would be empty)")
      next
    }
    
    colnames(M_loo) <- tp_levels
    rownames(M_loo) <- rownames(data_rep)
    
    # Subset to baseline features
    M_loo <- M_loo[rownames(M_loo) %in% names(baseline_clusters), , drop = FALSE]
    
    # Cluster LOO data
    loo <- cluster_k(M_loo, k = k_final)
    
    # Restrict to shared features
    g_inter <- intersect(names(baseline_clusters), names(loo$labels))
    
    if (length(g_inter) < 200) {
      message("  Overlap too small (n=", length(g_inter), "). Skipping.")
      next
    }
    
    lb <- baseline_clusters[g_inter]
    ll <- loo$labels[g_inter]
    
    # ARI before mapping
    ari_raw <- mclust::adjustedRandIndex(lb, ll)
    
    # Match LOO clusters to baseline
    mres <- match_clusters(
      Z_base[g_inter, , drop = FALSE], lb,
      loo$Z[g_inter, , drop = FALSE],  ll
    )
    
    ll_mapped <- as.integer(mres$ll_mapped)
    ari_mapped <- mclust::adjustedRandIndex(lb, ll_mapped)
    
    # Cluster retention
    same_cluster <- lb == ll_mapped
    retention_cluster <- mean(same_cluster, na.rm = TRUE)
    
    # Trajectory retention
    traj_base <- baseline_traj[g_inter]
    traj_loo  <- trajectory_map[as.character(ll_mapped)]
    traj_loo[is.na(traj_loo)] <- "Unassigned"
    
    same_traj <- traj_base == traj_loo
    retention_traj <- mean(same_traj, na.rm = TRUE)
    
    # Per-trajectory retention
    per_traj <- data.frame(
      Trajectory = traj_base,
      Kept       = same_traj,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::count(Trajectory, Kept) %>%
      tidyr::pivot_wider(
        names_from  = Kept,
        values_from = n,
        values_fill = 0
      ) %>%
      dplyr::rename(
        n_same   = `TRUE`,
        n_change = `FALSE`
      ) %>%
      dplyr::mutate(
        total     = n_same + n_change,
        retention = ifelse(total > 0, n_same / total, NA_real_)
      ) %>%
      dplyr::arrange(dplyr::desc(retention))
    
    # Save per-run outputs
    if (water$save_outputs) {
      
      safe_name <- gsub("[^A-Za-z0-9]+", "_", drop_col)
      run_dir <- file.path(table_dir, paste0("LOO_", safe_name))
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      
      write.csv(
        mres$mapping,
        file.path(run_dir, "cluster_mapping.csv"),
        row.names = FALSE
      )
      
      write.csv(
        data.frame(
          feature_id         = g_inter,
          BaselineCluster    = lb,
          LOOCluster_raw     = as.integer(ll),
          LOOCluster_mapped  = ll_mapped,
          BaselineTrajectory = traj_base,
          LOOTrajectory      = traj_loo,
          stringsAsFactors   = FALSE
        ),
        file.path(run_dir, "labels_baseline_vs_LOO.csv"),
        row.names = FALSE
      )
      
      write.csv(
        per_traj,
        file.path(run_dir, "per_trajectory_retention.csv"),
        row.names = FALSE
      )
    }
    
    loo_summaries[[drop_col]] <- data.frame(
      ReplicateDropped  = drop_col,
      DroppedTimepoint  = drop_tp,
      OverlapFeatures   = length(g_inter),
      ARI_raw           = ari_raw,
      ARI_mapped        = ari_mapped,
      Retention_cluster = retention_cluster,
      Retention_traj    = retention_traj,
      stringsAsFactors  = FALSE
    )
  }
  
  # ----------------------------
  # Summarise across iterations
  # ----------------------------
  summary_df <- do.call(rbind, loo_summaries)
  summary_df <- summary_df[order(summary_df$ARI_mapped, decreasing = TRUE), ]
  
  message("\n  Summary:")
  
  message(sprintf("    Mean ARI: %.3f",
                  mean(summary_df$ARI_mapped, na.rm = TRUE)))
  
  message(sprintf("    Cluster retention: %.1f%%",
                  100 * mean(summary_df$Retention_cluster, na.rm = TRUE)))
  
  message(sprintf("    Trajectory retention: %.1f%%",
                  100 * mean(summary_df$Retention_traj, na.rm = TRUE)))
  
  if (water$save_outputs) {
    water_save_output(
      summary_df,
      "water_loo_initial_summary",
      water,
      step = "loo_initial"
    )
  }
  
  # ----------------------------
  # Plot
  # ----------------------------
  summary_df$Cluster_pct    <- summary_df$Retention_cluster * 100
  summary_df$Trajectory_pct <- summary_df$Retention_traj    * 100
  
  loo_long <- summary_df %>%
    dplyr::select(ReplicateDropped, Cluster_pct, Trajectory_pct) %>%
    tidyr::pivot_longer(
      cols      = c(Cluster_pct, Trajectory_pct),
      names_to  = "Metric",
      values_to = "Percent"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        Cluster_pct    = paste0("Cluster retention (k = ", k_final, ")"),
        Trajectory_pct = "Trajectory retention"
      )
    )
  
  avg_lines <- loo_long %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(mean_pct = mean(Percent, na.rm = TRUE), .groups = "drop")
  
  loo_plot <- ggplot2::ggplot(
    loo_long,
    ggplot2::aes(x = ReplicateDropped, y = Percent, fill = Metric)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width    = 0.65
    ) +
    ggplot2::geom_hline(
      data     = avg_lines,
      ggplot2::aes(yintercept = mean_pct),
      linetype = "dashed",
      linewidth = 0.7,
      color    = "grey40"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "grey80",
        "grey40"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = "Initial LOO stability analysis",
      y     = "Retention (%)",
      x     = NULL,
      fill  = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position    = "top"
    )
  
  print(loo_plot)
  
  if (water$save_outputs) {
    ggplot2::ggsave(
      file.path(plot_dir, "loo_initial_retention.png"),
      loo_plot,
      width  = 10,
      height = 6,
      dpi    = 300
    )
    ggplot2::ggsave(
      file.path(plot_dir, "loo_initial_retention.pdf"),
      loo_plot,
      width  = 10,
      height = 6
    )
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$loo_initial <- list(
    summary    = summary_df,
    loo_plot   = loo_plot,
    k_final    = k_final,
    n_features = length(features)
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step               = "loo_initial",
    k_final            = k_final,
    n_iterations       = nrow(summary_df),
    mean_ari_mapped    = mean(summary_df$ARI_mapped,        na.rm = TRUE),
    mean_retention_cluster  = mean(summary_df$Retention_cluster, na.rm = TRUE),
    mean_retention_traj     = mean(summary_df$Retention_traj,    na.rm = TRUE),
    timestamp          = Sys.time()
  )
  
  return(water)
}