#' Secondary correlation analysis for WATER classic pipeline
#'
#' Re-computes correlation between features and FINAL integrated trajectory
#' labels after ambiguous and weak reclustering. This serves as a validation
#' step to assess trajectory coherence and identify any remaining weak fits.
#'
#' @param water WATER object
#' @param r_threshold Numeric; correlation threshold for classification.
#'   Default 0.8.
#'
#' @return Updated WATER object with secondary correlation results
#' @export

water_correlate_secondary <- function(water, r_threshold = 0.8) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$recluster_weak)) {
    stop("Run water_recluster_weak() before water_correlate_secondary()")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required")
  if (!requireNamespace("magrittr", quietly = TRUE)) stop("Package 'magrittr' required")
  
  `%>%` <- magrittr::`%>%`
  
  message("\n  [4.3] Secondary correlation")
  message("----------------------------------------")
  message(sprintf("  r_threshold: %.2f", r_threshold))
  
  # ----------------------------
  # Retrieve data
  # ----------------------------
  data_z <- water$label_classic$data_z
  integrated_df <- water$recluster_weak$integrated_df
  tp_levels <- water$timepoint_levels
  
  # ----------------------------
  # Apply initial trajectory labels
  # ----------------------------
  data_z$Trajectory <- integrated_df$Trajectory[
    match(data_z$feature_id, integrated_df$feature_id)
  ]
  
  data_z$Trajectory <- factor(data_z$Trajectory)
  
  # ----------------------------
  # Build trajectory templates
  # ----------------------------
  tmpl <- sapply(levels(data_z$Trajectory), function(lab) {
    idx <- which(data_z$Trajectory == lab)
    
    if (length(idx) == 0) {
      return(rep(NA, length(tp_levels)))
    }
    
    colMeans(
      as.matrix(data_z[idx, tp_levels, drop = FALSE]),
      na.rm = TRUE
    )
  })
  
  tmpl <- t(tmpl)
  
  # ----------------------------
  # Correlation
  # ----------------------------
  expr_mat <- as.matrix(data_z[, tp_levels, drop = FALSE])
  rownames(expr_mat) <- data_z$feature_id
  
  cor_all <- t(apply(expr_mat, 1, function(v) {
    stats::cor(t(tmpl), v, use = "pairwise.complete.obs")
  }))
  
  rownames(cor_all) <- data_z$feature_id
  colnames(cor_all) <- rownames(tmpl)
  
  # ----------------------------
  # Metrics
  # ----------------------------
  own_idx  <- match(as.character(data_z$Trajectory), rownames(tmpl))
  r_own    <- cor_all[cbind(seq_len(nrow(cor_all)), own_idx)]
  r_best   <- apply(cor_all, 1, max, na.rm = TRUE)
  lab_best <- rownames(tmpl)[apply(cor_all, 1, which.max)]
  
  data_z$r_own    <- r_own
  data_z$r_best   <- r_best
  data_z$lab_best <- lab_best
  
  # ----------------------------
  # Decision logic
  # ----------------------------
  data_z <- data_z %>%
    dplyr::mutate(
      decision = dplyr::case_when(
        r_own >= r_threshold ~ "Stable",
        r_own < r_threshold & r_best >= r_threshold ~ "Reassign_candidate",
        r_best < r_threshold ~ "Low_confidence"
      )
    )
  
  # ----------------------------
  # 🔥 FINAL TRAJECTORY ASSIGNMENT
  # ----------------------------
  data_z <- data_z %>%
    dplyr::mutate(
      Trajectory_final = dplyr::case_when(
        decision == "Stable" ~ as.character(Trajectory),
        decision == "Reassign_candidate" ~ as.character(lab_best),
        decision == "Low_confidence" ~ NA_character_
      )
    )
  
  data_z$Trajectory_final <- factor(data_z$Trajectory_final)
  
  # ----------------------------
  # Decision summary
  # ----------------------------
  decision_counts <- data_z %>%
    dplyr::count(decision) %>%
    dplyr::rename(n_features = n) %>%
    dplyr::mutate(percent = round(100 * n_features / sum(n_features), 1))
  
  message("\n  Decision counts:")
  print(as.data.frame(decision_counts), row.names = FALSE)
  
  # ----------------------------
  # FINAL trajectory summary
  # ----------------------------
  traj_summary <- data_z %>%
    dplyr::filter(!is.na(Trajectory_final)) %>%
    dplyr::group_by(Trajectory_final) %>%
    dplyr::summarise(
      mean_r   = mean(r_own, na.rm = TRUE),
      median_r = median(r_own, na.rm = TRUE),
      n_features = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_features))
  
  message("\n  Final trajectory summary:")
  print(as.data.frame(traj_summary), row.names = FALSE)
  
  # ----------------------------
  # Low confidence
  # ----------------------------
  low_conf_features <- data_z %>%
    dplyr::filter(decision == "Low_confidence") %>%
    dplyr::pull(feature_id)
  
  message(sprintf("\n  Low confidence features: %d", length(low_conf_features)))
  
  # ----------------------------
  # FINAL trajectory plot
  # ----------------------------
  data_long_traj <- data_z %>%
    dplyr::filter(!is.na(Trajectory_final)) %>%
    tidyr::pivot_longer(
      cols = tidyr::all_of(tp_levels),
      names_to = "Timepoint",
      values_to = "Z_Score"
    ) %>%
    dplyr::mutate(
      Timepoint = factor(Timepoint, levels = tp_levels),
      Trajectory = factor(Trajectory_final)
    )
  
  n_traj <- length(unique(data_long_traj$Trajectory))
  n_col  <- min(4, n_traj)
  
  traj_plot_final <- ggplot2::ggplot(
    data_long_traj,
    ggplot2::aes(x = Timepoint, y = Z_Score)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = feature_id),
      color = "grey",
      alpha = 0.2
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = Trajectory, color = Trajectory),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    ggplot2::facet_wrap(~ Trajectory, ncol = n_col) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0("Final WATER Trajectories (", n_traj, " groups)"),
      y = "Z-score"
    ) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  print(traj_plot_final)
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (water$save_outputs) {
    
    plot_dir  <- water_get_output_path(water, "correlate_secondary", "plots")
    table_dir <- water_get_output_path(water, "correlate_secondary", "tables")
    
    message("\n  Saving outputs:")
    message("    Tables: ", table_dir)
    message("    Plots:  ", plot_dir)
    
    water_save_output(
      data_z,
      "water_correlate_secondary_data",
      water,
      step = "correlate_secondary"
    )
    
    water_save_output(
      traj_summary,
      "water_final_trajectory_summary",
      water,
      step = "correlate_secondary"
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, "water_final_trajectories.png"),
      traj_plot_final, width = 12, height = 8, dpi = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, "water_final_trajectories.pdf"),
      traj_plot_final, width = 12, height = 8
    )
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$correlate_secondary <- list(
    data_z = data_z,
    final_assignments = data.frame(
      feature_id = data_z$feature_id,
      trajectory = data_z$Trajectory_final,
      stability_class = data_z$decision,
      stringsAsFactors = FALSE
    ),
    trajectory_summary = traj_summary,
    traj_plot_final = traj_plot_final
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step = "correlate_secondary",
    n_features = nrow(data_z),
    n_final_trajectories = n_traj,
    timestamp = Sys.time()
  )
  
  return(water)
}