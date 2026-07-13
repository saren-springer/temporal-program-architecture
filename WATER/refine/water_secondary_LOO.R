#' Secondary LOO stability analysis for WATER classic pipeline
#'
#' Performs leave-one-replicate-out stability analysis using FINAL
#' integrated trajectory assignments after weak/ambiguous reclustering.
#'
#' @param water WATER object
#' @param seed Integer; random seed. Default 7
#' @param r_threshold Numeric; correlation threshold (default 0.8)
#'
#' @return Updated WATER object with secondary LOO results
#' @export

water_secondary_LOO <- function(water, seed = 7, r_threshold = 0.8) {
  
  set.seed(seed)
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$recluster_weak)) {
    stop("Run water_recluster_weak() before water_secondary_LOO()")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required")
  if (!requireNamespace("magrittr", quietly = TRUE)) stop("Package 'magrittr' required")
  
  `%>%` <- magrittr::`%>%`
  
  # ----------------------------
  # Retrieve data
  # ----------------------------
  tp_levels <- water$timepoint_levels
  tp_labels <- water$timepoint_labels
  data_rep  <- water$data_norm
  
  integrated_df <- water$recluster_weak$integrated_df
  
  baseline_traj <- setNames(
    integrated_df$Trajectory,
    integrated_df$feature_id
  )
  
  features <- intersect(names(baseline_traj), rownames(data_rep))
  data_rep <- data_rep[features, , drop = FALSE]
  
  # ----------------------------
  # Helper: Z-score
  # ----------------------------
  rowz <- function(M) t(scale(t(M)))
  
  # ----------------------------
  # Build baseline matrix
  # ----------------------------
  M_base <- do.call(cbind, lapply(tp_levels, function(tp) {
    cols <- which(tp_labels == tp)
    rowMeans(data_rep[, cols, drop = FALSE])
  }))
  
  colnames(M_base) <- tp_levels
  rownames(M_base) <- features
  
  Z_base <- rowz(M_base)
  
  # ----------------------------
  # Build trajectory templates
  # ----------------------------
  build_templates <- function(Z, traj_labels) {
    labs <- unique(traj_labels)
    
    do.call(rbind, lapply(labs, function(lab) {
      idx <- names(traj_labels)[traj_labels == lab]
      colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
    })) |> `rownames<-`(labs)
  }
  
  templates_base <- build_templates(Z_base, baseline_traj)
  
  # ----------------------------
  # LOO loop
  # ----------------------------
  rep_cols <- colnames(data_rep)
  loo_summaries <- list()
  
  message("Running LOO across ", length(rep_cols), " replicates...")
  
  transition_list <- list()
  gene_shift_list <- list()
  
  for (drop_col in rep_cols) {
    
    message("  Dropping: ", drop_col)
    
    # ----------------------------
    # Recompute matrix
    # ----------------------------
    M_loo <- do.call(cbind, lapply(tp_levels, function(tp) {
      cols <- which(tp_labels == tp)
      cols <- cols[colnames(data_rep)[cols] != drop_col]
      
      if (length(cols) == 0) return(NULL)
      
      rowMeans(data_rep[, cols, drop = FALSE])
    }))
    
    if (is.null(M_loo) || ncol(M_loo) < length(tp_levels)) {
      message("  Skipping (missing timepoint)")
      next
    }
    
    colnames(M_loo) <- tp_levels
    rownames(M_loo) <- features
    
    Z_loo <- rowz(M_loo)
    
    # ----------------------------
    # Correlate genes to templates
    # ----------------------------
    cor_all <- t(apply(Z_loo, 1, function(v) {
      stats::cor(t(templates_base), v, use = "pairwise.complete.obs")
    }))
    
    rownames(cor_all) <- rownames(Z_loo)
    colnames(cor_all) <- rownames(templates_base)
    
    # ----------------------------
    # Assign best trajectory (NO LOGIC)
    # ----------------------------
    lab_best <- colnames(cor_all)[apply(cor_all, 1, which.max)]
    
    baseline_vec <- baseline_traj[rownames(Z_loo)]
    
    df_trans <- data.frame(
      feature_id = rownames(Z_loo),
      from = baseline_vec,
      to   = lab_best,
      stringsAsFactors = FALSE
    )
    
    # keep only true changes
    df_trans <- df_trans[df_trans$from != df_trans$to, ]
    
    transition_list[[drop_col]] <- df_trans
    
    df_gene <- data.frame(
      feature_id = rownames(Z_loo),
      baseline   = baseline_vec,
      reassigned = lab_best,
      shifted    = lab_best != baseline_vec,
      replicate  = drop_col,
      stringsAsFactors = FALSE
    )
    
    gene_shift_list[[drop_col]] <- df_gene
    
    # ----------------------------
    # Confidence (separate!)
    # ----------------------------
    r_best <- apply(cor_all, 1, max, na.rm = TRUE)
    low_conf <- r_best < r_threshold
    
    # ----------------------------
    # TRUE retention metric
    # ----------------------------
    same <- lab_best == baseline_vec
    
    pct_retention <- mean(same, na.rm = TRUE)
    pct_reassign  <- mean(!same, na.rm = TRUE)
    pct_lowconf   <- mean(low_conf, na.rm = TRUE)
    
    loo_summaries[[drop_col]] <- data.frame(
      ReplicateDropped = drop_col,
      Retention        = pct_retention,
      Reassign         = pct_reassign,
      Low_confidence   = pct_lowconf,
      stringsAsFactors = FALSE
    )
  }
  
  # ----------------------------
  # Combine results
  # ----------------------------
  summary_df <- do.call(rbind, loo_summaries)
  
  all_transitions <- do.call(rbind, transition_list)
  
  transition_counts <- all_transitions %>%
    dplyr::count(from, to, name = "n") %>%
    dplyr::arrange(desc(n))
  
  transition_probs <- transition_counts %>%
    dplyr::group_by(from) %>%
    dplyr::mutate(prob = n / sum(n)) %>%
    dplyr::ungroup()
  
  gene_shift_df <- do.call(rbind, gene_shift_list)
  
  gene_instability <- gene_shift_df %>%
    dplyr::group_by(feature_id) %>%
    dplyr::summarise(
      n_shifts = sum(shifted, na.rm = TRUE),
      n_tests  = dplyr::n(),
      shift_fraction = n_shifts / n_tests,
      baseline = dplyr::first(baseline),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(shift_fraction))
  
  gene_instability <- gene_instability %>%
    dplyr::mutate(
      stability_class = dplyr::case_when(
        shift_fraction == 0 ~ "Stable",
        shift_fraction < 0.25 ~ "Mostly stable",
        shift_fraction < 0.5 ~ "Moderately unstable",
        TRUE ~ "Highly unstable"
      )
    )
  
  message("\nSecondary LOO Summary:")
  message(sprintf("Mean Retention      = %.1f%%", 100 * mean(summary_df$Retention)))
  message(sprintf("Mean Reassign       = %.1f%%", 100 * mean(summary_df$Reassign)))
  message(sprintf("Mean Low confidence = %.1f%%", 100 * mean(summary_df$Low_confidence)))
  
  # ----------------------------
  # Plot
  # ----------------------------
  summary_long <- summary_df %>%
    tidyr::pivot_longer(
      cols = c(Retention, Reassign, Low_confidence),
      names_to = "Metric",
      values_to = "Fraction"
    ) %>%
    dplyr::mutate(Percent = Fraction * 100)
  
  p <- ggplot2::ggplot(summary_long,
                       ggplot2::aes(x = ReplicateDropped, y = Percent, fill = Metric)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Secondary LOO Stability (Final Trajectories)",
      y = "Percent",
      x = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  print(p)
  
  p_trans <- ggplot2::ggplot(transition_counts,
                             ggplot2::aes(x = from, y = to, fill = n)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Trajectory Transition Frequency (LOO)",
      x = "Baseline Trajectory",
      y = "Reassigned Trajectory",
      fill = "Count"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  print(p_trans)
  
  p_gene <- ggplot2::ggplot(gene_instability,
                            ggplot2::aes(x = shift_fraction)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Gene-level trajectory instability",
      x = "Fraction of LOO runs where trajectory changes",
      y = "Number of genes"
    )
  
  print(p_gene)
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (!is.null(water$save_outputs) && water$save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    # --------------------------
    # Resolve paths
    # --------------------------
    plot_dir  <- water_get_output_path(water, "loo_secondary", "plots")
    table_dir <- water_get_output_path(water, "loo_secondary", "tables")
    
    # --------------------------
    # Announce output directories
    # --------------------------
    message("\nSaving outputs to:")
    message("  Tables: ", table_dir)
    message("  Plots:  ", plot_dir)
    message("")
    
    # --------------------------
    # Save summary tables
    # --------------------------
    water_save_output(
      summary_df,
      "water_secondary_LOO_summary",
      water,
      step = "loo_secondary"
    )
    
    water_save_output(
      summary_long,
      "water_secondary_LOO_long",
      water,
      step = "loo_secondary"
    )
    
    # --------------------------
    # Save transition analysis
    # --------------------------
    water_save_output(
      transition_counts,
      "water_transition_counts",
      water,
      step = "loo_secondary"
    )
    
    water_save_output(
      transition_probs,
      "water_transition_probs",
      water,
      step = "loo_secondary"
    )
    
    message("Trajectory transition analysis saved.")
    
    # --------------------------
    # Save gene-level instability
    # --------------------------
    water_save_output(
      gene_instability,
      "water_gene_instability",
      water,
      step = "loo_secondary"
    )
    
    water_save_output(
      gene_shift_df,
      "water_gene_shift_long",
      water,
      step = "loo_secondary"
    )
    
    message("Gene-level instability metrics saved.")
    
    # --------------------------
    # Save plots
    # --------------------------
    ggplot2::ggsave(
      file.path(plot_dir, "secondary_LOO_stability.png"),
      p,
      width = 10, height = 6, dpi = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, "trajectory_transition_heatmap.png"),
      p_trans,
      width = 8, height = 6, dpi = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, "gene_instability_histogram.png"),
      p_gene,
      width = 8, height = 6, dpi = 300
    )
    
    message("LOO diagnostic plots saved.")
  }
  
  # ----------------------------
  # Store
  # ----------------------------
  water$loo_secondary <- list(
    summary = summary_df,
    plot    = p
  )
  
  water$loo_secondary$transitions <- list(
    counts = transition_counts,
    probs  = transition_probs,
    plot   = p_trans
  )
  
  water$loo_secondary$gene_instability <- gene_instability
  water$loo_secondary$gene_shift_long <- gene_shift_df
  
  water$history[[length(water$history) + 1]] <- list(
    step = "loo_secondary",
    mean_retention = mean(summary_df$Retention),
    mean_reassign  = mean(summary_df$Reassign),
    mean_lowconf   = mean(summary_df$Low_confidence),
    timestamp = Sys.time()
  )
  
  return(water)
}