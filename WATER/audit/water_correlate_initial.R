#' Initial correlation analysis for WATER classic pipeline
#'
#' Computes mean trajectory templates and correlates each feature against
#' all templates. Classifies features as Keep_original, Assign_new,
#' Weak, or Ambiguous based on correlation thresholds.
#'
#' @param water WATER object
#' @param r_threshold Numeric; correlation threshold for classification.
#'   Default 0.8.
#'
#' @return Updated WATER object with correlation results
#' @export

water_correlate_initial <- function(water, r_threshold = 0.8) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$label_classic)) {
    stop("Run water_label_classic() before water_correlate_initial()")
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
  
  message("Running initial correlation analysis (r_threshold = ", r_threshold, ")...")
  
  # ----------------------------
  # Retrieve data
  # ----------------------------
  data_z    <- water$label_classic$data_z
  tp_levels <- water$timepoint_levels
  
  # Ensure Trajectory is a factor
  data_z$Trajectory <- factor(data_z$Trajectory)
  
  # ----------------------------
  # Build trajectory templates
  # (mean z-score per trajectory per timepoint)
  # ----------------------------
  tmpl <- sapply(levels(data_z$Trajectory), function(lab) {
    idx <- which(data_z$Trajectory == lab)
    colMeans(
      as.matrix(data_z[idx, tp_levels, drop = FALSE]),
      na.rm = TRUE
    )
  })
  tmpl <- t(tmpl)  # rows = trajectories, cols = timepoints
  
  # ----------------------------
  # Correlate each feature against all templates
  # ----------------------------
  expr_mat <- as.matrix(data_z[, tp_levels, drop = FALSE])
  rownames(expr_mat) <- data_z$feature_id
  
  cor_all <- t(apply(expr_mat, 1, function(v) {
    stats::cor(t(tmpl), v, use = "pairwise.complete.obs")
  }))
  rownames(cor_all) <- data_z$feature_id
  colnames(cor_all) <- rownames(tmpl)
  
  # ----------------------------
  # Extract per-feature correlation metrics
  # ----------------------------
  own_idx  <- match(as.character(data_z$Trajectory), rownames(tmpl))
  r_own    <- cor_all[cbind(seq_len(nrow(cor_all)), own_idx)]
  r_best   <- apply(cor_all, 1, max, na.rm = TRUE)
  lab_best <- rownames(tmpl)[apply(cor_all, 1, which.max)]
  
  # ----------------------------
  # Attach to data_z
  # ----------------------------
  data_z$r_own     <- r_own
  data_z$r_best    <- r_best
  data_z$lab_best  <- lab_best
  
  # ----------------------------
  # Decision tree
  # ----------------------------
  data_z <- data_z %>%
    dplyr::mutate(
      decision = dplyr::case_when(
        # Keep: strong correlation to own trajectory and no better match
        r_own >= r_threshold & r_best == r_own ~
          "Keep_original",
        
        # Ambiguous: strong correlation to own AND another trajectory
        r_own >= r_threshold & r_best >= r_threshold & r_best > r_own ~
          "Ambiguous",
        
        # Weak: no strong correlation to any trajectory
        r_own < r_threshold & r_best < r_threshold ~
          "Weak",
        
        # Assign new: weak to own but strong to another
        r_own < r_threshold & r_best >= r_threshold ~
          "Assign_new"
      )
    )
  
  # ----------------------------
  # Apply Assign_new reassignments
  # ----------------------------
  data_z <- data_z %>%
    dplyr::mutate(
      Updated_trajectory = dplyr::case_when(
        decision == "Keep_original" ~ as.character(Trajectory),
        decision == "Assign_new"    ~ lab_best,
        decision == "Weak"          ~ "Weak",
        decision == "Ambiguous"     ~ "Ambiguous"
      )
    )
  
  # ----------------------------
  # Summary statistics
  # ----------------------------
  decision_counts <- data_z %>%
    dplyr::count(decision) %>%
    dplyr::rename(n_features = n) %>%
    dplyr::mutate(
      percent = round(100 * n_features / sum(n_features), 1)
    )
  
  collapsed_counts <- data_z %>%
    dplyr::mutate(
      collapsed = dplyr::case_when(
        decision == "Keep_original" ~ "Keep",
        decision == "Assign_new"    ~ "Reassign",
        decision %in% c("Ambiguous", "Weak") ~ "Re-cluster"
      )
    ) %>%
    dplyr::count(collapsed) %>%
    dplyr::rename(n_features = n) %>%
    dplyr::mutate(
      percent = round(100 * n_features / sum(n_features), 1)
    )
  
  # Per-trajectory breakdown
  traj_decision_counts <- data_z %>%
    dplyr::mutate(
      decision_collapsed = dplyr::case_when(
        decision == "Keep_original" ~ "Keep",
        decision == "Assign_new"    ~ "Reassign",
        decision %in% c("Ambiguous", "Weak") ~ "Re-cluster"
      )
    ) %>%
    dplyr::count(Trajectory, decision_collapsed) %>%
    dplyr::rename(n_features = n) %>%
    tidyr::complete(
      Trajectory,
      decision_collapsed,
      fill = list(n_features = 0)
    )
  
  # Confusion matrix
  trajectory_counts <- data_z %>%
    dplyr::group_by(Updated_trajectory, Trajectory) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop")
  
  confusion_matrix <- trajectory_counts %>%
    tidyr::pivot_wider(
      names_from  = Trajectory,
      values_from = n_features,
      values_fill = 0
    )
  
  # ----------------------------
  # Extract weak and ambiguous feature lists
  # ----------------------------
  weak_features <- data_z %>%
    dplyr::filter(decision == "Weak") %>%
    dplyr::pull(feature_id)
  
  ambiguous_features <- data_z %>%
    dplyr::filter(decision == "Ambiguous") %>%
    dplyr::pull(feature_id)
  
  # ----------------------------
  # Save outputs
  # ----------------------------
  if (water$save_outputs) {
    
    if (!exists("water_get_output_path") || !is.function(water_get_output_path)) {
      stop("water_get_output_path() not found.")
    }
    
    # --------------------------
    # Resolve path
    # --------------------------
    table_dir <- water_get_output_path(water, "correlate_initial", "tables")
    
    # --------------------------
    # Announce output directory
    # --------------------------
    message("\nSaving outputs to:")
    message("  Tables: ", table_dir)
    message("")
    
    # --------------------------
    # Save structured tables
    # --------------------------
    water_save_output(
      data_z,
      "water_correlate_initial_data",
      water,
      step = "correlate_initial"
    )
    
    water_save_output(
      decision_counts,
      "water_decision_counts",
      water,
      step = "correlate_initial"
    )
    
    water_save_output(
      collapsed_counts,
      "water_collapsed_counts",
      water,
      step = "correlate_initial"
    )
    
    water_save_output(
      traj_decision_counts,
      "water_per_trajectory_decisions",
      water,
      step = "correlate_initial"
    )
    
    water_save_output(
      confusion_matrix,
      "water_confusion_matrix",
      water,
      step = "correlate_initial"
    )
    
    # --------------------------
    # Save feature lists (raw txt)
    # --------------------------
    weak_path <- file.path(table_dir, "water_weak_features.txt")
    ambig_path <- file.path(table_dir, "water_ambiguous_features.txt")
    
    write.table(
      weak_features,
      weak_path,
      quote     = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    write.table(
      ambiguous_features,
      ambig_path,
      quote     = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    # Optional: explicit confirmation (nice touch)
    message("Saved feature lists:")
    message("  Weak features     → ", weak_path)
    message("  Ambiguous features → ", ambig_path)
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$correlate_initial <- list(
    data_z             = data_z,
    templates          = tmpl,
    cor_all            = cor_all,
    decision_counts    = decision_counts,
    collapsed_counts   = collapsed_counts,
    traj_decision_counts = traj_decision_counts,
    confusion_matrix   = confusion_matrix,
    weak_features      = weak_features,
    ambiguous_features = ambiguous_features,
    r_threshold        = r_threshold
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step               = "correlate_initial",
    r_threshold        = r_threshold,
    n_keep             = sum(data_z$decision == "Keep_original", na.rm = TRUE),
    n_assign_new       = sum(data_z$decision == "Assign_new",    na.rm = TRUE),
    n_weak             = length(weak_features),
    n_ambiguous        = length(ambiguous_features),
    timestamp          = Sys.time()
  )
  
  return(water)
}