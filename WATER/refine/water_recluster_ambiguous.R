#' Recluster ambiguous features from initial correlation analysis
#'
#' Subsets features classified as ambiguous, runs silhouette-guided
#' hierarchical clustering, and assigns trajectory labels. Follows the
#' same interactive workflow as water_classic().
#'
#' @param water WATER object
#' @param k_max Integer; maximum k for silhouette analysis. Defaults to
#'   n_timepoints * 7 + 2.
#' @param method Character; trajectory labeling method. One of "centroid"
#'   (default), "manual", or "none".
#' @param r_threshold Numeric; correlation threshold for centroid merging.
#'   Default 0.8.
#'
#' @return Updated WATER object with ambiguous reclustering results
#' @export

water_recluster_ambiguous <- function(
    water,
    k_max       = NULL,
    method      = "centroid",
    r_threshold = 0.8
) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$correlate_initial)) {
    stop("Run water_correlate_initial() before water_recluster_ambiguous()")
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required")
  }
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' required")
  }
  
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("Package 'factoextra' required")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required")
  }
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' required")
  }
  
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' required")
  }
  
  `%>%` <- magrittr::`%>%`
  
  # ----------------------------
  # Subset ambiguous features
  # ----------------------------
  ambiguous_features <- water$correlate_initial$ambiguous_features
  tp_levels          <- water$timepoint_levels
  
  if (length(ambiguous_features) == 0) {
    message("No ambiguous features found. Skipping reclustering.")
    return(water)
  }
  
  message("Reclustering ", length(ambiguous_features), " ambiguous features...")
  
  data_amb <- water$data_avg_norm[
    intersect(ambiguous_features, rownames(water$data_avg_norm)),
    tp_levels,
    drop = FALSE
  ]
  
  if (nrow(data_amb) == 0) {
    stop("No ambiguous features found in data_avg_norm")
  }
  
  # ----------------------------
  # Clean matrix
  # ----------------------------
  keep <- apply(data_amb, 1, function(x) all(is.finite(x))) &
    apply(data_amb, 1, sd) > 0
  
  data_amb <- data_amb[keep, , drop = FALSE]
  
  if (nrow(data_amb) < 2) {
    stop("Not enough ambiguous features after filtering")
  }
  
  # ----------------------------
  # Scale + distance
  # ----------------------------
  scaled_data <- scale(data_amb)
  diss        <- dist(scaled_data)
  
  n_timepoints <- ncol(data_amb)
  
  if (is.null(k_max)) {
    k_max <- n_timepoints * 7 + 2
  }
  
  # ----------------------------
  # Silhouette analysis
  # ----------------------------
  message("Running silhouette analysis across k = 2 to ", k_max, "...")
  message("This may take a few minutes for large datasets.")
  
  sil_plot <- factoextra::fviz_nbclust(
    scaled_data,
    FUN    = factoextra::hcut,
    method = "silhouette",
    k.max  = k_max,
    diss   = diss
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = "WATER Silhouette Analysis (Ambiguous features)",
      subtitle = paste0("k.max = ", k_max)
    )
  
  print(sil_plot)
  
  plot_dir <- water_get_output_path(water, "recluster_ambiguous", "plots")
  
  if (water$save_outputs) {
    png(file.path(plot_dir, "silhouette_ambiguous.png"),
        width = 1600, height = 1200, res = 150)
    print(sil_plot)
    dev.off()
  }
  
  # ----------------------------
  # Extract silhouette values
  # ----------------------------
  sil_df <- sil_plot$data
  
  if (!all(c("clusters", "y") %in% colnames(sil_df))) {
    stop("Unexpected silhouette plot structure from factoextra")
  }
  
  sil_df <- data.frame(
    k              = as.integer(as.character(sil_df$clusters)),
    avg_silhouette = sil_df$y
  )
  
  sil_vals <- sil_df$avg_silhouette
  ks       <- sil_df$k
  
  # ----------------------------
  # Select candidate k values
  # ----------------------------
  k_good <- ks[sil_vals > 0.5]
  
  drops          <- diff(sil_vals)
  drop_positions <- which(drops <= -0.075)
  k_tier1        <- ks[drop_positions]
  k_tier1        <- k_tier1[k_tier1 %in% k_good]
  
  k_tier2 <- if (length(k_good) > 0) max(k_good) else integer(0)
  
  k_recommended <- sort(unique(c(k_tier1, k_tier2)))
  
  message("Generating heatmaps for inspection...")
  
  # ----------------------------
  # hclust + heatmaps
  # ----------------------------
  hc <- hclust(diss, method = "ward.D2")
  
  color_palette <- colorRampPalette(
    c("#408493", "#94D2BD", "#E9D8A6", "#EE9C21", "#B73A27")
  )
  
  save_pheatmap_png <- function(x, filename) {
    png(filename, width = 1400, height = 3000, res = 150)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  grob_list <- list()
  
  if (length(k_recommended) > 0) {
    
    for (k in k_recommended) {
      message("Generating heatmap for k = ", k)
      
      ph <- pheatmap::pheatmap(
        data_amb,
        scale         = "row",
        cluster_rows  = TRUE,
        cluster_cols  = FALSE,
        color         = color_palette(100),
        cutree_rows   = k,
        show_rownames = FALSE,
        main          = paste0("k = ", k, " (sil = ",
                               round(sil_df$avg_silhouette[sil_df$k == k], 2), ")"),
        silent        = TRUE
      )
      
      grob_list[[as.character(k)]] <- ph$gtable
      
      if (water$save_outputs) {
        save_pheatmap_png(
          ph,
          file.path(plot_dir, paste0("heatmap_ambiguous_k_", k, ".png"))
        )
      }
    }
    
    n <- length(grob_list)
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow    = 2,
          ncol    = n,
          heights = grid::unit(c(0.05, 0.95), "npc")
        )
      )
    )
    
    grid::pushViewport(
      grid::viewport(layout.pos.row = 1, layout.pos.col = 1:n)
    )
    grid::grid.text(
      "Recommended k values for ambiguous reclustering",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    grid::upViewport()
    
    for (i in seq_along(grob_list)) {
      grid::pushViewport(
        grid::viewport(layout.pos.row = 2, layout.pos.col = i)
      )
      grid::grid.draw(grob_list[[i]])
      grid::upViewport()
    }
    grid::upViewport()
  }
  
  # ----------------------------
  # Prompt user for k_final
  # ----------------------------
  k_default <- if (length(k_recommended) > 0) max(k_recommended) else 2
  
  message("\nSelect number of clusters (k) for ambiguous reclustering.")
  message("Recommended k values: ", paste(k_recommended, collapse = ", "))
  message("Press Enter to use the default (k = ", k_default, ")")
  message("Or type a k value and press Enter: ")
  
  input   <- trimws(readline())
  k_final <- if (nchar(input) == 0) {
    message("Using default k = ", k_default)
    k_default
  } else {
    k_in <- suppressWarnings(as.integer(input))
    if (is.na(k_in)) {
      warning("Invalid input. Defaulting to k = ", k_default)
      k_default
    } else {
      k_in
    }
  }
  
  # ----------------------------
  # Cut tree from pheatmap
  # ----------------------------
  message("Cutting ambiguous tree at k = ", k_final)
  
  ph_final <- pheatmap::pheatmap(
    data_amb,
    scale         = "row",
    cluster_rows  = TRUE,
    cluster_cols  = FALSE,
    color         = color_palette(100),
    cutree_rows   = k_final,
    show_rownames = FALSE,
    main          = paste0("k = ", k_final),
    silent        = TRUE
  )
  
  clusters <- cutree(ph_final$tree_row, k = k_final)
  
  cluster_df <- data.frame(
    feature_id = names(clusters),
    Cluster    = factor(clusters),
    stringsAsFactors = FALSE
  )
  
  cluster_size_df <- cluster_df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop")
  
  message("Ambiguous cluster sizes:")
  print(as.data.frame(cluster_size_df))
  
  # ----------------------------
  # Row-scale for plotting
  # ----------------------------
  data_z            <- t(scale(t(data_amb)))
  data_z            <- as.data.frame(data_z)
  data_z$Cluster    <- factor(clusters[rownames(data_z)])
  data_z$feature_id <- rownames(data_z)
  
  # ----------------------------
  # Faceted trajectory plot
  # ----------------------------
  n_col <- min(5, k_final)
  
  data_long <- data_z %>%
    tidyr::pivot_longer(
      cols      = tidyr::all_of(tp_levels),
      names_to  = "Timepoint",
      values_to = "Z_Score"
    ) %>%
    dplyr::mutate(
      Timepoint = factor(Timepoint, levels = tp_levels)
    )
  
  traj_plot <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(x = Timepoint, y = Z_Score)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = feature_id),
      color = "grey", alpha = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = Cluster),
      fun = mean, geom = "line", linewidth = 1
    ) +
    ggplot2::facet_wrap(~ Cluster, ncol = n_col) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0("Ambiguous reclustering (k = ", k_final, ", ",
                     nrow(data_amb), " features)"),
      x = "", y = "Row Z-score"
    ) +
    ggplot2::theme(
      legend.position    = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )
  
  print(traj_plot)
  
  if (water$save_outputs) {
    ggplot2::ggsave(
      file.path(plot_dir, paste0("trajectory_ambiguous_k", k_final, ".png")),
      traj_plot, width = 12, height = 6, dpi = 300
    )
    ggplot2::ggsave(
      file.path(plot_dir, paste0("trajectory_ambiguous_k", k_final, ".pdf")),
      traj_plot, width = 12, height = 6
    )
  }
  
  # ----------------------------
  # Trajectory labeling method
  # ----------------------------
  method_specified <- !missing(method)
  
  if (!method_specified) {
    
    message("\nLabeling Ambiguous Trajectories.")
    message("How would you like to label ambiguous trajectories?")
    message("  1. Centroid clustering (default, r_threshold = ", r_threshold, ")")
    message("  2. Manual label_map")
    message("  3. None - keep all ", k_final, " clusters as separate trajectories")
    message("Press Enter for default (centroid), or type 1, 2, or 3: ")
    
    method_input <- trimws(readline())
    
    if (nchar(method_input) == 0 || method_input == "1") {
      method <- "centroid"
    } else if (method_input == "2") {
      method <- "manual"
    } else if (method_input == "3") {
      message("Warning: this will keep all ", k_final,
              " clusters as individual trajectories.")
      message("Are you sure? Type 'yes' to confirm, or press Enter for centroid: ")
      confirm <- trimws(readline())
      if (tolower(confirm) == "yes") {
        method <- "none"
      } else {
        message("Defaulting to centroid clustering.")
        method <- "centroid"
      }
    } else {
      message("Invalid input. Defaulting to centroid clustering.")
      method <- "centroid"
    }
  }
  
  # ----------------------------
  # Build centroids
  # ----------------------------
  centroids <- do.call(rbind, lapply(levels(cluster_df$Cluster), function(cl) {
    ids <- cluster_df$feature_id[cluster_df$Cluster == cl]
    colMeans(data_z[ids, tp_levels, drop = FALSE], na.rm = TRUE)
  }))
  rownames(centroids) <- levels(cluster_df$Cluster)
  
  # ----------------------------
  # Method: centroid
  # ----------------------------
  if (method == "centroid") {
    
    message("Computing centroid correlations (r_threshold = ", r_threshold, ")...")
    R        <- cor(t(centroids), method = "pearson")
    cl_names <- rownames(centroids)
    groups   <- list()
    
    for (i in seq_along(cl_names)) {
      current <- cl_names[i]
      placed  <- FALSE
      if (length(groups) > 0) {
        for (g in seq_along(groups)) {
          members <- groups[[g]]
          if (all(R[current, members] >= r_threshold)) {
            groups[[g]] <- c(members, current)
            placed <- TRUE
            break
          }
        }
      }
      if (!placed) groups[[length(groups) + 1]] <- current
    }
    
    trajectory_map <- setNames(rep(NA_character_, length(cl_names)), cl_names)
    for (i in seq_along(groups)) {
      trajectory_map[groups[[i]]] <- paste0("Ambiguous_Trajectory_", i)
    }
    
    message("Ambiguous clusters merged into ", length(groups), " trajectories:")
    print(trajectory_map)
    
    # ----------------------------
    # Method: manual
    # ----------------------------
  } else if (method == "manual") {
    
    message("\nManual Labeling for Ambiguous Features.")
    message("Cluster sizes:")
    print(as.data.frame(cluster_size_df), row.names = FALSE)
    
    template_df <- data.frame(
      Trajectory = "",
      Cluster    = as.integer(levels(cluster_df$Cluster)),
      n_features = cluster_size_df$n_features,
      stringsAsFactors = FALSE
    )
    
    template_path <- file.path(
      water_get_output_path(water, "recluster_ambiguous", "tables"),
      "ambiguous_label_map_template.csv"
    )
    
    write.csv(template_df, template_path, row.names = FALSE)
    message("\nTemplate CSV written to: ", template_path)
    message("Fill in the Trajectory column, save, then provide the path below.")
    message("Accepted formats: .csv, .txt, or .R file")
    message("\nFile path: ")
    
    path_input <- trimws(readline())
    
    if (nchar(path_input) == 0 || !file.exists(path_input)) {
      stop(
        "File not found or no path provided.\n",
        "Re-run water_recluster_ambiguous(water, method = 'manual')"
      )
    }
    
    label_map <- tryCatch({
      ext <- tolower(tools::file_ext(path_input))
      if (ext == "csv") {
        df <- read.csv(path_input, stringsAsFactors = FALSE)
        if (!all(c("Trajectory", "Cluster") %in% colnames(df))) {
          stop("CSV must have columns named 'Trajectory' and 'Cluster'")
        }
        split(as.integer(df$Cluster), df$Trajectory)
      } else if (ext %in% c("txt", "r")) {
        env <- new.env()
        source(path_input, local = env)
        if (!exists("label_map", envir = env)) {
          stop("File must contain an object named 'label_map'")
        }
        get("label_map", envir = env)
      } else {
        stop("Unsupported file type: .", ext)
      }
    }, error = function(e) {
      stop("Could not parse label_map file: ", e$message)
    })
    
    trajectory_map <- setNames(
      rep(names(label_map), lengths(label_map)),
      as.character(unlist(label_map))
    )
    
    missing <- setdiff(as.character(levels(cluster_df$Cluster)), names(trajectory_map))
    if (length(missing) > 0) {
      warning("Unassigned clusters labeled 'Unassigned': ",
              paste(missing, collapse = ", "))
      trajectory_map[missing] <- "Unassigned"
    }
    
    # ----------------------------
    # Method: none
    # ----------------------------
  } else if (method == "none") {
    trajectory_map <- setNames(
      paste0("Ambiguous_Cluster_", levels(cluster_df$Cluster)),
      as.character(levels(cluster_df$Cluster))
    )
  }
  
  # ----------------------------
  # Rename trajectories (centroid and none)
  # ----------------------------
  if (method %in% c("centroid", "none")) {
    
    message("\nRenaming Ambiguous Trajectories.")
    message("Auto-labeling produced ", length(unique(trajectory_map)),
            " trajectories with generic names.")
    message("Inspect the trajectory plots, then provide a rename file.")
    
    unique_trajs <- sort(unique(trajectory_map))
    
    rename_df <- data.frame(
      Generic_name     = unique_trajs,
      Biological_label = "",
      n_features       = sapply(unique_trajs, function(t) sum(trajectory_map == t)),
      stringsAsFactors = FALSE
    )
    
    rename_path <- file.path(
      water_get_output_path(water, "recluster_ambiguous", "tables"),
      "ambiguous_rename_template.csv"
    )
    
    write.csv(rename_df, rename_path, row.names = FALSE)
    message("\nRename template written to: ", rename_path)
    message("Fill in the Biological_label column, save, then")
    message("provide the path below and press Enter.")
    message("Accepted formats: .csv, .txt, or .R file")
    message("\nFile path: ")
    
    rename_input <- trimws(readline())
    
    if (nchar(rename_input) == 0 || !file.exists(rename_input)) {
      message("No rename file provided. Keeping generic trajectory names.")
    } else {
      
      rename_map <- tryCatch({
        ext <- tolower(tools::file_ext(rename_input))
        if (ext == "csv") {
          df <- read.csv(rename_input, stringsAsFactors = FALSE)
          if (!all(c("Generic_name", "Biological_label") %in% colnames(df))) {
            stop("CSV must have columns 'Generic_name' and 'Biological_label'")
          }
          setNames(df$Biological_label, df$Generic_name)
        } else if (ext %in% c("txt", "r")) {
          env <- new.env()
          source(rename_input, local = env)
          if (!exists("rename_map", envir = env)) {
            stop("File must contain an object named 'rename_map'")
          }
          get("rename_map", envir = env)
        } else {
          stop("Unsupported file type: .", ext)
        }
      }, error = function(e) {
        warning("Could not parse rename file: ", e$message,
                "\nKeeping generic trajectory names.")
        NULL
      })
      
      if (!is.null(rename_map)) {
        renamed <- rename_map[trajectory_map]
        valid   <- !is.na(renamed) & nchar(trimws(renamed)) > 0
        trajectory_map[valid] <- renamed[valid]
        message("Ambiguous trajectories renamed successfully.")
      }
    }
  }
  
  # ----------------------------
  # Apply trajectory labels
  # ----------------------------
  cluster_df$Trajectory <- trajectory_map[as.character(cluster_df$Cluster)]
  cluster_df$Trajectory[is.na(cluster_df$Trajectory)] <- "Unassigned"
  
  data_z$Trajectory <- cluster_df$Trajectory[
    match(data_z$feature_id, cluster_df$feature_id)
  ]
  
  traj_size_df <- cluster_df %>%
    dplyr::group_by(Trajectory) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_features))
  
  message("Ambiguous trajectory sizes:")
  print(as.data.frame(traj_size_df), row.names = FALSE)
  
  # ----------------------------
  # Simplified trajectory plot
  # ----------------------------
  data_long_traj <- data_z %>%
    tidyr::pivot_longer(
      cols      = tidyr::all_of(tp_levels),
      names_to  = "Timepoint",
      values_to = "Z_Score"
    ) %>%
    dplyr::mutate(
      Timepoint  = factor(Timepoint, levels = tp_levels),
      Trajectory = factor(Trajectory)
    )
  
  n_traj <- length(unique(data_long_traj$Trajectory))
  n_col  <- min(3, n_traj)
  
  traj_plot_labeled <- ggplot2::ggplot(
    data_long_traj,
    ggplot2::aes(x = Timepoint, y = Z_Score)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = feature_id),
      color = "grey", alpha = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = Trajectory, color = Trajectory),
      fun = mean, geom = "line", linewidth = 1
    ) +
    ggplot2::facet_wrap(~ Trajectory, ncol = n_col) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0("Ambiguous trajectories (", n_traj, " groups)"),
      x = "", y = "Z-score"
    ) +
    ggplot2::theme(
      legend.position    = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )
  
  print(traj_plot_labeled)
  
  if (water$save_outputs) {
    ggplot2::ggsave(
      file.path(plot_dir, "trajectory_ambiguous_labeled.png"),
      traj_plot_labeled, width = 10, height = 6, dpi = 300
    )
    ggplot2::ggsave(
      file.path(plot_dir, "trajectory_ambiguous_labeled.pdf"),
      traj_plot_labeled, width = 10, height = 6
    )
    
    water_save_output(
      cluster_df, "water_ambiguous_assignments",
      water, step = "recluster_ambiguous"
    )
    
    water_save_output(
      traj_size_df, "water_ambiguous_trajectory_sizes",
      water, step = "recluster_ambiguous"
    )
    
    water_save_output(
      data.frame(
        Cluster    = names(trajectory_map),
        Trajectory = unname(trajectory_map),
        stringsAsFactors = FALSE
      ),
      "water_ambiguous_trajectory_map",
      water, step = "recluster_ambiguous"
    )
  }
  
  # ----------------------------
  # Integrate with initial trajectories
  # ----------------------------
  corr_data  <- water$correlate_initial$data_z
  initial_df <- water$label_classic$cluster_df
  
  integrated_df <- data.frame(
    feature_id = initial_df$feature_id,
    stringsAsFactors = FALSE
  )
  
  integrated_df$decision <- corr_data$decision[
    match(integrated_df$feature_id, corr_data$feature_id)
  ]
  
  integrated_df$Trajectory <- dplyr::case_when(
    integrated_df$decision == "Keep_original" ~
      as.character(corr_data$Trajectory[
        match(integrated_df$feature_id, corr_data$feature_id)]),
    integrated_df$decision == "Assign_new" ~
      corr_data$lab_best[
        match(integrated_df$feature_id, corr_data$feature_id)],
    integrated_df$decision == "Ambiguous" ~
      cluster_df$Trajectory[
        match(integrated_df$feature_id, cluster_df$feature_id)],
    integrated_df$decision == "Weak" ~
      "Pending_weak_reclustering",
    TRUE ~ "Unassigned"
  )
  
  n_pending    <- sum(integrated_df$Trajectory == "Pending_weak_reclustering",
                      na.rm = TRUE)
  n_unassigned <- sum(integrated_df$Trajectory == "Unassigned", na.rm = TRUE)
  
  message("\nIntegration Summary")
  message("Keep_original: ",
          sum(integrated_df$decision == "Keep_original", na.rm = TRUE))
  message("Assign_new:    ",
          sum(integrated_df$decision == "Assign_new",    na.rm = TRUE))
  message("Ambiguous:     ",
          sum(integrated_df$decision == "Ambiguous",     na.rm = TRUE))
  message("Weak (pending): ", n_pending)
  if (n_unassigned > 0) {
    warning(n_unassigned, " features could not be assigned.")
  }
  
  if (water$save_outputs) {
    water_save_output(
      integrated_df,
      "water_integrated_trajectories_pre_weak",
      water,
      step = "recluster_ambiguous"
    )
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$recluster_ambiguous <- list(
    k_final           = k_final,
    hclust            = ph_final$tree_row,
    cluster_df        = cluster_df,
    cluster_sizes     = cluster_size_df,
    data_z            = data_z,
    centroids         = centroids,
    trajectory_map    = trajectory_map,
    trajectory_sizes  = traj_size_df,
    silhouette        = sil_df,
    k_recommended     = k_recommended,
    method            = method,
    r_threshold       = if (method == "centroid") r_threshold else NULL,
    traj_plot         = traj_plot,
    traj_plot_labeled = traj_plot_labeled,
    integrated_df     = integrated_df
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step           = "recluster_ambiguous",
    k_final        = k_final,
    method         = method,
    n_features     = nrow(data_amb),
    n_trajectories = length(unique(cluster_df$Trajectory)),
    n_pending_weak = n_pending,
    timestamp      = Sys.time()
  )
  
  return(water)
}