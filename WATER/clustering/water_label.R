#' Assign trajectory labels to WATER clusters
#'
#' Collapses clusters into named trajectories using one of three methods:
#' centroid-based correlation (default), manual label_map, or none.
#'
#' @param water WATER object
#' @param label_map Named list mapping trajectory names to cluster numbers
#'   (e.g. list("T1 High" = c(1, 3))). Only used when method = "manual".
#' @param method Character; one of "centroid" (default), "manual", or "none"
#' @param r_threshold Numeric; correlation threshold for centroid merging.
#'   Default 0.8. Only used when method = "centroid".
#'
#' @return Updated WATER object with trajectory assignments
#' @export

water_label <- function(
    water,
    label_map   = NULL,
    method      = "centroid",
    r_threshold = 0.8
) {
  
  # ----------------------------
  # Checks
  # ----------------------------
  if (is.null(water$cut_classic)) {
    stop("Run water_cut() before water_label_classic()")
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
  
  # Override method if label_map provided directly
  if (!is.null(label_map)) {
    method <- "manual"
  }
  
  if (!method %in% c("centroid", "manual", "none")) {
    stop("method must be one of 'centroid', 'manual', or 'none'")
  }
  
  message("Labeling trajectories using method: ", method)
  
  # ----------------------------
  # Retrieve stored objects
  # ----------------------------
  data_z     <- water$cut_classic$data_z
  cluster_df <- water$cut_classic$cluster_df
  k_final    <- water$cut_classic$k_final
  tp_levels  <- water$timepoint_levels
  
  # ----------------------------
  # Build cluster centroids
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
    
    message("Computing pairwise centroid correlations (r_threshold = ",
            r_threshold, ")...")
    
    R <- cor(t(centroids), method = "pearson")
    
    cl_names <- rownames(centroids)
    groups <- list()
    
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
      
      if (!placed) {
        groups[[length(groups) + 1]] <- current
      }
    }
    
    trajectory_map <- setNames(
      rep(NA_character_, length(cl_names)),
      cl_names
    )
    
    for (i in seq_along(groups)) {
      trajectory_map[groups[[i]]] <- paste0("Trajectory_", i)
    }
    
    message("Clusters merged into ", length(groups), " trajectories:")
    print(trajectory_map)
    
    # ----------------------------
    # Plot trajectories BEFORE rename prompt
    # ----------------------------
    data_z$Trajectory <- trajectory_map[as.character(cluster_df$Cluster)]
    data_z$Trajectory[is.na(data_z$Trajectory)] <- "Unassigned"
    
    data_long <- data_z %>%
      tidyr::pivot_longer(
        cols      = tidyr::all_of(tp_levels),
        names_to  = "Timepoint",
        values_to = "Z_Score"
      ) %>%
      dplyr::mutate(
        Timepoint  = factor(Timepoint, levels = tp_levels),
        Trajectory = factor(Trajectory)
      )
    
    n_traj <- length(unique(data_long$Trajectory))
    n_col  <- min(3, n_traj)
    
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
        ggplot2::aes(group = Trajectory, color = Trajectory),
        fun       = mean,
        geom      = "line",
        linewidth = 1
      ) +
      ggplot2::facet_wrap(~ Trajectory, ncol = n_col) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0("WATER Trajectories (", n_traj, " groups)"),
        x     = "",
        y     = "Z-score"
      ) +
      ggplot2::theme(
        legend.position    = "none",
        strip.text         = ggplot2::element_text(face = "bold")
      )
    
    print(traj_plot)
    
    # ----------------------------
    # Method: manual
    # ----------------------------
  } else if (method == "manual") {
    
    if (is.null(label_map)) {
      
      # Print cluster plot and sizes for reference
      print(water$cut_classic$trajectory_plot)
      message("\n--- Manual Trajectory Labeling ---")
      message("Cluster sizes:")
      print(as.data.frame(water$cut_classic$cluster_sizes), row.names = FALSE)
      
      # Write template CSV
      template_df <- data.frame(
        Trajectory = "",
        Cluster    = as.integer(levels(cluster_df$Cluster)),
        n_features = water$cut_classic$cluster_sizes$n_features,
        stringsAsFactors = FALSE
      )
      
      template_path <- file.path(
        water_get_output_path(water, "label_classic", "tables"),
        "label_map_template.csv"
      )
      
      write.csv(template_df, template_path, row.names = FALSE)
      message("\nTemplate CSV written to: ", template_path)
      message("Fill in the Trajectory column, save the file, then")
      message("provide the path below and press Enter.")
      message("Accepted formats: .csv, .txt, or .R file")
      message("\nFile path: ")
      
      path_input <- trimws(readline())
      
      if (nchar(path_input) == 0 || !file.exists(path_input)) {
        stop(
          "File not found or no path provided.\n",
          "Re-run water_label_classic(water, method = 'manual') and ",
          "provide a valid file path."
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
          stop(
            "Unsupported file type: .", ext, "\n",
            "Provide a .csv, .txt, or .R file."
          )
        }
        
      }, error = function(e) {
        stop("Could not parse label_map file: ", e$message)
      })
    }
    
    # Build cluster -> trajectory lookup
    trajectory_map <- setNames(
      rep(names(label_map), lengths(label_map)),
      as.character(unlist(label_map))
    )
    
    # Check all clusters are covered
    all_clusters <- as.character(levels(cluster_df$Cluster))
    missing      <- setdiff(all_clusters, names(trajectory_map))
    
    if (length(missing) > 0) {
      warning(
        "The following clusters are not in label_map and will be ",
        "labeled 'Unassigned': ",
        paste(missing, collapse = ", ")
      )
      trajectory_map[missing] <- "Unassigned"
    }
    
    # ----------------------------
    # Method: none
    # ----------------------------
  } else if (method == "none") {
    
    trajectory_map <- setNames(
      paste0("Cluster_", levels(cluster_df$Cluster)),
      as.character(levels(cluster_df$Cluster))
    )
  }
  
  # ----------------------------
  # Rename centroid trajectories
  # ----------------------------
  if (method %in% c("centroid", "none")) {
    
    message("\nRenaming Trajectories.")
    message("Auto-labeling produced ", length(unique(trajectory_map)),
            " trajectories with generic names.")
    message("Inspect the trajectory plot above, then provide a rename file.")
    message("The file should map generic names to biological labels.")
    
    # Write rename template
    unique_trajs <- sort(unique(trajectory_map))
    
    rename_df <- data.frame(
      Generic_name   = unique_trajs,
      Biological_label = "",
      n_features     = sapply(unique_trajs, function(t) {
        sum(trajectory_map == t)
      }),
      stringsAsFactors = FALSE
    )
    
    rename_path <- file.path(
      water_get_output_path(water, "label_classic", "tables"),
      "trajectory_rename_template.csv"
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
        # Apply rename - only replace non-empty labels
        renamed <- rename_map[trajectory_map]
        valid   <- !is.na(renamed) & nchar(trimws(renamed)) > 0
        trajectory_map[valid] <- renamed[valid]
        message("Trajectories renamed successfully.")
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
  
  n_unassigned <- sum(cluster_df$Trajectory == "Unassigned")
  if (n_unassigned > 0) {
    warning(n_unassigned, " features could not be assigned a trajectory.")
  }
  
  # ----------------------------
  # Trajectory size summary
  # ----------------------------
  traj_size_df <- cluster_df %>%
    dplyr::group_by(Trajectory) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_features))
  
  message("\nTrajectory sizes:")
  print(as.data.frame(traj_size_df), row.names = FALSE)
  
  # ----------------------------
  # Simplified trajectory plot
  # ----------------------------
  data_long <- data_z %>%
    tidyr::pivot_longer(
      cols      = tidyr::all_of(tp_levels),
      names_to  = "Timepoint",
      values_to = "Z_Score"
    ) %>%
    dplyr::mutate(
      Timepoint  = factor(Timepoint, levels = tp_levels),
      Trajectory = factor(Trajectory)
    )
  
  n_traj <- length(unique(data_long$Trajectory))
  n_col  <- min(3, n_traj)
  
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
      ggplot2::aes(group = Trajectory, color = Trajectory),
      fun       = mean,
      geom      = "line",
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ Trajectory, ncol = n_col) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0("WATER Trajectories (", n_traj, " groups, ",
                     nrow(cluster_df), " features)"),
      x     = "",
      y     = "Z-score"
    ) +
    ggplot2::theme(
      legend.position    = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
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
    plot_dir  <- water_get_output_path(water, "label_classic", "plots")
    table_dir <- water_get_output_path(water, "label_classic", "tables")
    
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
      file.path(plot_dir, "trajectory_plot.png"),
      traj_plot,
      width  = 10,
      height = 6,
      dpi    = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, "trajectory_plot.pdf"),
      traj_plot,
      width  = 10,
      height = 6
    )
    
    # --------------------------
    # Save tables
    # --------------------------
    water_save_output(
      cluster_df,
      "water_trajectory_assignments",
      water,
      step = "label_classic"
    )
    
    water_save_output(
      traj_size_df,
      "water_trajectory_sizes",
      water,
      step = "label_classic"
    )
    
    traj_map_df <- data.frame(
      Cluster    = names(trajectory_map),
      Trajectory = unname(trajectory_map),
      stringsAsFactors = FALSE
    )
    
    water_save_output(
      traj_map_df,
      "water_trajectory_map",
      water,
      step = "label_classic"
    )
  }
  
  # ----------------------------
  # Store results
  # ----------------------------
  water$label_classic <- list(
    method           = method,
    r_threshold      = if (method == "centroid") r_threshold else NULL,
    trajectory_map   = trajectory_map,
    cluster_df       = cluster_df,
    data_z           = data_z,
    trajectory_sizes = traj_size_df,
    centroids        = centroids,
    trajectory_plot  = traj_plot
  )
  
  water$history[[length(water$history) + 1]] <- list(
    step           = "label_classic",
    method         = method,
    r_threshold    = if (method == "centroid") r_threshold else NULL,
    n_trajectories = length(unique(cluster_df$Trajectory)),
    timestamp      = Sys.time()
  )
  
  return(water)
}