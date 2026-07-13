run_condition_comparison <- function(
    cond1_name,
    cond2_name,
    cond1_dir,
    cond2_dir,
    out_dir,
    mut_colors = NULL,
    run_GO = TRUE
) {
  
  message("\n========================================")
  message("WATER Condition Comparison")
  message("========================================")
  
  # ----------------------------
  # REQUIRED PACKAGES
  # ----------------------------
  pkgs <- c("dplyr", "tidyr", "readr", "ggplot2", "readxl", "tibble", "patchwork")
  if (run_GO) pkgs <- c(pkgs, "clusterProfiler", "org.Mm.eg.db")
  
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package required but not installed: ", p)
    }
  }
  
  `%>%` <- dplyr::`%>%`
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(mut_colors)) {
    message("No mut_colors provided -> using default colors")
    mut_colors <- tibble::tribble(
      ~Mut_Traj, ~spaghetti, ~mean,
      "default", "grey60", "black"
    )
  }
  
  # ----------------------------
  # HELPERS
  # ----------------------------
  safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  
  read_expression <- function(cond_dir) {
    file <- file.path(cond_dir, "normalize", "tables", "water_data_normalized_avg.xlsx")
    
    if (!file.exists(file)) {
      stop("Normalized averaged expression file not found: ", file)
    }
    
    df <- readxl::read_excel(file)
    
    expr <- df %>%
      dplyr::rename(feature_id = 1) %>%
      tibble::column_to_rownames("feature_id") %>%
      as.matrix()
    
    if (is.null(rownames(expr))) {
      stop("Expression matrix missing rownames")
    }
    
    expr
  }
  
  read_nonDA <- function(cond_dir) {
    file <- file.path(cond_dir, "filter", "tables", "water_is_non_da.xlsx")
    
    if (!file.exists(file)) {
      message("NonDA file not found: ", file)
      return(character())
    }
    
    df <- readxl::read_excel(file)
    colnames(df) <- gsub("[^a-z0-9_]", "", tolower(colnames(df)))
    
    if (!all(c("feature_id", "is_non_da") %in% colnames(df))) {
      stop("NonDA file must contain columns: Feature_ID, is_non_da")
    }
    
    df %>%
      dplyr::filter(is_non_da == TRUE) %>%
      dplyr::pull(feature_id)
  }
  
  zscore_matrix <- function(mat) {
    z <- t(scale(t(mat)))
    z[is.na(z)] <- 0
    z
  }
  
  run_GO_terms <- function(genes, out_dir, name_prefix) {
    if (!run_GO) return(NULL)
    
    OrgDb <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
    
    map <- tryCatch(
      clusterProfiler::bitr(
        unique(genes),
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = OrgDb
      ),
      error = function(e) NULL
    )
    
    if (is.null(map) || nrow(map) == 0) return(NULL)
    
    res <- list()
    
    for (ont in c("BP", "MF", "CC")) {
      ego <- tryCatch(
        clusterProfiler::enrichGO(
          gene = unique(map$ENTREZID),
          OrgDb = OrgDb,
          keyType = "ENTREZID",
          ont = ont,
          readable = TRUE
        ),
        error = function(e) NULL
      )
      
      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
      
      ego_s <- tryCatch(
        clusterProfiler::simplify(ego),
        error = function(e) ego
      )
      
      ego_df <- as.data.frame(ego_s)
      readr::write_csv(
        ego_df,
        file.path(out_dir, paste0(name_prefix, "_", ont, "_GO.csv"))
      )
      ego_df$Ontology <- ont
      res[[ont]] <- ego_df
    }
    
    if (length(res) == 0) return(NULL)
    dplyr::bind_rows(res)
  }
  
  make_wt_only_plot <- function(df_plot, wt_traj, n_genes, wt_color) {
    
    p_wt <- ggplot2::ggplot(
      dplyr::filter(df_plot, Group == cond1_name),
      ggplot2::aes(Timepoint, Zscore, group = feature_id)
    ) +
      ggplot2::geom_line(
        color = wt_color,
        alpha = 0.10,
        linewidth = 0.5
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "line",
        linewidth = 1.5,
        color = wt_color
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0(cond1_name, "\nExpression"),
        y = "z-score",
        x = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    
    p_mut <- ggplot2::ggplot(
      dplyr::filter(df_plot, Group == cond2_name),
      ggplot2::aes(Timepoint, Zscore, group = feature_id)
    ) +
      ggplot2::geom_line(
        color = "grey70",
        alpha = 0.10,
        linewidth = 0.5
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "line",
        linewidth = 1.5,
        color = "grey50"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0(cond2_name, "\nExpression"),
        y = "z-score",
        x = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    
    patchwork::wrap_plots(p_wt, p_mut, ncol = 2) +
      patchwork::plot_annotation(
        title = paste0(cond1_name, "-only ", wt_traj, " genes"),
        subtitle = paste("n =", n_genes, "genes")
      )
  }
  
  make_wt_only_dominant_shift_plot <- function(
    df_plot,
    wt_traj,
    dominant_shift,
    n_genes,
    n_shift,
    wt_color,
    mut_spaghetti,
    mut_mean
  ) {
    
    p_wt <- ggplot2::ggplot(
      dplyr::filter(df_plot, Group == cond1_name),
      ggplot2::aes(Timepoint, Zscore, group = feature_id)
    ) +
      ggplot2::geom_line(
        color = wt_color,
        alpha = 0.10,
        linewidth = 0.5
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "line",
        linewidth = 1.5,
        color = wt_color
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0(cond1_name, "\nExpression"),
        y = "z-score",
        x = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    
    p_mut <- ggplot2::ggplot(
      dplyr::filter(df_plot, Group == cond2_name),
      ggplot2::aes(Timepoint, Zscore, group = feature_id)
    ) +
      # all mutant genes in grey
      ggplot2::geom_line(
        color = "grey75",
        alpha = 0.10,
        linewidth = 0.5
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "line",
        linewidth = 1.2,
        color = "grey55"
      ) +
      # highlight dominant mutant shift only
      ggplot2::geom_line(
        data = dplyr::filter(df_plot, Group == cond2_name, Highlight),
        color = mut_spaghetti,
        alpha = 0.25,
        linewidth = 0.6
      ) +
      ggplot2::stat_summary(
        data = dplyr::filter(df_plot, Group == cond2_name, Highlight),
        fun = mean,
        geom = "line",
        linewidth = 1.5,
        color = mut_mean
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0(cond2_name, "\nExpression"),
        y = "z-score",
        x = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    
    patchwork::wrap_plots(p_wt, p_mut, ncol = 2) +
      patchwork::plot_annotation(
        title = paste0(cond1_name, "-only ", wt_traj, " genes"),
        subtitle = paste0(
          "n = ", n_genes, " genes | dominant shift: ",
          wt_traj, " -> ", dominant_shift,
          " (n = ", n_shift, ")"
        )
      )
  }
  
  # ----------------------------
  # LOAD DATA
  # ----------------------------
  message("\n[1/7] Loading WATER outputs")
  
  w1 <- readRDS(file.path(cond1_dir, "final", "objects", "water_final_results.rds"))
  w2 <- readRDS(file.path(cond2_dir, "final", "objects", "water_final_results.rds"))
  
  traj1 <- w1$trajectories %>%
    dplyr::select(feature_id, trajectory) %>%
    dplyr::mutate(trajectory = as.character(trajectory)) %>%
    dplyr::rename(Cond1_Traj = trajectory)
  
  traj2 <- w2$trajectories %>%
    dplyr::select(feature_id, trajectory) %>%
    dplyr::mutate(trajectory = as.character(trajectory)) %>%
    dplyr::rename(Cond2_Traj = trajectory)
  
  message("\nLoading normalized averaged expression")
  
  expr1 <- read_expression(cond1_dir)
  expr2 <- read_expression(cond2_dir)
  
  if (!identical(colnames(expr1), colnames(expr2))) {
    stop("Condition expression matrices do not share identical timepoint columns")
  }
  
  expr1_z <- zscore_matrix(expr1)
  expr2_z <- zscore_matrix(expr2)
  
  # ----------------------------
  # MASTER GENE LIST
  # ----------------------------
  message("\n[2/7] Building master gene list")
  
  all_genes <- Reduce(union, list(
    traj1$feature_id,
    traj2$feature_id,
    rownames(expr1),
    rownames(expr2)
  ))
  
  master_df <- data.frame(feature_id = all_genes, stringsAsFactors = FALSE) %>%
    dplyr::left_join(traj1, by = "feature_id") %>%
    dplyr::left_join(traj2, by = "feature_id")
  
  # ----------------------------
  # FILL MISSING
  # ----------------------------
  message("\n[3/7] Assigning missing states")
  
  nonDA1 <- read_nonDA(cond1_dir)
  nonDA2 <- read_nonDA(cond2_dir)
  
  master_df$Cond1_Traj <- ifelse(
    is.na(master_df$Cond1_Traj),
    ifelse(master_df$feature_id %in% nonDA1, "NonDA", "Not_Exp"),
    master_df$Cond1_Traj
  )
  
  master_df$Cond2_Traj <- ifelse(
    is.na(master_df$Cond2_Traj),
    ifelse(master_df$feature_id %in% nonDA2, "NonDA", "Not_Exp"),
    master_df$Cond2_Traj
  )
  
  # ----------------------------
  # TRAJECTORY OVERLAPS (Panel A)
  # ----------------------------
  message("\n[4/7] Generating trajectory overlaps")
  
  traj_levels <- sort(unique(c(master_df$Cond1_Traj, master_df$Cond2_Traj)))
  
  panelA_df <- lapply(traj_levels, function(traj_name) {
    data.frame(
      Trajectory = traj_name,
      Cond1_only = sum(master_df$Cond1_Traj == traj_name & master_df$Cond2_Traj != traj_name),
      Shared     = sum(master_df$Cond1_Traj == traj_name & master_df$Cond2_Traj == traj_name),
      Cond2_only = sum(master_df$Cond2_Traj == traj_name & master_df$Cond1_Traj != traj_name)
    )
  }) %>% dplyr::bind_rows()
  
  panelA_long <- panelA_df %>%
    tidyr::pivot_longer(
      -Trajectory,
      names_to = "Category",
      values_to = "Count"
    )
  
  pA <- ggplot2::ggplot(
    panelA_long,
    ggplot2::aes(Trajectory, Count, fill = Category)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      y = "Number of genes",
      x = NULL,
      fill = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    file.path(out_dir, "Trajectory_overlap_between_conditions.pdf"),
    pA,
    width = 10,
    height = 5
  )
  
  readr::write_csv(
    panelA_df,
    file.path(out_dir, "Trajectory_overlap_between_conditions.csv")
  )
  
  # ----------------------------
  # BUILD GLOBAL Z TABLE
  # ----------------------------
  message("\n[5/7] Building z-score table")
  
  combined_z <- dplyr::bind_rows(
    data.frame(feature_id = rownames(expr1_z), expr1_z, Group = cond1_name),
    data.frame(feature_id = rownames(expr2_z), expr2_z, Group = cond2_name)
  ) %>%
    tidyr::pivot_longer(
      -c(feature_id, Group),
      names_to = "Timepoint",
      values_to = "Zscore"
    ) %>%
    dplyr::mutate(
      Timepoint = factor(Timepoint, levels = colnames(expr1_z))
    ) %>%
    dplyr::left_join(master_df, by = "feature_id")
  
  transitions <- master_df %>%
    dplyr::count(Cond1_Traj, Cond2_Traj, name = "n")
  
  readr::write_csv(
    transitions,
    file.path(out_dir, "trajectory_transition_counts.csv")
  )
  
  # ----------------------------
  # TRANSITION PLOTS + GO
  # ----------------------------
  message("\n[6/7] Generating transition plots and GO")
  
  GO_MASTER <- list()
  WT_ONLY_GO_MASTER <- list()
  
  baseline_levels <- sort(unique(master_df$Cond1_Traj))
  baseline_levels <- baseline_levels[!baseline_levels %in% c("NonDA", "Not_Exp")]
  
  for (wt in baseline_levels) {
    
    wt_dir <- file.path(out_dir, safe(wt))
    dir.create(wt_dir, recursive = TRUE, showWarnings = FALSE)
    
    # --------------------------
    # WT-only plot for this baseline trajectory
    # --------------------------
    wt_only_genes <- master_df %>%
      dplyr::filter(Cond1_Traj == wt, Cond2_Traj != wt) %>%
      dplyr::pull(feature_id) %>%
      unique()
    
    if (length(wt_only_genes) > 0) {
      
      readr::write_csv(
        tibble::tibble(feature_id = wt_only_genes),
        file.path(wt_dir, "WT_only_genes.csv")
      )
      
      df_wt_only <- combined_z %>%
        dplyr::filter(feature_id %in% wt_only_genes)
      
      wt_cols <- mut_colors %>%
        dplyr::filter(Mut_Traj == wt)
      
      wt_color <- if (nrow(wt_cols) > 0) wt_cols$mean[1] else "#4daf4a"
      
      # ----------------------------------
      # Plot 1: WT-only genes in WT and MUT
      # ----------------------------------
      p_wt_only <- make_wt_only_plot(
        df_plot = df_wt_only,
        wt_traj = wt,
        n_genes = length(wt_only_genes),
        wt_color = wt_color
      )
      
      ggplot2::ggsave(
        file.path(wt_dir, "WT_only_plot.pdf"),
        p_wt_only,
        width = 8,
        height = 4.5
      )
      
      # ----------------------------------
      # Determine dominant mutant shift
      # exclude same trajectory, NonDA, Not_Exp
      # ----------------------------------
      dominant_shift_df <- master_df %>%
        dplyr::filter(
          Cond1_Traj == wt,
          Cond2_Traj != wt,
          !Cond2_Traj %in% c("NonDA", "Not_Exp")
        ) %>%
        dplyr::count(Cond2_Traj, sort = TRUE)
      
      if (nrow(dominant_shift_df) > 0) {
        
        dominant_shift <- dominant_shift_df$Cond2_Traj[1]
        dominant_n <- dominant_shift_df$n[1]
        
        dominant_genes <- master_df %>%
          dplyr::filter(
            Cond1_Traj == wt,
            Cond2_Traj == dominant_shift
          ) %>%
          dplyr::pull(feature_id) %>%
          unique()
        
        dominant_cols <- mut_colors %>%
          dplyr::filter(Mut_Traj == dominant_shift)
        
        if (nrow(dominant_cols) == 0) {
          dominant_cols <- mut_colors %>%
            dplyr::filter(Mut_Traj == "default")
          
          if (nrow(dominant_cols) == 0) {
            dominant_cols <- data.frame(
              spaghetti = "grey60",
              mean = "black"
            )
          }
        }
        
        df_wt_only_highlight <- df_wt_only %>%
          dplyr::mutate(
            Highlight = feature_id %in% dominant_genes
          )
        
        p_wt_only_shift <- make_wt_only_dominant_shift_plot(
          df_plot = df_wt_only_highlight,
          wt_traj = wt,
          dominant_shift = dominant_shift,
          n_genes = length(wt_only_genes),
          n_shift = dominant_n,
          wt_color = wt_color,
          mut_spaghetti = dominant_cols$spaghetti[1],
          mut_mean = dominant_cols$mean[1]
        )
        
        ggplot2::ggsave(
          file.path(wt_dir, "WT_only_dominant_shift_plot.pdf"),
          p_wt_only_shift,
          width = 8,
          height = 4.5
        )
        
        readr::write_csv(
          tibble::tibble(feature_id = dominant_genes),
          file.path(wt_dir, "WT_only_dominant_shift_genes.csv")
        )
      }
      
      go_wt_only <- run_GO_terms(
        wt_only_genes,
        wt_dir,
        paste0(safe(wt), "_WT_only")
      )
      
      if (!is.null(go_wt_only)) {
        go_wt_only$Cond1_Traj <- wt
        go_wt_only$Cond2_Traj <- "WT_only"
        WT_ONLY_GO_MASTER[[wt]] <- go_wt_only
      }
    }
    
    # --------------------------
    # Transition plots for this baseline trajectory
    # --------------------------
    mut_targets <- transitions %>%
      dplyr::filter(Cond1_Traj == wt) %>%
      dplyr::pull(Cond2_Traj) %>%
      unique()
    
    for (mut in mut_targets) {
      
      genes <- master_df %>%
        dplyr::filter(Cond1_Traj == wt, Cond2_Traj == mut) %>%
        dplyr::pull(feature_id) %>%
        unique()
      
      if (length(genes) == 0) next
      
      combo_dir <- file.path(wt_dir, paste0(safe(wt), "to", safe(mut)))
      dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
      
      readr::write_csv(
        tibble::tibble(feature_id = genes),
        file.path(combo_dir, "genes.csv")
      )
      
      cols <- mut_colors %>%
        dplyr::filter(Mut_Traj == mut)
      
      if (nrow(cols) == 0) {
        cols <- mut_colors %>%
          dplyr::filter(Mut_Traj == "default")
        
        if (nrow(cols) == 0) {
          cols <- data.frame(spaghetti = "grey60", mean = "grey20")
        }
      }
      
      df_plot <- combined_z %>%
        dplyr::filter(feature_id %in% genes)
      
      p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(Timepoint, Zscore, group = feature_id)
      ) +
        ggplot2::geom_line(
          data = dplyr::filter(df_plot, Group == cond1_name),
          ggplot2::aes(color = Group),
          alpha = 0.10,
          linewidth = 0.5
        ) +
        ggplot2::stat_summary(
          data = dplyr::filter(df_plot, Group == cond1_name),
          ggplot2::aes(group = Group),
          fun = mean,
          geom = "line",
          linewidth = 1.5,
          color = "grey60",
          alpha = 1
        ) +
        ggplot2::geom_line(
          data = dplyr::filter(df_plot, Group == cond2_name),
          ggplot2::aes(color = Group),
          alpha = 0.25,
          linewidth = 0.6
        ) +
        ggplot2::stat_summary(
          data = dplyr::filter(df_plot, Group == cond2_name),
          ggplot2::aes(group = Group),
          fun = mean,
          geom = "line",
          linewidth = 1.5,
          color = cols$mean[1],
          alpha = 1
        ) +
        ggplot2::scale_color_manual(
          values = c(
            setNames("grey70", cond1_name),
            setNames(cols$spaghetti[1], cond2_name)
          )
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste(wt, "->", mut),
          subtitle = paste("n =", length(genes)),
          y = "Z-score",
          x = NULL,
          color = NULL
        ) +
        ggplot2::theme(
          legend.position = "none",
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.y = ggplot2::element_blank()
        )
      
      ggplot2::ggsave(
        file.path(combo_dir, "trajectory_plot.pdf"),
        p,
        width = 4,
        height = 4
      )
      
      go_transition <- run_GO_terms(
        genes,
        combo_dir,
        paste0(safe(wt), "to", safe(mut))
      )
      
      if (!is.null(go_transition)) {
        go_transition$Cond1_Traj <- wt
        go_transition$Cond2_Traj <- mut
        GO_MASTER[[paste(wt, mut, sep = "_")]] <- go_transition
      }
    }
  }
  
  # ----------------------------
  # SAVE MASTER OUTPUTS
  # ----------------------------
  message("\n[7/7] Saving master tables")
  
  readr::write_csv(
    master_df,
    file.path(out_dir, "trajectory_comparison_table.csv")
  )
  
  if (length(GO_MASTER) > 0) {
    GO_ALL <- dplyr::bind_rows(GO_MASTER)
    readr::write_csv(
      GO_ALL,
      file.path(out_dir, "ALL_WT_MUT_GO_ENRICHMENT.csv")
    )
  }
  
  if (length(WT_ONLY_GO_MASTER) > 0) {
    GO_WT_ONLY <- dplyr::bind_rows(WT_ONLY_GO_MASTER)
    readr::write_csv(
      GO_WT_ONLY,
      file.path(out_dir, "ALL_WT_ONLY_GO_ENRICHMENT.csv")
    )
  }
  
  message("\n========================================")
  message("Condition comparison complete")
  message("========================================")
}