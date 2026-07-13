run_water <- function(
    data,
    replicates,
    log2fc = NULL,
    preprocess_params = list(),
    r_threshold = 0.8,
    run_refine = TRUE,
    run_loo = TRUE
) {
  
  message("\n========================================")
  message("Running WATER")
  message("========================================")
  
  default_preprocess <- list(
    modality = "RNA",
    output_dir = "water_output",
    save_outputs = TRUE,
    tpm_threshold = 1,
    feature_list = NULL
  )
  
  preprocess_args <- modifyList(default_preprocess, preprocess_params)
  
  # ----------------------------
  # Load modules
  # ----------------------------
  source("preprocess/water_preprocess.R")
  source("clustering/water_cluster_classic.R")
  source("audit/water_audit_classic.R")
  source("refine/water_refine_classic.R")
  
  # ----------------------------
  # STEP 1: Preprocess
  # ----------------------------
  message("\n[1/4] Preprocessing")
  message("----------------------------------------")
  
  water <- do.call(
    water_preprocess,
    c(
      list(
        data = data,
        replicates = replicates,
        log2fc_data = log2fc
      ),
      preprocess_args
    )
  )
  
  # ----------------------------
  # STEP 2: Clustering
  # ----------------------------
  message("\n[2/4] Clustering")
  message("----------------------------------------")
  
  water <- water_cluster_classic(water)
  
  # ----------------------------
  # STEP 3: Audit
  # ----------------------------
  message("\n[3/4] Audit")
  message("----------------------------------------")
  
  water <- water_audit_classic(water)
  
  # ----------------------------
  # STEP 4: Refinement
  # ----------------------------
  message("\n[4/4] Refinement")
  message("----------------------------------------")
  
  if (run_refine) {
    water <- water_refine_classic(
      water,
      r_threshold = r_threshold,
      run_loo = run_loo
    )
  } else {
    message("Skipping refinement")
  }
  
  # ----------------------------
  # FINAL OUTPUTS
  # ----------------------------
  message("\n[Final] Preparing outputs")
  message("----------------------------------------")
  
  # ✅ USE FINAL TRAJECTORIES (FIXED)
  final_df <- water$correlate_secondary$final_assignments
  
  if (is.null(final_df)) {
    stop("Final trajectory assignments not found")
  }
  
  # ----------------------------
  # Merge expression
  # ----------------------------
  expr <- water$data_avg_norm
  
  expr_df <- data.frame(
    feature_id = rownames(expr),
    expr,
    stringsAsFactors = FALSE
  )
  
  final_table <- merge(
    final_df,
    expr_df,
    by = "feature_id",
    all.x = TRUE
  )
  
  rownames(final_table) <- NULL
  
  # ----------------------------
  # Save final table
  # ----------------------------
  table_dir <- water_get_output_path(water, "final", "tables")
  
  write.csv(
    final_table,
    file.path(table_dir, "water_final_gene_table.csv"),
    row.names = FALSE
  )
  
  message(sprintf("  Table: %s",
                  file.path(table_dir, "water_final_gene_table.csv")))
  
  # ----------------------------
  # Save final plot
  # ----------------------------
  plot_dir <- water_get_output_path(water, "final", "plots")
  
  traj_plot <- water$correlate_secondary$traj_plot_final
  
  if (is.null(traj_plot)) {
    stop("Final trajectory plot not found")
  }
  
  ggplot2::ggsave(
    file.path(plot_dir, "water_final_trajectories.png"),
    traj_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(plot_dir, "water_final_trajectories.pdf"),
    traj_plot,
    width = 12,
    height = 8
  )
  
  message(sprintf("  Plot: %s",
                  file.path(plot_dir, "water_final_trajectories.png")))
  
  # ----------------
  # Save final RDS 
  # ----------------
  obj_dir <- water_get_output_path(water, "final", "objects")
  
  water_final <- list(
    trajectories = final_df,
    timepoints   = water$timepoint_levels
  )
  
  t0 <- Sys.time()
  
  saveRDS(
    water_final,
    file.path(obj_dir, "water_final_results.rds"),
    compress = FALSE   # 🔥 critical for speed
  )
  
  t1 <- Sys.time()
  
  message(sprintf("  RDS saved (%.2f sec)",
                  as.numeric(difftime(t1, t0, units = "secs"))))
  
  message("\n========================================")
  message("WATER COMPLETE")
  message("========================================")
  
  return(water)
}