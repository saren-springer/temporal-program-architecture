#' Run WATER audit pipeline
#'
#' Runs the initial LOO stability analysis followed by the initial
#' correlation analysis. Together these assess the quality of initial
#' trajectory assignments and identify features requiring refinement.
#'
#' @param water WATER object
#' @param r_threshold Numeric; correlation threshold for classification.
#'   Default 0.8.
#' @param seed Integer; random seed for LOO reproducibility. Default 7.
#'
#' @return Updated WATER object with audit results
#' @export

water_audit_classic <- function(
    water,
    r_threshold = 0.8,
    seed        = 7
) {
  
  # Load required WATER functions
  source("audit/water_initial_LOO.R")
  source("audit/water_correlate_initial.R")
  
  # ----------------------------
  # Step 1: Initial LOO
  # ----------------------------
  message("  [3.1] LOO stability analysis")
  water <- water_loo_initial(water, seed = seed)
  
  # ----------------------------
  # Step 2: Initial correlation
  # ----------------------------
  message("  [3.2] Correlation-based audit")
  water <- water_correlate_initial(water, r_threshold = r_threshold)
  
  # ----------------------------
  # Audit summary
  # ----------------------------
  message("  [3.3] Audit summary")
  
  loo_summary <- water$loo_initial$summary
  decision_counts <- water$correlate_initial$decision_counts
  
  mean_ari <- mean(loo_summary$ARI_mapped, na.rm = TRUE)
  mean_cluster_ret <- mean(loo_summary$Retention_cluster, na.rm = TRUE)
  mean_traj_ret <- mean(loo_summary$Retention_traj, na.rm = TRUE)
  
  get_val <- function(label) {
    decision_counts$n_features[decision_counts$decision == label]
  }
  
  get_pct <- function(label) {
    decision_counts$percent[decision_counts$decision == label]
  }
  
  message(sprintf("    LOO ARI (mean): %.3f", mean_ari))
  message(sprintf("    Cluster retention: %.1f%%", 100 * mean_cluster_ret))
  message(sprintf("    Trajectory retention: %.1f%%", 100 * mean_traj_ret))
  
  message("")
  message(sprintf("    Keep original: %d (%.1f%%)",
                  get_val("Keep_original"),
                  get_pct("Keep_original")))
  
  message(sprintf("    Reassigned:    %d (%.1f%%)",
                  get_val("Assign_new"),
                  get_pct("Assign_new")))
  
  message(sprintf("    Weak:          %d (%.1f%%)",
                  length(water$correlate_initial$weak_features),
                  get_pct("Weak")))
  
  message(sprintf("    Ambiguous:     %d (%.1f%%)",
                  length(water$correlate_initial$ambiguous_features),
                  get_pct("Ambiguous")))
  
  return(water)
}