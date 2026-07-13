water_refine_classic <- function(
    water,
    r_threshold = 0.8,
    run_loo = TRUE
) {
  
  # ----------------------------
  # Load required functions
  # ----------------------------
  source("refine/water_recluster_ambiguous.R")
  source("refine/water_recluster_weak.R")
  source("refine/water_correlate_secondary.R")
  source("refine/water_secondary_LOO.R")
  
  # ----------------------------
  # Step 1: Reclustering ambiguous
  # ----------------------------
  message("  [4.1] Reclustering ambiguous features")
  water <- water_recluster_ambiguous(water)
  
  # ----------------------------
  # Step 2: Reclustering weak
  # ----------------------------
  message("  [4.2] Reclustering weak features")
  water <- water_recluster_weak(water)
  
  # ----------------------------
  # Step 3: Secondary correlation
  # ----------------------------
  message("  [4.3] Secondary correlation validation")
  water <- water_correlate_secondary(
    water,
    r_threshold = r_threshold
  )
  
  # ----------------------------
  # Step 4: Secondary LOO (optional)
  # ----------------------------
  if (run_loo) {
    message("  [4.4] Secondary LOO stability analysis")
    water <- water_secondary_LOO(
      water,
      r_threshold = r_threshold
    )
  } else {
    message("\nSkipping LOO analysis (run_loo = FALSE)")
  }
  
  return(water)
}