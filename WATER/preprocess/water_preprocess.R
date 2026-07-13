water_preprocess <- function(
    data,
    replicates,
    log2fc_data = NULL,
    modality = "unknown",
    output_dir = NULL,
    save_outputs = TRUE,
    tpm_threshold = 1,
    feature_list = NULL,
    log2fc_threshold = 1,
    average_replicates = TRUE
) {
  
  # Load required WATER functions
  source("preprocess/water_get_output_path.R")
  source("preprocess/water_save_output.R")
  
  source("preprocess/water_initialize.R")
  source("preprocess/water_filter_expression.R")   
  source("preprocess/water_filter_features.R")
  source("preprocess/water_normalize.R")
  source("preprocess/water_PCA.R")
  
  # ----------------------------
  # Step 1: Initialize
  # ----------------------------
  message("\n  [1.1] Initialize")
  water <- water_initialize(
    data = data,
    replicates = replicates,
    log2fc_data = log2fc_data,
    modality = modality,
    output_dir = output_dir,
    save_outputs = save_outputs
  )
  
  # ----------------------------
  # Step 2: Expression filtering
  # ----------------------------
  message("\n  [1.2] Expression filtering")
  water <- water_filter_expression(
    water,
    tpm_threshold = tpm_threshold
  )
  
  # ----------------------------
  # Step 3: Feature processing
  # ----------------------------
  message("\n  [1.3] Feature filtering")
  water <- water_filter_features(
    water,
    feature_list = feature_list,
    log2fc_threshold = log2fc_threshold,
    average_replicates = average_replicates
  )
  
  # ----------------------------
  # Step 4: Normalize (replicates + averaged)
  # ----------------------------
  message("\n  [1.4] Normalization")
  water <- water_normalize(water)
  
  # ----------------------------
  # Step 5: PCA (replicates)
  # ----------------------------
  message("\n  [1.5] PCA")
  water <- water_PCA(water)
  
  message("\nPreprocessing complete")
  
  return(water)
}