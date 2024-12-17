#' Check EAF for Exposure and Outcome
#' @description This function checks the effect allele frequency (EAF) in both exposure and outcome data 
#' within a harmonized dataset. EAF cannot be added using reID because reID is only available in the reference panel. 
#' The effect allele (EA) in the reference panel may not necessarily be the EA for the reID in your dataset. 
#' This can only be addressed by comparing the EA and non-effect allele (NEA) in the reference panel 
#' versus the EA and NEA in your own dataset.
#' @param dat A harmonized dataset containing the genetic summary statistics.
#' @param pop Population for reference panel, default is "EUR".
#'
#' @return The input dataset with EAF information checked for both exposure and outcome. Messages are displayed for missing values.
#' @export
#'
#' @examples
#' # Assuming 'dat' is a harmonized dataset
#' check_eaf(dat, pop = "EUR")
#'
check_eaf <- function(dat, pop = "EUR") {
  required_cols <- c("eaf.exposure", "eaf.outcome")
  
  # Check if required columns exist in the dataset
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing in the dataset:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check EAF for exposure
  missing_exposure_eaf <- which(is.na(dat$eaf.exposure))
  if (length(missing_exposure_eaf) == 0) {
    message("\n\nAll SNPs in the exposure data have complete EAF. No action is needed.\n\n")
  } else {
    message(sprintf("%d SNPs in the exposure data are missing EAF values.\n", length(missing_exposure_eaf)))
    message("Consider using a different GWAS dataset or providing a reference panel for imputation.\n")
  }
  
  # Check EAF for outcome
  missing_outcome_eaf <- which(is.na(dat$eaf.outcome))
  if (length(missing_outcome_eaf) == 0) {
    message("\n\nAll SNPs in the outcome data have complete EAF. No action is needed.\n\n")
  } else {
    message(sprintf("%d SNPs in the outcome data are missing EAF values.\n", length(missing_outcome_eaf)))
    message("Consider using a different GWAS dataset or providing a reference panel for imputation.\n")
  }
  
  return(dat)
}
