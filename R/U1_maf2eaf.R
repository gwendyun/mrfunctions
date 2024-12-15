#' @title U1_maf2eaf
#' @description Convert MAF to EAF in a dataframe, checking for and renaming relevant columns.
#'
#' @param df A dataframe containing MAF or EAF columns.
#'
#' @return A dataframe with EAF column, preserving the original MAF column.
#' @export
#'
#' @examples
#' df <- U1_maf2eaf(df)
#'
#'
U1_maf2eaf <- function(df) {
  require(dplyr)
  # https://github.com/sinarueeger/GWAS.utils/blob/master/R/eaf2maf.R
  # Check for MAF columns and rename if present
  if (any(c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf") %in% colnames(df))) {
    df <- rename_col(df, patterns = c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf"), format = "maf")
  }

  # Check for EAF columns and rename if present
  if (any(c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency") %in% colnames(df))) {
    df <- rename_col(df, patterns = c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency"), format = "eaf")
  }

  # If no MAF column found, return message
  if (!"maf" %in% colnames(df)) {
    message('No MAF column found, please check if there is an EAF column available or this will be added by U4_add_eaf function')
    return(NULL)
  }

  # Ensure MAF column is numeric
  df$maf <- as.numeric(df$maf)

  # Calculate EAF from MAF
  df$eaf <- 1 - df$maf

  # Return the dataframe with the new EAF column, preserving the original MAF column
  return(df)
}
