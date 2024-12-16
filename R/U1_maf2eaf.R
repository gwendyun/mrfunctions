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

  rename_col <- function(df ,
                         patterns ,
                         format ){
    B <- colnames( df )

    for(i in 1:length(B)){
      col <- B[i]
      for( y in 1:length(patterns) ){
        if( tolower(patterns[y]) == tolower(B[i])){

          if (!B[i] == format) {
            cat("Replacing", B[i], "with", format, "\n")
          }

          B[i] <- format
        }
      }
    }

    colnames(df) <- B

    if (length(unique(colnames(df))) != length(colnames(df))) {
      stop(paste(format, "column name is duplicated, please use `help(U1_Clean_data)` to check the documentation and provide only one!"), call. = FALSE)
    }

    return(df)
  }

  if(is.data.frame(df)){ df = as.data.frame(df) }

  # Check for MAF columns and rename if present
  if (any(c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf") %in% colnames(df))) {
    df <- rename_col(df, patterns = c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf"), format = "maf")
  }

  # Check for EAF columns and rename if present
  if (any(c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency") %in% colnames(df))) {
    df <- rename_col(df, patterns = c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency"), format = "eaf")
  }

  if (!"maf" %in% colnames(df)) {
    message('No MAF column found, please check if there is an EAF column available or this will be added by U4_add_eaf function')
    return(NULL)
  }

  if(is.data.frame(df)){ df = as.data.frame(df) }

  df$maf <- as.numeric(df$maf)

  if ("eaf" %in% colnames(df)) {
    df$eaf <- ifelse(df$eaf > 0.5, df$eaf, ifelse(df$maf < 0.5, df$maf, 1 - df$maf))
  } else {
    df$eaf <- ifelse(df$maf < 0.5, df$maf, 1 - df$maf)
  }

  # Return the dataframe with the new EAF column, preserving the original MAF column
  return(df)
}
