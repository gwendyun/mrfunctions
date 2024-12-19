#' @title U1_maf2eaf
#' @description Convert MAF to EAF in a dataframe, checking for and renaming relevant columns.
#'
#' @param dat A dataframe containing MAF or EAF columns.
#'
#' @return A dataframe with EAF column, preserving the original MAF column.
#' @export
#'
#' @examples
#' dat <- U1_maf2eaf(dat)
#'
#'
U1_maf2eaf <- function(dat, pop = "European") {
  suppressMessages(require('dplyr'))
  suppressMessages(require('progress'))
  
  rename_col <- function(df, patterns, format) {
    B <- colnames(df)
    for (i in 1:length(B)) {
      col <- B[i]
      for (y in 1:length(patterns)) {
        if (tolower(patterns[y]) == tolower(B[i])) {
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
  
  if (is.data.frame(dat)) { dat = as.data.frame(dat) }
  
  # Check for MAF columns and rename if present
  if (any(c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf") %in% colnames(dat))) {
    dat <- rename_col(dat, patterns = c("MinorAlleleFrequency", "minor_allele_frequency", "MAF", "Freq", "freq", "alt_freq", 'freq_minor', 'Freq_Minor_Allele', "maf"), format = "maf")
    dat <- rename_col(dat, patterns = c("snp", "rsid", "rsids", "snpid", "rnpid", "rs", "variant_id"), format = "SNP")
    
    # if maf exists, then add eaf_ref from EU 1000G
    if ('maf' %in% colnames(dat)) {
      message(paste0('Using 1000G ', pop, ' population as default.'))
      
      dat <- mrfunctions::find_ea_ref_ncb(dat, pop = pop)  # Pass the pop parameter here

      # If still NA, fill with 0.5
      index <- which(is.na(dat$eaf))
      if (length(index) > 0) {
        message(paste0("There are still ", length(index), " SNPs without EAF. Filling with 0.5."))
        dat$eaf_ref[index] <- 0.5
      } else {
        message("All SNPs in eaf_ref have been filled.")
      }
    }
  } else if (any(c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency") %in% colnames(dat))) {
    dat <- rename_col(dat, patterns = c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency"), format = "eaf")
  } else {
    message('No MAF and eaf column found. It is really important to have this. This cannot be added just using 1000G.')
    return(NULL)
  }
  
  return(dat)
}
