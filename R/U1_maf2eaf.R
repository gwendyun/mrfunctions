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
    
    # if maf exist, then add eaf_ref from EU 1000G 
    if ('maf' %in% colnames(df)) {

      message ('1000G use EU pop as default')

      find_eaf <- function(SNP, pop = "EUR") {
        server <- "http://rest.ensembl.org/variation/Homo_sapiens/%s?content-type=application/json;pops=1"
        server <- sprintf(server, SNP)
        res <- httr::GET(server)
        httr::stop_for_status(res)
        res_pop <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(res)))$populations
        eaf <- subset(res_pop, population == paste0("1000GENOMES:phase_3:", pop))
        if (nrow(eaf) == 0) {
          eaf.final <- NA
        } else {
          eaf.final <- eaf$frequency[1][[1]]
        }
        return(eaf.final)
      }

      df$eaf_ref <- sapply(df$SNP, find_eaf)

      
      if (all(df$maf < 0.5, na.rm = TRUE) & all(df$eaf_ref < 0.5, na.rm = TRUE)) {
        df$eaf <- df$maf
        message('compare maf with eaf_ref. both are less than 0.5. eaf=maf')
      } else {
        message('compare maf with eaf_ref. effect_allele is a major allele, df$eaf = 1 - df$maf')
        return(df)
      }
    }
  } else if (any(c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency") %in% colnames(df))) {
    df <- rename_col(df, patterns = c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency"), format = "eaf")
  } else {
    # no eaf and maf 
    message('No MAF and eaf column found. it is really important to have this. this can not be added just using 1000G')
    return(NULL)
  }

  # Return the dataframe with the new EAF column, preserving the original MAF column
  return(df)
}
