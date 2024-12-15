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
  
  rename_col <- function(df ,
                         patterns ,
                         format ){
    B <- colnames( df )
    
    for(i in 1:length(B)){
      col <- B[i]
      for( y in 1:length(patterns) ){
        if( tolower(patterns[y]) == tolower(B[i])){
          
          if(!B[i]==format){cat("将",B[i],"替换为", format,"\n")}
          
          B[i] <- format
          
        }
      }
    }
    
    colnames( df ) <-  B
    
    if(length( unique(colnames(df) ) ) != length(colnames( df))){
      stop(paste(format,"列名重复，请用`help(U1_Clean_data)`查看说明文档，只提供一个!\n"), call. = FALSE)
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