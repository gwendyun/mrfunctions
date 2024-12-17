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
U1_maf2eaf <- function(dat) {
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
      message('1000G use EU pop as default')
      
      # Function to fetch EAF from Ensembl for each SNP
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
      
      # Function to fill eaf_ref column using find_eaf for each SNP
      fill_eaf_ref <- function(dat, pop = "EUR") {
        # Find the rows where eaf_ref is NA
        index <- which(is.na(dat$eaf_ref))
        
        if (length(index) > 0) {
          message(paste0("There are ", length(index), " SNPs without EAF in eaf_ref."))
        }
        
        # Get the unique SNPs that need to be filled
        index_SNP <- unique(dat$SNP[index])
        
        # Create a progress bar for the query process
        pb <- progress::progress_bar$new(total = length(index_SNP))
        
        # Iterate through the SNPs and fill in the eaf_ref column
        for (i in 1:length(index_SNP)) {
          dat$eaf_ref[which(dat$SNP == index_SNP[i])] <- find_eaf(index_SNP[i], pop)
          pb$tick()
        }
        
        # If still NA after attempting to fill, fill with 0.5
        index <- which(is.na(dat$eaf_ref))
        if (length(index) > 0) {
          message(paste0("There are still ", length(index), " SNPs without EAF in eaf_ref. Filling with 0.5."))
          dat$eaf_ref[index] <- 0.5
        } else {
          message("All SNPs in eaf_ref have been filled.")
        }
        
        return(dat)
      }
      
      # Call the fill_eaf_ref function
      dat <- fill_eaf_ref(dat, pop = "EUR")
      
      # Handle assignment of eaf column
      if (all(dat$maf < 0.5, na.rm = TRUE) & all(dat$eaf_ref < 0.5, na.rm = TRUE)) {
        dat$eaf <- dat$maf
        message('Compare maf with eaf_ref: both are less than 0.5. eaf=maf')
      } else {
        message('Compare maf with eaf_ref: effect allele is a major allele, eaf = 1 - maf')
        dat$eaf <- 1 - dat$maf
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
