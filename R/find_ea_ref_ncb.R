#' @title find_ea_ref_ncb
#' @description Compare your maf with 1000G as reference panel.
#'
#' @param df A dataframe containing MAF columns.
#' @param pop The population to use for EAF reference, default is "European".
#'
#' @return A dataframe with ea_ref column from 1000G.
#' @export
#'
#' @examples
#' dat <- find_ea_ref_ncb(dat)
#'
#'
find_ea_ref_ncb <- function(df, pop = "European") {
  suppressMessages(require('rvest'))
  suppressMessages(require('dplyr'))

  tryCatch({
    df$sample_size <- NA
    df$ref_allele <- NA
    df$alt_allele <- NA
    df$eaf <- NA
    df$alt_base <- NA
    df$alt_freq <- NA

    for (i in 1:nrow(df)) {
      # Construct the URL for the SNP
      url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$SNP[i])

      # Read the HTML content of the SNP page
      page <- read_html(url)

      # Extract allele frequency data from the page
      allele_data <- tryCatch({
        page %>%
          html_nodes("table") %>%
          html_table(fill = TRUE) %>%
          .[[2]]  # Get the second table
      }, error = function(e) {
        message("Error: ", e)
        df$ref_allele[i] <- NA
        df$alt_allele[i] <- NA
        return(NULL)
      })

      if (is.null(allele_data)) next

      # Extract data for the European population
      ea_row <- subset(allele_data, Study == '1000Genomes' & Population == "Europe")

      if (nrow(ea_row) > 0) {
        # Extract relevant values
        ref_allele <- ea_row$`Ref Allele`
        alt_allele <- ea_row$`Alt Allele`
        ea_row <- ea_row %>% 
          mutate(
            alt_base = sub("=.*", "", `Alt Allele`), 
            alt_freq = as.numeric(sub(".*=", "", `Alt Allele`))
          )
        # Fill the df dataframe with the extracted values
        df$ref_allele[i] <- ref_allele
        df$alt_allele[i] <- alt_allele
        df$alt_base[i] <- ea_row$alt_base
        df$alt_freq[i] <- ea_row$alt_freq

        # Check if effect_allele matches alt_allele
        if (df$effect_allele[i] == df$alt_base[i]) {
          if (df$alt_freq[i] < 0.5 && df$maf[i] < 0.5) {
            df$eaf[i] <- df$maf[i]
          } else if (df$alt_freq[i] > 0.5) {
            df$eaf[i] <- 1 - df$maf[i]
          } else if (df$alt_freq[i] < 0.5) {
            df$eaf[i] <- df$maf[i]
          }
        } else {
          df$eaf[i] <- 1 - df$alt_freq[i]
          if (df$eaf[i] < 0.5 && df$maf[i] < 0.5) {
            df$eaf[i] <- df$maf[i]
          } else if (df$alt_freq[i] > 0.5) {
            df$eaf[i] <- 1 - df$maf[i]
          } else if (df$alt_freq[i] < 0.5) {
            df$eaf[i] <- df$maf[i]
          }
        }
      } else {
        message("No data found for SNP: ", df$SNP[i])
        df$ref_allele[i] <- NA
        df$alt_allele[i] <- NA
      }
    }
    return(df)
  }, error = function(e) {
    message("Error: ", e)
    return(NA)
  })
}