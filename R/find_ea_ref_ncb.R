#' @title find_ea_ref_ncb
#' @description Compare your maf with 1000G as reference panel.
#'
#' @param dat A dataframe containing MAF columns.
#'
#' @return A dataframe with ea_ref column from 1000G.
#' @export
#'
#' @examples
#' dat <- find_ea_ref_ncb(dat)
#'
#'
find_ea_ref_ncb <- function(df, pop="European") {
  suppressMessages(require('rvest'))
  suppressMessages(require('dplyr'))

  tryCatch({
    df$sample_size <- NA
    df$ref_allele <- NA
    df$alt_allele <- NA

    for (i in 1:nrow(df)) {
      # Construct the URL for the SNP
      url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$SNP[i])

      # Read the HTML content of the SNP page
      page <- read_html(url)

      # Extract allele frequency data from the page
      allele_data <- page %>%
        html_nodes("table") %>%
        html_table(fill = TRUE) %>%
        .[[2]]  # Get the second table

      # Extract data for the European population
      ea_row <- subset(allele_data, Study == '1000Genomes' & Population == "Europe")


      if (nrow(ea_row) > 0) {
        # Extract relevant values
        sample_size <- ea_row$`Sample Size`
        ref_allele <- ea_row$`Ref Allele`
        alt_allele <- ea_row$`Alt Allele`

        # Fill the df dataframe with the extracted values
        df$sample_size[i] <- sample_size
        df$ref_allele[i] <- ref_allele
        df$alt_allele[i] <- alt_allele
      } else {
        message("No data found for SNP: ", df$SNP[i])
      }
    }
    return(df)
  }, error = function(e) {
    message("Error: ", e)
    return(NA)
  })
}
