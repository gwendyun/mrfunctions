#' @title dat_to_mr_input
#' @description TwoSampleMR package transfer to mrinput from MendelianRandomization package
#'
#' @param dat
#'
#' @return mr_input_obj
#' @export
#'
#' @examples mr_input_obj <- dat_to_mr_input(dat)
#'
#'
dat_to_mr_input <- function(dat) {
  suppressMessages(require('MendelianRandomization'))

  if (!all(c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome", "SNP") %in% colnames(dat))) {
    stop("Data does not have the required columns for MRInput creation.")
  }

  mr_input_obj <- mr_input(
    bx = as.numeric(dat[["beta.exposure"]]),
    bxse = as.numeric(dat[["se.exposure"]]),
    by = as.numeric(dat[["beta.outcome"]]),
    byse = as.numeric(dat[["se.outcome"]]),
    exposure = unique(dat[["exposure"]])[1], 
    outcome = unique(dat[["outcome"]])[1],   
    snps = as.character(dat[["SNP"]]),
    effect_allele = as.character(dat[["effect_allele.exposure"]]),
    other_allele = as.character(dat[["other_allele.exposure"]]),
    eaf = as.numeric(dat[["eaf.exposure"]])
  )

  return(mr_input_obj)
}
