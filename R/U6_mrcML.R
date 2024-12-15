#' @title U6_mrcML
#' @description further analysis for univariable cML from MendelianRandomization package
#'
#' @param dat
#'
#' @return mrcML_res
#' @export
#'
#' @examples mrcML_res <- dat %>%
#'  dat_to_mr_input() %>%
#'  U6_mrcML(dat = ., random_seed = 314, Alpha = 0.05, num_pert = 200, maxit = 100)
#'
#'
U6_mrcML <- function(mr_input_obj, dat, random_seed = 314, Alpha = 0.05, num_pert = 200, maxit = 100) {
  require('MendelianRandomization')

  mrcML <- mr_cML(
    object = mr_input_obj,
    MA = TRUE,
    DP = TRUE,
    K_vec = 0:(length(mr_input_obj@betaX) - 2),
    random_start = 0,
    num_pert = num_pert,
    random_start_pert = 0,
    maxit = maxit,
    random_seed = random_seed,
    n = min(c(dat$samplesize.exposure, dat$samplesize.outcome), na.rm = TRUE),
    Alpha = Alpha
  )

  # Convert results into a data frame
  mrcML_res <- data.frame(
    Exposure = mrcML@Exposure,
    Outcome = mrcML@Outcome,
    Estimate = mrcML@Estimate,
    StdError = mrcML@StdError,
    Pvalue = mrcML@Pvalue,
    CI_Lower = mrcML@CILower,
    CI_Upper = mrcML@CIUpper,
    GOF1_p = mrcML@GOF1_p,
    GOF2_p = mrcML@GOF2_p,
    SNPs = mrcML@SNPs,
    Alpha = mrcML@Alpha,
    MA = mrcML@MA,
    DP = mrcML@DP
  )

  return(mrcML_res)
}

