#' @title U6_mrcML
#' @description further analysis for univariable cML from MendelianRandomization package
#'
#' @param dat
#'
#' @return mrcML_res
#' @export
#'
#' @examples mrcML_res <- U6_mrcML(dat, res, random_seed = 314, Alpha = 0.05, num_pert = 200, maxit = 100)
#'
#'
U6_mrcML <- function(dat, res, random_seed = 314, Alpha = 0.05, num_pert = 200, maxit = 100) {
  suppressMessages(require('MendelianRandomization'))
  suppressMessages(require('mrfunctions'))
  
  mr_input_obj <- mrfunctions::dat_to_mr_input(dat)

  # Run the mr_cML function
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
    exposure_cML = mrcML@Exposure,
    outcome_cML = mrcML@Outcome,
    b_cML = mrcML@Estimate,
    se_cML = mrcML@StdError,
    Pvalue = mrcML@Pvalue,
    lo_ci_cML = mrcML@CILower,
    up_ci_cML = mrcML@CIUpper,
    GOF1_p = mrcML@GOF1_p,
    GOF2_p = mrcML@GOF2_p,
    nSNP_cML = mrcML@SNPs,
    Alpha = mrcML@Alpha,
    MA = mrcML@MA,
    DP = mrcML@DP
  )

  # Calculate OR and log-transformed OR
  mrcML_res$or_cML <- exp(mrcML_res$b_cML)
  mrcML_res$or_lci95_cML <- exp(mrcML_res$lo_ci_cML)
  mrcML_res$or_uci95_cML <- exp(mrcML_res$up_ci_cML)

  mrcML_res$`logb_cML` <- 0.693 * mrcML_res$`b_cML`
  mrcML_res$`loglo_ci_cML` <- 0.693 * mrcML_res$`lo_ci_cML`
  mrcML_res$`logup_ci_cML` <- 0.693 * mrcML_res$`up_ci_cML`

  mrcML_res$`logor_cML` <- exp(mrcML_res$`logb_cML`)
  mrcML_res$`logor_lci95_cML` <- exp(mrcML_res$`loglo_ci_cML`)
  mrcML_res$`logor_uci95_cML` <- exp(mrcML_res$`logup_ci_cML`)

  mrcML_res$`logBeta (95% CI)_cML` <- sprintf("%.3f (%.3f to %.3f)",
                                        mrcML_res$`logb_cML`,
                                        mrcML_res$`loglo_ci_cML`,
                                        mrcML_res$`logup_ci_cML`)

  mrcML_res$`logOR (95% CI)_cML` <- sprintf("%.3f (%.3f to %.3f)",
                                      mrcML_res$`logor_cML`,
                                      mrcML_res$`logor_lci95_cML`,
                                      mrcML_res$`logor_uci95_cML`)

  # Merge results into the existing res data frame
  res <- merge(
    res, 
    mrcML_res, 
    by.x = c("exposure_Inverse variance weighted", "outcome_Inverse variance weighted"), # Columns in 'res'
    by.y = c("exposure_cML", "outcome_cML"), # Columns in 'mrcML_res'
    all.x = TRUE # Retain all rows from 'res'
  )

  return(res)
}
