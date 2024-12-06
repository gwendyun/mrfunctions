#' @title per doubling in genetic liability for xxx
#' @description https://arxiv.org/abs/1804.05545
#'
#' @param res
#'
#' @return log res
#' @export
#'
#' @examples res <- log_res(res)
#'
#'
log_res <- function(res) {
  #IVW
  res$`logb_Inverse variance weighted` <- 0.693 * res$`b_Inverse variance weighted`
  res$`loglo_ci_Inverse variance weighted` <- 0.693 * res$`lo_ci_Inverse variance weighted`
  res$`logup_ci_Inverse variance weighted` <- 0.693 * res$`up_ci_Inverse variance weighted`

  res$`logor_Inverse variance weighted` <- exp(res$`logb_Inverse variance weighted`)
  res$`logor_lci95_Inverse variance weighted` <- exp(res$`loglo_ci_Inverse variance weighted`)
  res$`logor_uci95_Inverse variance weighted` <- exp(res$`logup_ci_Inverse variance weighted`)

  res$`logBeta (95% CI)_Inverse variance weighted` <- sprintf("%.3f (%.3f to %.3f)",
                                                              res$`logb_Inverse variance weighted`,
                                                              res$`loglo_ci_Inverse variance weighted`,
                                                              res$`logup_ci_Inverse variance weighted`)

  res$`logOR (95% CI)_Inverse variance weighted` <- sprintf("%.3f (%.3f to %.3f)",
                                                            res$`logor_Inverse variance weighted`,
                                                            res$`logor_lci95_Inverse variance weighted`,
                                                            res$`logor_uci95_Inverse variance weighted`)


  res$`logb_MR Egger` <- 0.693 * res$`b_MR Egger`
  res$`loglo_ci_MR Egger` <- 0.693 * res$`lo_ci_MR Egger`
  res$`logup_ci_MR Egger` <- 0.693 * res$`up_ci_MR Egger`

  res$`logor_MR Egger` <- exp(res$`logb_MR Egger`)
  res$`logor_lci95_MR Egger` <- exp(res$`loglo_ci_MR Egger`)
  res$`logor_uci95_MR Egger` <- exp(res$`logup_ci_MR Egger`)

  res$`logBeta (95% CI)_MR Egger` <- sprintf("%.3f (%.3f to %.3f)",
                                             res$`logb_MR Egger`,
                                             res$`loglo_ci_MR Egger`,
                                             res$`logup_ci_MR Egger`)

  res$`logOR (95% CI)_MR Egger` <- sprintf("%.3f (%.3f to %.3f)",
                                           res$`logor_MR Egger`,
                                           res$`logor_lci95_MR Egger`,
                                           res$`logor_uci95_MR Egger`)


  res$`logb_Weighted median` <- 0.693 * res$`b_Weighted median`
  res$`loglo_ci_Weighted median` <- 0.693 * res$`lo_ci_Weighted median`
  res$`logup_ci_Weighted median` <- 0.693 * res$`up_ci_Weighted median`

  res$`logor_Weighted median` <- exp(res$`logb_Weighted median`)
  res$`logor_lci95_Weighted median` <- exp(res$`loglo_ci_Weighted median`)
  res$`logor_uci95_Weighted median` <- exp(res$`logup_ci_Weighted median`)

  res$`logBeta (95% CI)_Weighted median` <- sprintf("%.3f (%.3f to %.3f)",
                                                    res$`logb_Weighted median`,
                                                    res$`loglo_ci_Weighted median`,
                                                    res$`logup_ci_Weighted median`)

  res$`logOR (95% CI)_Weighted median` <- sprintf("%.3f (%.3f to %.3f)",
                                                  res$`logor_Weighted median`,
                                                  res$`logor_lci95_Weighted median`,
                                                  res$`logor_uci95_Weighted median`)


  res$`logb_Simple mode` <- 0.693 * res$`b_Simple mode`
  res$`loglo_ci_Simple mode` <- 0.693 * res$`lo_ci_Simple mode`
  res$`logup_ci_Simple mode` <- 0.693 * res$`up_ci_Simple mode`

  res$`logor_Simple mode` <- exp(res$`logb_Simple mode`)
  res$`logor_lci95_Simple mode` <- exp(res$`loglo_ci_Simple mode`)
  res$`logor_uci95_Simple mode` <- exp(res$`logup_ci_Simple mode`)

  res$`logBeta (95% CI)_Simple mode` <- sprintf("%.3f (%.3f to %.3f)",
                                                res$`logb_Simple mode`,
                                                res$`loglo_ci_Simple mode`,
                                                res$`logup_ci_Simple mode`)

  res$`logOR (95% CI)_Simple mode` <- sprintf("%.3f (%.3f to %.3f)",
                                              res$`logor_Simple mode`,
                                              res$`logor_lci95_Simple mode`,
                                              res$`logor_uci95_Simple mode`)


  res$`logb_Weighted mode` <- 0.693 * res$`b_Weighted mode`
  res$`loglo_ci_Weighted mode` <- 0.693 * res$`lo_ci_Weighted mode`
  res$`logup_ci_Weighted mode` <- 0.693 * res$`up_ci_Weighted mode`

  res$`logor_Weighted mode` <- exp(res$`logb_Weighted mode`)
  res$`logor_lci95_Weighted mode` <- exp(res$`loglo_ci_Weighted mode`)
  res$`logor_uci95_Weighted mode` <- exp(res$`logup_ci_Weighted mode`)

  res$`logBeta (95% CI)_Weighted mode` <- sprintf("%.3f (%.3f to %.3f)",
                                                  res$`logb_Weighted mode`,
                                                  res$`loglo_ci_Weighted mode`,
                                                  res$`logup_ci_Weighted mode`)

  res$`logOR (95% CI)_Weighted mode` <- sprintf("%.3f (%.3f to %.3f)",
                                                res$`logor_Weighted mode`,
                                                res$`logor_lci95_Weighted mode`,
                                                res$`logor_uci95_Weighted mode`)

  res$logb_MR_PRESSO <- 0.693 * res$b_MR_PRESSO
  res$loglo_ci_MR_PRESSO <- 0.693 * res$lo_ci_MR_PRESSO
  res$logup_ci_MR_PRESSO <- 0.693 * res$up_ci_MR_PRESSO
  
  res$logor_MR_PRESSO <- exp(res$logb_MR_PRESSO)
  res$logor_lci95_MR_PRESSO <- exp(res$loglo_ci_MR_PRESSO)
  res$logor_uci95_MR_PRESSO <- exp(res$logup_ci_MR_PRESSO)
  
  res$`logBeta (95% CI)_MR_PRESSO` <- sprintf("%.3f (%.3f to %.3f)",
                                               res$logb_MR_PRESSO,
                                               res$loglo_ci_MR_PRESSO,
                                               res$logup_ci_MR_PRESSO)
  
  res$`logOR (95% CI)_MR_PRESSO` <- sprintf("%.3f (%.3f to %.3f)",
                                             res$logor_MR_PRESSO,
                                             res$logor_lci95_MR_PRESSO,
                                             res$logor_uci95_MR_PRESSO)

library(data.table)
# 指定需要调整顺序的列名
cols_to_move <- c("outliers", 
                  "Beta (95% CI)_MR_PRESSO", 
                  "OR (95% CI)_MR_PRESSO", 
                  "logBeta (95% CI)_MR_PRESSO", 
                  "logOR (95% CI)_MR_PRESSO", 
                  "logb_MR_PRESSO", 
                  "loglo_ci_MR_PRESSO", 
                  "logup_ci_MR_PRESSO", 
                  "logor_MR_PRESSO", 
                  "logor_lci95_MR_PRESSO", 
                  "logor_uci95_MR_PRESSO", 
                  "PRESSO_Distortion_beta", 
                  "PRESSO_Distortion_Pvalue")

# 找出不存在的列，并创建这些列填充为 NA
missing_cols <- setdiff(cols_to_move, names(res))
if (length(missing_cols) > 0) {
  for (col in missing_cols) {
    res[[col]] <- NA
  }
}

# 调整列顺序，将目标列移到最后
setcolorder(res, c(setdiff(names(res), cols_to_move), cols_to_move))

  return(res)
}
