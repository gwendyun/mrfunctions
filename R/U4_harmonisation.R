#' U4_harmonisation merges exposure and outcome data into dat
#'
#' @param exposure_dat  Instrumental variables extracted by U2_extract_instruments
#' @param outcome_dat  Outcomes extracted by U3_extract_outcomes_data
#' @param action Default is 2, refer to the description of harmonise_data in TwoSampleMR
#'
#' @return Harmonised data
#' @export
#'
#' @examples
#' dat <- U4_harmonisation(IVs, Outs)
#'
#' \dontrun{
#'   IVs <- U2_extract_instruments(c("ieu-a-31", "ieu-a-300", "ieu-a-2"))
#'   Outs <- U3_extract_outcomes_data(
#'     outcome = c("ieu-a-297", "ieu-a-298", "ieu-a-299"),
#'     exposure_iv = IVs
#'   )
#'   dat <- U4_harmonisation(IVs, Outs) %>%
#'     U4_add_eaf() %>%
#'     U4_check_samplesize()
#'
#'
#' #### Check if the sample size is correct ######################################
#'
#' # If there are too many IDs and they are not fully displayed, you can use View(dat) to check "description_exposure" and "description_outcome"
#'
#' #### Exposure part ###########################################################
#'
#' # If the sample size of Inflammatory bowel disease is incorrect, it needs to be modified
#'
#' # Modify total sample size
#' dat$samplesize.exposure[which(dat$exposure == "Inflammatory bowel disease")] <- 12000
#'
#' # Modify ncase
#' dat$ncase.exposure[which(dat$exposure == "Inflammatory bowel disease")] <- 1200
#'
#' # Modify ncontrol
#' dat$ncontrol.exposure[which(dat$exposure == "Inflammatory bowel disease")] <- 10800
#'
#' # Modify name
#' dat$exposure[which(dat$exposure == "Inflammatory bowel disease")] <- "IBD"
#'
#' #### Outcome part ############################################################
#'
#' # If HDL cholesterol is actually a binary variable, it needs to be modified
#'
#' # Add ncase
#' dat$ncase.outcome[which(dat$outcome == "HDL cholesterol")] <- 1200
#' # Add ncontrol
#' dat$ncontrol.outcome[which(dat$outcome == "HDL cholesterol")] <- 10800
#'
#'
#' # For the outcome name Alzheimer's disease, if you want to modify the first Alzheimer's disease
#'
#' dat$outcome[which(dat$id.outcome == "ieu-a-297")] <- "AD"
#'
#'
#' #### Final confirmation ######################################################
#'
#' dat <- U4_check_samplesize(dat)
#'
#' ###### If you want to remove SNPs with F value less than 10 ##################
#'
#' dat <- subset(dat, singleF > 10)
#' print(subset(dat, singleF <= 10)$SNP)
#'
#' ###### steiger_filtering #####################################################
#'
#' dat <- dat %>%
#'   TwoSampleMR::steiger_filtering() %>%
#'   {removed_snps <- dplyr::filter(., steiger_dir != TRUE) %>% dplyr::pull(SNP);
#'   cat("Removed SNPs:\n", removed_snps, "\n\nTotal number of removed SNPs:", length(removed_snps), "\n");
#'   dplyr::filter(., steiger_dir == TRUE)}
#'
#'}
U4_harmonisation <- function(exposure_dat,
                             outcome_dat,
                             action = 2){

  message("Harmonising...")
  dat <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                     outcome_dat = outcome_dat,
                                     action = action)
  dat <- plyr::ddply(dat,
                     c("id.exposure", "id.outcome"),
                     function(x1) {
                       x <- subset(x1, mr_keep) %>%
                         dplyr::distinct(SNP, .keep_all = TRUE)
                     })
  dat <- plyr::ddply(dat, "id.exposure",
                     function(dat) {

                       if (length(unique(dat$samplesize.exposure)) > 1) {
                         message("Exposure", dat$exposure[1], "has multiple different sample sizes, only the largest sample size will be used for subsequent calculations")
                       }

                       if ("ncase.exposure" %in% colnames(dat) & "ncontrol.exposure" %in% colnames(dat)) {

                         dat$ncase.exposure <- as.numeric(dat$ncase.exposure)
                         dat$ncontrol.exposure <- as.numeric(dat$ncontrol.exposure)

                         samplesize.exposure <- as.numeric(max(dat$samplesize.exposure))

                         if (all(!is.na(dat$ncase.exposure)) & all(!is.na(dat$ncontrol.exposure))) {

                           dat$type.exposure <- 'Binary'
                           dat$ncase.exposure <- dat$ncase.exposure[which(dat$samplesize.exposure == samplesize.exposure)[1]]
                           dat$ncontrol.exposure <- dat$ncontrol.exposure[which(dat$samplesize.exposure == samplesize.exposure)[1]]
                           dat$samplesize.exposure <- as.numeric(max(dat$samplesize.exposure))

                         } else {
                           dat$type.exposure <- "Continuous"
                           dat$samplesize.exposure <- as.numeric(max(dat$samplesize.exposure))
                         }

                       } else {
                         dat$samplesize.exposure <- as.numeric(max(dat$samplesize.exposure))
                         dat$type.exposure <- "Continuous"
                       }

                       return(dat)

                     })

  dat <- plyr::ddply(dat, c("id.outcome"),
                     function(dat) {

                       if (length(unique(dat$samplesize.outcome)) > 1) {
                         message("Outcome", dat$outcome[1], "has multiple different sample sizes, only the largest sample size will be used for subsequent calculations")
                       }

                       if ("ncase.outcome" %in% colnames(dat) & "ncontrol.outcome" %in% colnames(dat)) {

                         dat$ncase.outcome <- as.numeric(dat$ncase.outcome)
                         dat$ncontrol.outcome <- as.numeric(dat$ncontrol.outcome)

                         samplesize.outcome <- as.numeric(max(dat$samplesize.outcome))

                         if (all(!is.na(dat$ncase.outcome)) & all(!is.na(dat$ncontrol.outcome))) {

                           dat$type.outcome <- 'Binary'
                           dat$ncase.outcome <- dat$ncase.outcome[which(dat$samplesize.outcome == samplesize.outcome)[1]]
                           dat$ncontrol.outcome <- dat$ncontrol.outcome[which(dat$samplesize.outcome == samplesize.outcome)[1]]
                           dat$ratio <- dat$ncase.outcome / dat$ncontrol.outcome
                           dat$samplesize.outcome <- as.numeric(max(dat$samplesize.outcome))

                         } else {
                           dat$type.outcome <- "Continuous"
                           dat$samplesize.outcome <- as.numeric(max(dat$samplesize.outcome))
                         }

                       } else {
                         dat$samplesize.outcome <- as.numeric(max(dat$samplesize.outcome))
                         dat$type.outcome <- "Continuous"
                       }
                       return(dat)
                     })

  if (any(is.na(dat$samplesize.exposure))) { warning("There are NA values in the total sample size of the exposure, please add them manually") }

  if (any(is.na(dat$samplesize.outcome))) { warning("There are NA values in the total sample size of the outcome, please add them manually") }

  return(dat)

}
