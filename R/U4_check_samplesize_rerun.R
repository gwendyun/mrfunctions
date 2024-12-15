#' @title check samplesize
#'
#' @description Check if the sample size is correct, if not, it needs to be modified. If the name is not satisfactory, it also needs to be modified. Do remember rerun again after correction.
#' @param dat dat after U4_harmonise_data
#'
#' @return dat
#' @export
#'
#' @examples
#'
#' help(U4_harmonise_data)
#'
#'
#'
U4_check_samplesize_rerun<- function(dat){

  dat1 <- plyr::ddply(dat,"id.exposure",
                      function(dat) {

                        if( length( unique(dat$samplesize.exposure) )>1 ){
                          message("Exposure",dat$exposure[1],"has multiple different sample sizes, only the largest sample size is selected for subsequent calculations")
                        }

                        if( ( "ncase.exposure" %in% colnames(dat)) & ( "ncontrol.exposure" %in% colnames(dat) ) ){

                          dat$ncase.exposure<-  as.numeric(dat$ncase.exposure)
                          dat$ncontrol.exposure<-  as.numeric(dat$ncontrol.exposure)

                          samplesize.exposure <- as.numeric(  max(dat$samplesize.exposure) )

                          if( all(!is.na(dat$ncase.exposure))  & all( !is.na(dat$ncontrol.exposure)) ){

                            dat$type.exposure<- 'Binary'
                            dat$ncase.exposure<-dat$ncase.exposure[which(dat$samplesize.exposure == samplesize.exposure)[1] ]
                            dat$ncontrol.exposure<-dat$ncontrol.exposure[which(dat$samplesize.exposure == samplesize.exposure)[1] ]
                            dat$samplesize.exposure <- as.numeric(  max(dat$samplesize.exposure) )

                          }else{
                            dat$type.exposure<-"Continuous"
                            dat$samplesize.exposure <- as.numeric(  max(dat$samplesize.exposure) )
                          }

                        }else{
                          dat$samplesize.exposure <- as.numeric(  max(dat$samplesize.exposure) )
                          dat$type.exposure<-"Continuous"
                        }

                        return(dat)

                      } )

  for (i in 1:nrow(dat1) ) {
    dat1$description_exposure[i] <- paste0( "Exposure is ", dat1$exposure[i], ". Total sample size is ",dat1$samplesize.exposure[i],". Inferred as ",ifelse( dat1$type.exposure[i]=="Continuous", "continuous variable"  , paste0("binary variable, case sample size is ",dat1$ncase.exposure[i],", control sample size is ",dat1$ncontrol.exposure[i] )  ) )
  }

  message( "Please check if the following exposure descriptions are incorrect. If so, please use help(U4_harmonise_data) to see how to modify the sample size and other information.\n" )
  print( unique(dat1$description_exposure) )

  dat2 <- plyr::ddply(dat1,c("id.outcome"),
                      function(dat) {

                        if( length( unique(dat$samplesize.outcome) )>1 ){
                          message("Outcome",dat$outcome[1],"has multiple different sample sizes, only the largest sample size is selected for subsequent calculations")
                        }

                        if( ( "ncase.outcome" %in% colnames(dat)) & ( "ncontrol.outcome" %in% colnames(dat) ) ){

                          dat$ncase.outcome<-  as.numeric(dat$ncase.outcome)
                          dat$ncontrol.outcome<-  as.numeric(dat$ncontrol.outcome)

                          samplesize.outcome <- as.numeric(  max(dat$samplesize.outcome) )

                          if( all(!is.na(dat$ncase.outcome))  & all( !is.na(dat$ncontrol.outcome)) ){

                            dat$type.outcome<- 'Binary'
                            dat$ncase.outcome<-dat$ncase.outcome[which(dat$samplesize.outcome == samplesize.outcome)[1] ]
                            dat$ncontrol.outcome<-dat$ncontrol.outcome[which(dat$samplesize.outcome == samplesize.outcome)[1] ]
                            dat$samplesize.outcome <- as.numeric(  max(dat$samplesize.outcome) )

                          }else{
                            dat$type.outcome<-"Continuous"
                            dat$samplesize.outcome <- as.numeric(  max(dat$samplesize.outcome) )
                          }

                        }else{
                          dat$samplesize.outcome <- as.numeric(  max(dat$samplesize.outcome) )
                          dat$type.outcome<-"Continuous"
                        }
                        return(dat)
                      })


  for (i in 1:nrow(dat2) ) {
    dat2$description_outcome[i] <- paste0( "Outcome is ", dat2$outcome[i], ". Total sample size is ",dat2$samplesize.outcome[i],". Inferred as ",ifelse( dat2$type.outcome[i]=="Continuous", "continuous variable", paste0("binary variable, case sample size is ",dat2$ncase.outcome[i],", control sample size is ",dat2$ncontrol.outcome[i] )  ) )
  }

  message( "\n\nPlease check if the following outcome descriptions are incorrect. If so, please use help(U4_harmonise_data) to see how to modify the outcome name, sample size and other information.\n" )
  print( unique(dat2$description_outcome) )

  return( dat2 )

}
