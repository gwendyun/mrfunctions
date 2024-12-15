#' @title U1_Clean_data_exp
#' @description clean exp. exp does not require that many of snp
#'
#' @param df Data frame
#' @param continuous If continuous variable, will not ask for case sample size
#' @param save Default FALSE, if TRUE will save as local file, filename will be the name of id
#'
#' @return Returns cleaned data
#'
#' @examples df <- U1_Clean_data_exp(df,continuous=FALSE,save=FALSE)
#'
#'
#'
U1_Clean_data_exp <- function(df,continuous=FALSE,save=FALSE){
  require(Oneclick)

  if(nrow(df)<=20){
    warning('The number of SNPs is too small, it may not be suitable for exp analysis')
  }
  df<- Clean_data(df=df,continuous=continuous)
  class(df)<- c("df" , class(df) )

  if(save){
    arrow::write_parquet(df, paste0(df$id[1], ".parquet" ) , compression = "gzip")
  }

  return(df)
}


Clean_data = function(df,continuous=FALSE){

  rename_col <- function(df ,
                         patterns ,
                         format ){
    B <- colnames( df )

    for(i in 1:length(B)){
      col <- B[i]
      for( y in 1:length(patterns) ){
        if( tolower(patterns[y]) == tolower(B[i])){

          if(!B[i]==format){cat("Replace",B[i],"with", format,"\n")}

          B[i] <- format

        }
      }
    }

    colnames( df ) <-  B

    if(length( unique(colnames(df) ) ) != length(colnames( df))){
      stop(paste(format,"Column name is duplicated, please use `help(U1_Clean_data)` to view the documentation, only provide one!\n"), call. = FALSE)
    }

    return(df)
  }

  # Convert to data frame if it is a data frame
  if(is.data.frame(df)){ df = as.data.frame(df) }


  df<-rename_col(df,patterns=c("snp","rsid","rsids","snpid","rnpid","rs","variant_id"),format="SNP")
  df<-rename_col(df,patterns=c("chr","#chrom","chromosome"),format="chr")
  df<-rename_col(df,patterns=c("pos","position","base_pair_location"),format="pos")
  df<-rename_col(df,patterns=c("effect_allele","ea","alt", "alts","Tested_Allele","Alternate.Allele"),format="effect_allele")
  df<-rename_col(df,patterns=c("other_allele","oa","ref","Reference.Allele","NEA"),format="other_allele")
  df<-rename_col(df,patterns=c("beta","b","Effect","LogOR"),format="beta")
  df<-rename_col(df,patterns=c("se","sebeta","standard error","standard_error","StdErr","StdErrLogOR"),format="se")
  df<-rename_col(df,patterns=c("pval","p","p_value","pvalue"),format="pval")
  df<-rename_col(df,patterns=c("z","zscore"),format="z")
  df<-rename_col(df,patterns=c("eaf","FREQ", "af_alt", "FREQ1","effect_allele_frequency","Freq_Tested_Allele","Alternate.Allele.Frequency"),format="eaf")
  df<-rename_col(df,patterns=c("samplesize","n","sample_size","TotalSampleSize"),format="samplesize")
  df<-rename_col(df,patterns=c("ncase","n_cases", "ncases", "n_case"),format="ncase")
  df<-rename_col(df,patterns=c("ncontrol","n_controls","ncontrols", "n_control","Ntotal" ),format="ncontrol")
  df<-rename_col(df,patterns=c("or","odds ratio","odds_ratio"),format="or")
  df<-rename_col(df,patterns=c("gene","nearest_genes"),format="gene")
  df<-rename_col(df,patterns=c("info"),format="info")

  message('usually. effect_allele = Allele1; other_allele = Allele2. but please check manually. and use df <- dplyr::rename(df, other_allele = Allele2)')


  if(!"SNP" %in% colnames( df ) ){
    warning("SNP is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide an SNP column! \n")
  }else{
    ratio = sum(grepl("^rs\\d+$", df$SNP))/length(df$SNP)
    if( ratio < 0.5 ){ warning("A large number or all of the SNPs are not in the format rs123456, please use",cli::style_underline(cli::col_br_red("help(U1_add_SNP)")),"to view the documentation, match the SNP column!\n")  }
  }

  if(!"chr" %in% colnames( df ) | !"pos" %in% colnames( df ) ){
    cat("No chr or pos, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide a chr or pos column! Or use",cli::style_underline(cli::col_br_red("help(U1_add_chr_pos)")),"to view the documentation to match chr and pos, it does not affect the analysis if not available \n")
  }

  if( "pos" %in% colnames(df) ){
    df$pos<-as.numeric( df$pos ) }


  # beta part
  if( "or" %in% colnames( df ) & !"beta" %in% colnames( df ) ){
    cat("No beta provided, calculated based on or, formula is beta = log(or) \n")
    df$beta<-as.numeric(log(df$or))
  }else if("se" %in% colnames( df ) &  "pval" %in% colnames( df ) & !"beta" %in% colnames( df )  ){
    cat("No beta provided, calculated based on se and pval, formula is beta = se * sqrt(qchisq(pval,1,lower.tail=F) ) \n")
    df$beta<-as.numeric(  df$se * sqrt(qchisq(df$pval,1,lower.tail=F) )         )
  }else if(!"beta" %in% colnames( df )){
    warning("beta is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide a beta column! \n")
  }

  # se part
  if("beta" %in% colnames( df ) &  "pval" %in% colnames( df ) & !"se" %in% colnames( df )){
    cat("No se provided, calculated based on beta and pval, formula is se = beta / sqrt(qchisq(pval,1,lower.tail=F) ) \n")
    df$se<-as.numeric(  df$beta / sqrt(qchisq(df$pval,1,lower.tail=F) )         )
  }else if(!"se" %in% colnames( df )){
    warning("se is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide an se column! \n")
  }

  if(!"effect_allele" %in% colnames( df ) ){
    warning("effect_allele is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide an effect_allele column! \n")
  }

  if(!"other_allele" %in% colnames( df ) ){
    warning("other_allele is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide an other_allele column! \n")
  }


  if( "pval" %in% colnames(df) ){
    if (!is.numeric(df$pval)){
      cat("Convert pval to numeric\n")
      df$pval<-as.numeric(df$pval) }
  }else if(("beta" %in% colnames(df)) & ("se" %in% colnames(df))  & !("pval" %in% colnames(df))   ){
    cat("No pval provided, calculated based on beta and se, formula is pval<- 2*pnorm(abs(beta/se),lower.tail=FALSE) \n")
    df$beta<-as.numeric(df$beta)
    df$se<-as.numeric(df$se)
    df$pval<- 2*pnorm(abs(df$beta/df$se),lower.tail=FALSE)
    df$pval<-as.numeric(df$pval)
  }else if( "z" %in% colnames(df)  & !("pval" %in% colnames(df)) ){
    cat("No pval provided, cannot be calculated based on beta and se, calculated based on z, formula is pval<- 2*pnorm(abs(z),lower.tail=FALSE) \n")
    df$z<-as.numeric(df$z)
    df$pval<- 2*pnorm(abs(df$z),lower.tail=FALSE)
    df$pval<-as.numeric(df$pval)
  }else if(!("pval" %in% colnames(df))){
    warning("pval is very important, please use",cli::style_underline(cli::col_br_red("help(U1_Clean_data)")),"to view the documentation, provide a pval column! \n")
  }

  # if( !("z" %in% colnames(df)) & ("pval" %in% colnames(df)) ){
  #   cat("No z provided, calculated based on pval, formula is z <- sqrt(qchisq(pval,1,lower.tail=F)) \n")
  #   df$z<- sqrt(qchisq(pval,1,lower.tail=F))
  # }

  if(!"Phenotype" %in% colnames(df) ){
    df$Phenotype<- readline("Enter phenotype name, e.g., Body Mass Index, then press Enter:")
  }

  if(!"id" %in% colnames(df) ){
    df$id<- readline("Create a unique ID for your phenotype, e.g., Oneclick-0001, then press Enter:")
  }

  if(!continuous){
    if( (!"ncase" %in% colnames(df)) |  (!"ncontrol" %in% colnames(df))       ){
      YES_NO<- yesno::yesno2("Is it case-control data, e.g., disease with cases and controls? Select 1 for yes, 2 for no", yes = "Yes", no = "No")
      if( (!"ncase" %in% colnames(df)) & YES_NO    ){df$ncase<- as.numeric( readline("Enter number of cases, e.g., 4560, then press Enter:")  )}
      if( (!"ncontrol" %in% colnames(df)) & YES_NO  ){ df$ncontrol<- as.numeric( readline("Enter number of controls, e.g., 78200, then press Enter:")  )}
    }
  }
  if(("ncase" %in% colnames(df)) & ("ncontrol" %in% colnames(df)) & !"samplesize" %in% colnames(df) ){
    cat("No samplesize provided, calculated based on ncase and ncontrol, formula is samplesize<- ncase + ncontrol \n")
    df$samplesize<- as.numeric(df$ncase + df$ncontrol  )
  }else if( !"samplesize" %in% colnames(df) ){
    df$samplesize<- as.numeric( readline("Enter total sample size, e.g., 478000, then press Enter:")  )
  }

  return(df)

}



prepare_Munge<-function(){

  options(timeout = 1200)

  if (!require("BiocManager")) install.packages("BiocManager")

  # https://www.bioconductor.org/about/mirrors/
  options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")

  if (!require("MungeSumstats")){
    BiocManager::install("MungeSumstats")
  }


  # 38 reference group
  if (!require("SNPlocs.Hsapiens.dbSNP155.GRCh38")) {
    message( "The first installation of the SNPlocs.Hsapiens.dbSNP155.GRCh38 package takes a long time" )
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
  }
  if (!require("BSgenome.Hsapiens.NCBI.GRCh38")){
    message( "The first installation of the BSgenome.Hsapiens.NCBI.GRCh38 package takes a long time" )
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38") }

  # 37 reference group
  if (!require("SNPlocs.Hsapiens.dbSNP155.GRCh37")){
    message( "The first installation of the SNPlocs.Hsapiens.dbSNP155.GRCh37 package takes a long time" )
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
  }
  if (!require("BSgenome.Hsapiens.1000genomes.hs37d5")){
    message( "The first installation of the BSgenome.Hsapiens.1000genomes.hs37d5 package takes a long time" )
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
  }

  options(timeout = 60)


  library(MungeSumstats)
}
