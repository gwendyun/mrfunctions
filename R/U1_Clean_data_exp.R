#' @title U1_Clean_data_exp
#' @description Clean exposure data. Exposure does not require many SNPs. https://gitee.com/OneClickAnalyses/Oneclick/blob/master/R/U1_Clean_data.R
#'
#' @param df Dataframe
#' @param continuous If TRUE, will not ask for case sample size.
#' @param save Default is FALSE. If TRUE, will save the dataframe as a local file with the name of the id.
#'
#' @return Returns the cleaned data.
#'
#' @examples
#' df <- U1_Clean_data_exp(df)
#'
U1_Clean_data_exp <- function(df, continuous = FALSE, save = FALSE) {
  if (nrow(df) <= 30) {
    warning('Too few SNPs, may not be suitable for exp analysis')
  }
  df <- Clean_data(df = df, continuous = continuous)
  class(df) <- c("df", class(df))

  if (save) {
    arrow::write_parquet(df, paste0(df$id[1], ".parquet"), compression = "gzip")
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

  # dataframe
  if(is.data.frame(df)){ df = as.data.frame(df) }


  df <- rename_col(df, patterns = c("snp", "rsid", "rsids", "snpid", "rnpid", "rs", "variant_id"), format = "SNP")
  df <- rename_col(df, patterns = c("chr", "#chrom", "chromosome"), format = "chr")
  df <- rename_col(df, patterns = c("pos", "position", "base_pair_location"), format = "pos")
  df <- rename_col(df, patterns = c("effect_allele", "ea", "alt", "alts", "Tested_Allele", "Alternate.Allele"), format = "effect_allele")
  df <- rename_col(df, patterns = c("other_allele", "oa", "ref", "Reference.Allele", "NEA"), format = "other_allele")
  df <- rename_col(df, patterns = c("beta", "b", "Effect", "LogOR"), format = "beta")
  df <- rename_col(df, patterns = c("se", "sebeta", "standard error", "standard_error", "StdErr", "StdErrLogOR"), format = "se")
  df <- rename_col(df, patterns = c("pval", "p", "p_value", "pvalue"), format = "pval")
  df <- rename_col(df, patterns = c("z", "zscore"), format = "z")
  df <- rename_col(df, patterns = c("eaf", "FREQ", "af_alt", "FREQ1", "effect_allele_frequency", "Freq_Tested_Allele", "Alternate.Allele.Frequency"), format = "eaf")
  df <- rename_col(df, patterns = c("samplesize", "n", "sample_size", "TotalSampleSize"), format = "samplesize")
  df <- rename_col(df, patterns = c("ncase", "n_cases", "ncases", "n_case"), format = "ncase")
  df <- rename_col(df, patterns = c("ncontrol", "n_controls", "ncontrols", "n_control", "Ntotal"), format = "ncontrol")
  df <- rename_col(df, patterns = c("or", "odds ratio", "odds_ratio"), format = "or")
  df <- rename_col(df, patterns = c("gene", "nearest_genes"), format = "gene")
  df <- rename_col(df, patterns = c("info"), format = "info")

  message('usually. effect_allele = Allele1; other_allele = Allele2. but please check manually. and use df <- dplyr::rename(df, other_allele = Allele2)')


  if (!"SNP" %in% colnames(df)) {
    warning("SNP is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide an SNP column! \n")
  } else {
    ratio <- sum(grepl("^rs\\d+$", df$SNP)) / length(df$SNP)
    if (ratio < 0.5) {
      warning("A large number or all of the SNPs are not in the format rs123456, please use", cli::style_underline(cli::col_br_red("help(U1_add_SNP)")), "to check the documentation and match the SNP column! \n")
    }
  }


  if (!"chr" %in% colnames(df) | !"pos" %in% colnames(df)) {
    cat("No chr or pos, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide a chr or pos column! Or use", cli::style_underline(cli::col_br_red("help(U1_add_chr_pos)")), "to check the documentation and match chr and pos. If not available, it will not affect the analysis.\n")
  }

  if ("pos" %in% colnames(df)) {
    df$pos <- as.numeric(df$pos)
  }

  # Beta part
  if ("or" %in% colnames(df) & !"beta" %in% colnames(df)) {
    cat("No beta provided, calculating from or, formula: beta = log(or)\n")
    df$beta <- as.numeric(log(df$or))
  } else if ("se" %in% colnames(df) & "pval" %in% colnames(df) & !"beta" %in% colnames(df)) {
    cat("No beta provided, calculating from se and pval, formula: beta = se * sqrt(qchisq(pval, 1, lower.tail = F))\n")
    df$beta <- as.numeric(df$se * sqrt(qchisq(df$pval, 1, lower.tail = FALSE)))
  } else if (!"beta" %in% colnames(df)) {
    warning("Beta is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide a beta column! \n")
  }

  # SE part
  if ("beta" %in% colnames(df) & "pval" %in% colnames(df) & !"se" %in% colnames(df)) {
    cat("No se provided, calculating from beta and pval, formula: se = beta / sqrt(qchisq(pval, 1, lower.tail = F))\n")
    df$se <- as.numeric(df$beta / sqrt(qchisq(df$pval, 1, lower.tail = FALSE)))
  } else if (!"se" %in% colnames(df)) {
    warning("SE is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide an se column! \n")
  }

  if (!"effect_allele" %in% colnames(df)) {
    warning("Effect allele is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide an effect_allele column! \n")
  }

  if (!"other_allele" %in% colnames(df)) {
    warning("Other allele is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide an other_allele column! \n")
  }

  if ("pval" %in% colnames(df)) {
    if (!is.numeric(df$pval)) {
      cat("Converting pval to numeric\n")
      df$pval <- as.numeric(df$pval)
    }
  } else if ("beta" %in% colnames(df) & "se" %in% colnames(df) & !"pval" %in% colnames(df)) {
    cat("No pval provided, calculating from beta and se, formula: pval <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)\n")
    df$beta <- as.numeric(df$beta)
    df$se <- as.numeric(df$se)
    df$pval <- 2 * pnorm(abs(df$beta / df$se), lower.tail = FALSE)
    df$pval <- as.numeric(df$pval)
  } else if ("z" %in% colnames(df) & !"pval" %in% colnames(df)) {
    cat("No pval provided, calculating from z, formula: pval <- 2 * pnorm(abs(z), lower.tail = FALSE)\n")
    df$z <- as.numeric(df$z)
    df$pval <- 2 * pnorm(abs(df$z), lower.tail = FALSE)
    df$pval <- as.numeric(df$pval)
  } else if (!"pval" %in% colnames(df)) {
    warning("Pval is very important, please use", cli::style_underline(cli::col_br_red("help(U1_Clean_data)")), "to check the documentation and provide a pval column! \n")
  }

  if (!"Phenotype" %in% colnames(df)) {
    df$Phenotype <- readline("Enter the phenotype name, e.g., Body Mass Index, and press Enter:")
  }

  if (!"id" %in% colnames(df)) {
    df$id <- readline("Create a unique ID for your phenotype, e.g., Oneclick-0001, and press Enter:")
  }

  if (!continuous) {
    if (!"ncase" %in% colnames(df) | !"ncontrol" %in% colnames(df)) {
      YES_NO <- yesno::yesno2("Is it case-control data, e.g., disease with cases and controls? Yes select 1, No select 2", yes = "Yes", no = "No")
      if (!"ncase" %in% colnames(df) & YES_NO) {
        df$ncase <- as.numeric(readline("Enter the number of cases, e.g., 4560, and press Enter:"))
      }
      if (!"ncontrol" %in% colnames(df) & YES_NO) {
        df$ncontrol <- as.numeric(readline("Enter the number of controls, e.g., 78200, and press Enter:"))
      }
    }
  }
  if ("ncase" %in% colnames(df) & "ncontrol" %in% colnames(df) & !"samplesize" %in% colnames(df)) {
    cat("No samplesize provided, calculating from ncase and ncontrol, formula: samplesize <- ncase + ncontrol\n")
    df$samplesize <- as.numeric(df$ncase + df$ncontrol)
  } else if (!"samplesize" %in% colnames(df)) {
    df$samplesize <- as.numeric(readline("Enter the total sample size, e.g., 478000, and press Enter:"))
  }

  return(df)
}

prepare_Munge <- function() {

  options(timeout = 1200)

  if (!require("BiocManager")) install.packages("BiocManager")

  # https://www.bioconductor.org/about/mirrors/
  options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

  if (!require("MungeSumstats")) {
    BiocManager::install("MungeSumstats")
  }

  # GRCh38 reference group
  if (!require("SNPlocs.Hsapiens.dbSNP155.GRCh38")) {
    message("Installing SNPlocs.Hsapiens.dbSNP155.GRCh38 package for the first time may take a while")
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
  }
  if (!require("BSgenome.Hsapiens.NCBI.GRCh38")) {
    message("Installing BSgenome.Hsapiens.NCBI.GRCh38 package for the first time may take a while")
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
  }

  # GRCh37 reference group
  if (!require("SNPlocs.Hsapiens.dbSNP155.GRCh37")) {
    message("Installing SNPlocs.Hsapiens.dbSNP155.GRCh37 package for the first time may take a while")
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
  }
  if (!require("BSgenome.Hsapiens.1000genomes.hs37d5")) {
    message("Installing BSgenome.Hsapiens.1000genomes.hs37d5 package for the first time may take a while")
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
  }

  options(timeout = 60)

  library(MungeSumstats)
}
