% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/U1_Clean_data_exp.R
\name{U1_Clean_data_exp}
\alias{U1_Clean_data_exp}
\title{U1_Clean_data_exp}
\usage{
U1_Clean_data_exp(df, continuous = FALSE, save = FALSE)
}
\arguments{
\item{df}{Data frame}

\item{continuous}{If continuous variable, will not ask for case sample size}

\item{save}{Default FALSE, if TRUE will save as local file, filename will be the name of id}
}
\value{
Returns cleaned data
}
\description{
clean exp. exp does not require that many of snp
}
\examples{
df <- U1_Clean_data_exp(df,continuous=FALSE,save=FALSE)



}
