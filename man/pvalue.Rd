\name{pvalue}
\alias{pvalue}
\title{P value}
\description{
This function calculates p-value under numerical distribution}
\usage{
pvalue(a, dist)
}
\arguments{
  \item{a}{Numerical value}
  \item{dist}{Numerical distribution}
}
\details{
This function needs the next packages: KernSmooth and Hmisc}

\author{
Helena Brunel <helena.brunelatupc.edu>}

\examples{
require("MEET")
dist<-rnorm(10000,0,1)
pvalue(a=3.298,dist)}

