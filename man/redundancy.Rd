\name{redundancy}
\alias{redundancy}

\title{To calculate redundancy
}
\description{
This function calculates the redundancy of random variable X with N possibles states
}
\usage{
redundancy(HX, HXmax, Herror, finite, correction = TRUE)
}
\arguments{
  \item{HX}{A entropy vector}
  \item{HXmax}{Maximum entropy}
  \item{Herror}{Error entropy}
  \item{finite}{Finite Sample Size}
  \item{correction}{If TRUE, considers correction of the finite Sample size}}
\details{
This function depends the correction.redundancy}

\author{
Joan Maynou <joan.maynouatupc.edu>
}
\seealso{
correction.redundancy, correction.entropy, entropy.Shannon,entropy.Renyi}
\examples{
require("MEET")
data(BackgroundOrganism)
data(iicc)
q<-iicc$q
Factortrans<-iicc$Transcriptionfactor
HXmax<-iicc$HXmax
correction<-correction.entropy(q,p=nrow(Factortrans),long=1,iicc)
Herror<-slot(correction,"Herror")
pmX<-probability(Factortrans,Prob)
HX<-switch(iicc$classentropy,"Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))
test<-redundancy(HX,HXmax,Herror,finite=nrow(Factortrans),correction=TRUE)}

