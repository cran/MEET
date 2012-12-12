\name{diffInstructions}
\alias{diffInstructions}
\title{The measurement of the variation of the total redundancy
}
\description{
This function measures the variation of the total redundancy.}
\usage{
diffInstructions(training.set, HX, HXmax, Herror, Redundancia_corregida)
}
\arguments{
  \item{training.set}{A set of aligned nucleotide sequences}
  \item{HX}{Entropy}
  \item{HXmax}{Maximum entropy}
  \item{Herror}{Entropy error.  Correction of the Finite Sample Size Effect}
  \item{Redundancia_corregida}{ Redundancy correction of the Finite Sample Size Effect}
}
\details{This function depends on detector_1erOrdre_diff}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\examples{
require("MEET")
data(iicc)
data(BackgroundOrganism)
q<-iicc$q
training.set<-iicc$Transcriptionfactor
correction<-correction.entropy(q,p=nrow(training.set),long=1,iicc)
Herror<-slot(correction,"Herror")
HXmax<-iicc$HXmax
prob<-probability(training.set,Prob)
HX<-entropy.Shannon(prob)
Redundancia_corregida<-CalculRedundancy(training.set,q,iicc)
test<-diffInstructions (training.set,HX,HXmax,Herror,Redundancia_corregida)}

