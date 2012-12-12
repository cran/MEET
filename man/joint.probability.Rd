\name{joint.probability}
\alias{joint.probability}
\title{Joint Probability}
\description{This function calculates the joint probability of each base-couple  among two positions of the training.set}
\usage{joint.probability(training.set,Prob,Probtrans)}

\arguments{
  \item{training.set}{ A set of aligned nucleotide sequence}
  \item{Prob}{ Background probability of each base within genome}
  \item{Probtrans}{Background probability of correlation among bases within genome}
}
\author{Joan Maynou <joan.maynouatupc.edu>}
\seealso{ probability}
\examples{
require("MEET")
data(TranscriptionFactor)
data(BackgroundOrganism)
training.set<-TranscriptionFactor
Probtrans<-probability.couple(Prob)
pmXY<-joint.probability(training.set, Prob,Probtrans)
}





