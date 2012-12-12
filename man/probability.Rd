\name{probability}
\alias{probability}
\title{
Probability
}
\description{
This function calculates the probability of random varible with N states considering missing values.
}
\usage{
probability(wind, Prob)
}
\arguments{
  \item{wind}{A set of nucleotide aligned sequences}
  \item{Prob}{Probability of each base within genome}
  }

\author{
Joan Maynou <joan.maynouatupc.edu>}

\seealso{
joint.probability, probIdem, probtransicio
}
\examples{
require("MEET")
data(TranscriptionFactor)
data(BackgroundOrganism)
wind<-TranscriptionFactor
probability(wind,Prob)
}

