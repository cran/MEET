\name{entropy.Renyi}
\alias{entropy.Renyi}
\title{
Renyi Entropy 
}
\description{
This function calculates the Renyi's entropy of the variable X with N possible states
}
\usage{
entropy.Renyi(wind, q)
}

\arguments{
  \item{wind}{Probability matrix (4xm).}
  \item{q}{Renyi Order}
}
\details{ wind parameter is calculated by means of probability function}

\author{
Joan Maynou <joan.maynouatupc.edu>}
\seealso{
entropy.max, entropy.conjunta, probability, entropy.corrected }
\examples{
 require("MEET")
 data(TranscriptionFactor)
 data(RenyiOrder)
 data(BackgroundOrganism)
 training.set<-TranscriptionFactor
 wind<-probability(training.set,Prob)
 test<-entropy.Renyi(wind,q)}

