\name{entropy.joint}
\alias{entropy.joint}
\title{ To calculate joint entropy}
\description{
This function calculates the joint entropy between two variables X and Y with N possible state. 
}
\usage{
entropy.joint(pmXY,q,iicc)
}

\arguments{
  \item{pmXY}{ Joint probability}
  \item{q}{Renyi Order}
  \item{iicc}{ A list of options}
}
\details{If q is equal 1, Renyi joint entropy converges to Shannon joint entropy}
\author{Joan Maynou <joan.maynouatupc.edu>}
\seealso{entropy}
\examples{
require("MEET")
data(BackgroundOrganism)
data(iicc)
training.set<-iicc$Transcriptionfactor
Probtrans<-probability.couple(Prob)
pmXY<-joint.probability(training.set, Prob, Probtrans)
test<-entropy.joint(pmXY,q=iicc$q,iicc)}

