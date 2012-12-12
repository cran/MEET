\name{divergence.Shannon}
\alias{divergence.Shannon}
\title{
Divergencia.Shannon: Mutual Information
}
\description{
This function calculates Mutual Information (Renyi Order equal 1) by means of Kullback-Leibler divergence}
\usage{
divergence.Shannon(training.set, H, HXY,correction)
}
\arguments{
  \item{training.set}{A set of aligned nucleotide sequences}
  \item{H}{Entropy}
  \item{HXY}{Joint Entropy}
  \item{correction}{Correction of the Finite Sample Size Effect}
}
\details{Renyi Order has to be equal 1.}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
divergence.Renyi, PredictDivergence, kfold.divergence
}
\examples{
require("MEET")
data(TranscriptionFactor)
data(BackgroundOrganism)
data(iicc)
q<-1
training.set<-TranscriptionFactor
correction<-correction.entropy(q,p=nrow(training.set),long=1,iicc)
HXmax<-entropy.Shannon(as.matrix(Prob))
pmX<-probability(training.set,Prob)
Probtrans<-probability.couple(Prob)
H<-entropy.Shannon(pmX)
pmXY<-joint.probability(training.set, Prob, Probtrans)
HXY<-entropy.joint(pmXY,q,iicc)
divergence.Shannon(training.set,H,HXY,correction)}

