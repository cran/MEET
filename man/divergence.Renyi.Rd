\name{divergence.Renyi}
\alias{divergence.Renyi}
\title{Renyi divergence}
\description{
This function calculates parametric divergence (Renyi Order different 1)
}
\usage{
divergence.Renyi(training.set, pmX, pmXY, q, correction)
}
\arguments{
  \item{training.set}{ A set of aligned nucleotide sequences}
  \item{pmX}{ Relative frequency of a nucleotide at a motif position (position independency model) as an estimation of the probability of this fact.}
  \item{pmXY}{To extend to pmX to include of correlated positions (position dependency model)}
 
  \item{q}{Renyi Order}
  \item{correction}{Correction of the Finite Sample Size Effect}
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
divergence.Shannon, PredictDivergence
}
\examples{
require("MEET")
data(iicc)
data(BackgroundOrganism)
training.set<-iicc$Transcriptionfactor
q<-iicc$q<-0.5
correction<-correction.entropy(q,p=nrow(training.set),long=1,iicc)
HXmax<-iicc$HXmax
pmX<-probability(training.set,Prob)
Probtrans<-probability.couple(Prob)
H<-entropy.Renyi(pmX,q)
pmXY<-joint.probability(training.set, Prob, Probtrans)
HXY<-entropy.joint(pmXY,q,iicc)
divergence.Renyi(training.set,pmX,pmXY,q,correction)
}

