\name{detector_2nOrdre}
\alias{detector_2nOrdre}
\title{Detection of Transcription Factor Binding Sites Through Parametric PredictDivergence}
\description{
This detection algorithm is based on Information Theory.  Specifically, this method uses a parametric divergence.  This algorithm evaluates the variation on the total Renyi entropy of a set of sequences assuming correlation between positions in the bindng sequence.  When a candidate sequence is assumed to be a a true binding site belonging to the set.The measurement of the variation of the total redundancy when the candidate sequence is added to the set has been computed by using the difference between the redundancy profile.
}
\usage{
detector_2nOrdre(training.set, val.set, iicc)
}

\arguments{
  \item{training.set}{A set of aligned nucleotide sequences
}
  \item{val.set}{A candidate sequence}
  \item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend) and the initial conditions for the divergence method (PredictDivergence, correction 1rOrdre, Exterior product entropic profile (Mperfil), maximun entropy, entropy and Renyi order)}
}
\details{Options parameter has to contain the initial conditions for the divergence method: divergence matrix (D), maximum entropy (HXmax),  correction entropy and redundancy from the Finite Sample Size Effect, exterior product entropic profile (Mperfil), entropy (HX) and the Renyi order.}

\author{
Joan Maynou <joan.maynouatupc.edu>
}
\seealso{
detector_1rOrdre_diff
}
\examples{
require("MEET")
data(BackgroundOrganism)
data(iicc)
training.set<-iicc$Transcriptionfactor
val.set<-sample(c('A','T','C','G'),ncol(training.set), replace=TRUE,Prob)
iicc<-detector_2nOrdre_init(training.set, val.set, iicc)
out<-detector_2nOrdre(training.set, val.set, iicc)
}

