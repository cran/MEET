\name{PredictDivergence}
\alias{PredictDivergence}
\title{
A set of functions for detection of Transcription Factor Binding Sites by means of Divergence}
\description{
PredictDivergence contains a set of functions for detection of Transcription Factor Binding sites through parametric divergence.  This algorithm evaluates the variation on the total Renyi entropy of a set of sequences assuming correlation between positions in the bindng sequence.  
}
\usage{
PredictDivergence(iicc)
}
\arguments{
  \item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
}
\details{ This function integrates the Mutual Information (Renyi order equal 1) and parametric divergence (Renyi order different 1)
}
\references{ J. Maynou, M. Vallverdu, F. Claria, J.J. Gallardo-Chacon, P. Caminal and A. Perera," Transcription Factor Binding Site Detection through Position Cross-Mutual Information variability analysis". 31st Annual International Conference of the IEEE Engineering in Medicine and Biology Society.
}
\author{
Joan Maynou <joan.maynouatupc.edu>
}
\seealso{
PredictEntropy, PredictMEME, PredictMDscan, PredictPCA and PredictMATCH
}
\examples{
require("MEET")
data(HomoDivergence)
iicc[["AP1"]]$threshold<-0.01
PredictDivergence(iicc[["AP1"]])
}

