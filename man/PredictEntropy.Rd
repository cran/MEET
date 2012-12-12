\name{PredictEntropy}
\alias{PredictEntropy}
\title{
PredictEntropy: Detection of Transcription Factor Binding Sites by means of Renyi entropy
}
\description{Detection of transcription factor binding sites through parametric uncertainty measurement (Renyi entropy). This detection algorithm evaluates the variation on the total Renyi entropy of a set of sequences when a candidate sequence is assumed to be a true binding site belonging to the set.}
\usage{
PredictEntropy(iicc)
}

\arguments{
  \item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
}
\details{
This function contains a set of function for the detection of transcription factor binding sites:correction.entropy.R,
  correction.redundancy.R, entropy.R entropy.max.R, entropy.corrected.R, probability.R, CalculRedundancy.R, diff.instructions.R, redundancy.R, ROC.R, detector_1rOrdre_diff.R, pvalue.R}

\author{
Joan Maynou <joan.maynouatupc.edu>}

\seealso{PredictDivergence, PredictMEME, PredictMDscan, PredictPCA, PredictMATCH}
\examples{
require("MEET")
data(iicc)
PredictEntropy(iicc)}

