\name{PredictPCA}
\alias{PredictPCA}
\title{ Q-residuals detection of TFBS, using a principal components model}
\description{Performs the detection of a modeled motif in a chormosomics sequence}
\usage{PredictPCA(iicc)}

\arguments{
  \item{iicc}{options described in the MEET program}
}
\details{
The specific options for this program are the model constructed for the motif, the background probability for the studied organism and the DNA query sequence
}
\value{
 output:list of the TFBS found in the studied DNA sequence: Sequence (BS found), pvalue (Pvalue of the given BS) and position (Position within the studied DNA sequence).
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}

\seealso{
ModelPCA, biol.PCA, MEET, detection}
\examples{
require("MEET")
data(iicc)
data(TranscriptionFactor)
iicc$method<-"PCA"
iicc$pvalor<-0.1
iicc$parameterIdeal<-2
iicc$model<-ModelPCA(iicc)
DetectedSequences<-PredictPCA(iicc)
}