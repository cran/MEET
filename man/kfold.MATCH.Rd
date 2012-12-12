\name{kfold.MATCH}
\alias{kfold.MATCH}
\title{MATCH validation process}
\description{
For a vector of Core cut values this function calculates returns a matrix of Similarities and labels indicating if a sequence position is a binding site or not. To calculate the Similartiy a leave-one-out cross validation model is used. With the utput of this function a ROC curve can be calculated for each Core Similarity and the results can be compared.
}
\usage{
kfold.MATCH(iicc, Seqin)
}
\arguments{
  \item{iicc}{List of options described in the MEET program}
  \item{Seqin}{DNA sequence }
}

\references{
A.E. Kel , E. Gossling , I. Reuter , E. Cheremushkin , O.V. Kel-Margoulis , and E. Wingender
MATCHTM: a tool for searching transcription factor binding sites in DNA sequences 
Nucl. Acids Res. 31: 3576-3579.}
\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}
\seealso{
Match, MEET
}
\examples{
require("MEET")
require("seqinr")
data(iicc)
data(TranscriptionFactor)
pathMEET <- system.file("sequences", package = "MEET")
iicc$method<-"MATCH"
iicc$vector<-c(0.5, 0.8)
kfold.MATCH(iicc, Seqin = paste(pathMEET, "AP1.fa", sep = "/"))

}
