\name{MImemory}
\alias{MImemory}
\title{Library of PredictDivergence values
}
\description{ This function calculates divergence values from all combinations possibles of a couple of nucleotides given a matrix of TFBS sequences. The divergence values is calculated only by TFBS positions correlated. }
\usage{
MImemory(iicc, training.set)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package}
  \item{training.set}{ A set of nucleotide sequences}
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
MIread}

\examples{
require("MEET")
data(iicc)
iicc <- detector_2nOrdre_init(training.set=iicc$Transcriptionfactor, val.set=val.set, iicc)
MImemory(iicc,training.set=iicc$Transcriptionfactor)
}

