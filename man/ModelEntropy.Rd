\name{ModelEntropy}
\alias{ModelEntropy}
\title{To create Model Entropy}
\description{ From the best parameter, this function, ModelEntropy, calculates the redundancy profile of the set of aligned nucleotide sequences.}
\usage{
ModelEntropy(iicc)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package}
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}
\examples{
require("MEET")
data(iicc)
iicc$parametersIdeal<-1
ModelEntropy(iicc)
}

