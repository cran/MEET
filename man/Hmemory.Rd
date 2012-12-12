\name{Hmemory}
\alias{Hmemory}
\title{Library of entropy values
}
\description{ This function calculates entropy  values from all combinations possibles of nucleotides given a matrix of TFBS sequences. }
\usage{
Hmemory(iicc, training.set)
}

\arguments{
\item{iicc}{A set of initial conditions for the MEET-package}
  \item{training.set}{ A set of nucleotide sequences}
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
Hread}

\examples{
require("MEET")
data(iicc)
factor<-iicc$Transcriptionfactor
Hmemory(iicc,factor)
}

