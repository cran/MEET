\name{ModelDivergence}
\alias{ModelDivergence}
\title{To create Model Divergence}
\description{ From the best parameter, this function, ModelDivergence, calculates the divergence matrix of the set of aligned nucleotide sequences.}
\usage{
ModelDivergence(iicc)
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
iicc$parametersIdeal<-0.5
ModelDivergence(iicc)
}

