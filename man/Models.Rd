\name{Models}
\alias{Models}
\title{To create Detection Model}
\description{ From the best parameter, this function,ConstructModel, calculates the detection model of the set of aligned nucleotide sequences. There are different kinds of models: Entropy, Divergence, PCA, MEME, MDscan, MATCH and transMEME}
\usage{
Models(iicc)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package}
}

\author{
Joan Maynou <joan.maynouatupc.edu> and Erola Pairo  <epeiroatibec.pcb.ub.es>
}
\examples{
require("MEET")
data(iicc)
iicc$parametersIdeal<-1
Models(iicc)
}

