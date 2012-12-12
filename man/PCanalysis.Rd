\name{PCanalysis}
\alias{PCanalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PC analysis on numerical DNA sequences
}
\description{
Performs a PC analysis of numerical DNA sequences, using the pcaMethods package and projects the studied DNA sequence to the subspace. The Q-residuals (euclidean distance between the sequences and the modeled subspace) of the DNA studied sequence are returned as output.
}
\usage{
PCanalysis(TFBS, nPCs, Sequences)
}

\arguments{
  \item{TFBS}{
DNA numerical matrix of known TFBS
}
  \item{nPCs}{
Number of principal components used to build the model
}
  \item{Sequences}{
numerical matrix of DNA sequences to study (converted as a numerical matrix)
}
}

\value{
residus: Euclidean distance from the studied sequences to the modeled motif
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}

\seealso{
PCA, kfold.PCA
}
