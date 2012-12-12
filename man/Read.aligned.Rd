\name{Read.aligned}
\alias{Read.aligned}
\title{
Read nucleotide sequences}
\description{
This function reads nucleotide sequences in format .fasta.
}
\usage{
Read.aligned(file)
}
\arguments{
  \item{file}{ A set of nucleotide sequences}}
\details{
The file has to be in format .fasta}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
ReadSequence, ReadTF
}
\examples{
require("MEET")
pathMEET <- system.file("sequences", package = "MEET")
Read.aligned(file=paste(pathMEET, "AP1.fa", sep = "/"))
}
