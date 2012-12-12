\name{align.clustalw}
\alias{align.clustalw}
\title{Multiple sequence alignment by means of ClustalW}
\description{ DNA sequences are aligned by means of ClustalW (Multiple Sequence Alignment).}
\usage{
align.clustalw(filein, fileout = "Sq.fa", call)
}
\arguments{
  \item{filein}{ A set of nucleotide sequences in FASTA format.}
  \item{fileout}{ Output file in FASTA format.}
  \item{call}{string of characters needed to call Clustalw from the working directory}
}
\details{This function needs aaMI-package}
\value{ Fileout is a file in FASTA format with aligned nucleotide sequences.}
\author{Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>.}
\seealso{align.MEME, align.MUSCLE}


