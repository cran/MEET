\name{align.MEME}
\alias{align.MEME}
\title{
Multiple sequence alignment by means of MEME.}
\description{ DNA sequences are aligned by means of MEME Version 4.4.0. (Multiple Expectation-Maximization for Motif Elicitation)}
\usage{align.MEME(filein, fileout = "Sq.fa", iicc)}
\arguments{
  \item{filein}{ A set of nucleotide sequences in FASTA format.}
  \item{fileout}{ Output file in FASTA format}
  \item{iicc}{ A list of argument input.}
}
\details{ This function needs aaMI-package.
	  This funtions works with meme<=4.3.0}
\value{Output is a file in FASTA format with aligned nucleotide sequences.}
\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>
}
\seealso{align.clustalw, align.MUSCLE}


