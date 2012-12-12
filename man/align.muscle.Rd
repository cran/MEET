\name{align.muscle}
\alias{align.muscle}
\title{
Multiple sequence alignment by means of Muscle (MUltiple Sequence Comparison by Log-Expectation)
}
\description{
 DNA sequences are aligned by means of Muscle Version 3.8. (Multiple Sequences Alignment)
}
\usage{
align.muscle(filein , fileout = "Sq.fa", gapopen = gapopen, maxiters = maxiters, gapextend = gapextend, call)
}
\arguments{
  \item{filein}{A set of nucleotide sequences in FASTA format.}
  \item{fileout}{Output file in FASTA format.}
  \item{gapopen}{Gap open score.}
  \item{maxiters}{Maximum number of iterations.}
  \item{gapextend}{Gap extend score.}
  \item{call}{string of characters needed to call Muscle from the working directory}

}
\details{
Gapopen and gapextend must be negative. Output is a file in FASTA format with aligned nucleotide sequences. This function needs aaMI-package.
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
align.clustalw, align.MEME
}


