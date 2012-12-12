\name{CreateConsensus}
\alias{CreateConsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Consensus Sequence for a DNA motif
}
\description{
Using seqinr R package, constructs a consensus DNA sequences from a set of aligned Sequences
}
\usage{
CreateConsensus(alignedSequences,iicc,filein)
}

\arguments{
	\item{alignedSequences}{Aligned DNA sequences}
	\item{iicc}{Initial options}
	\item{filein}{Fasta file with a set of aligned sequences with MEME}
}
\details{
package "seqinr" is required
}
\value{
SeqCons: Consensus sequence (if gaps a "Na" is returned)
}
\references{
https://r-forge.r-project.org/projects/seqinr/
}
\author{
Erola Pairo <epeiroatibec.ub.pcb.es>
}


\seealso{
align.clustalw, align.MEME, read.aligned, align.muscle
}
\examples{
require("MEET")
data(TranscriptionFactor)
data(iicc)
pathMEET <- system.file("sequences", package = "MEET")
SeqCons<-CreateConsensus(TranscriptionFactor,iicc,filein=NULL)
}

