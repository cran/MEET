\name{ReadSequence}
\alias{ReadSequence}
\title{Convert a DNA sequence in a numerical DNA matrix}
\description{Convert a DNA symbolical sequence to a DNA numerical matrix of length m, using a conversion where each nucleotide is placed at the vertex of a numerical tetrahedron}
\usage{ReadSequence(Seq, m,background,convertDNA)}
\arguments{
  \item{Seq}{Vector containing a DNA symbolical sequence}
  \item{m}{length of the matrix}
  \item{background}{Probability background}
  \item{convertDNA}{Numerical representation of each nucleotide}
}
\value{matriuSeq:numerical matrix representing de DNA sequence. }
\author{Erola Pairo <epeiroatibec.pcb.ub.es>}
\seealso{convertDNA, numericalDNA, readTF}
