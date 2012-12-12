\name{CalculPSSM}
\alias{CalculPSSM}
\title{
 Position Specific Scoring Matrices from a set of aligned sequences
}
\description{ Calculate a PSSM of a set of DNA aligned sequences, taking into account background probabilities. The output is the logodds matrix.
}
\usage{
CalculPSSM(matriu, Prob)
}
\arguments{
  \item{matriu}{ Set of DNA aligned sequences}
  \item{Prob}{ Probability distribution of the nucleotides in the background model}
}
\value{
logodds: logodds matrix calculated as a log2 of the frequency matrix
}
\references{
  Gary D. Stormo. DNA binding sites: representation and discovery. Bioinformatics (2000) 16(1): 16-23 doi:10.1093/bioinformatics/16.1.16
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}
\seealso{
CalculInformation
}
\examples{
data(TranscriptionFactor)
data(BackgroundOrganism)
CalculPSSM(matriu=TranscriptionFactor, Prob=Prob)
}
