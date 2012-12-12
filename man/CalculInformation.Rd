\name{CalculInformation}
\alias{CalculInformation}
\title{ Information content in each position of a set of aligned DNA sequences}
\description{
Using as an input the set of aligned DNA sequences, information content in each position is calculated, taking into account background probability for each nucleotide. The first row of the returned sequence is the information in each postion, and the others are the loggods matrix.
}
\usage{CalculInformation(matriu, Prob)}
\arguments{
  \item{matriu}{Aligned DNA sequences}
  \item{Prob}{Background Probability for the studied individual}
}
\value{
matriu is a matrix with a first row that is the information content in each position and the next ones the logodds matrix using information content.
}

\author{Erola Pairo <epeiroatibec.pcb.ub.es>}
\seealso{kfold.MATCH, CalculScores}
\examples{
data(TranscriptionFactor)
data(BackgroundOrganism)
CalculInformation(matriu=TranscriptionFactor, Prob=Prob)}


