\name{CalculScores}
\alias{CalculScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calcul Score of a Sequence, using a loggods matrix}
\description{ Giving a logodds matrix, the Score of a sequence is calculated as the sum of the logodds of each nucleotide at each position of the sequence (equivalent to calculating the probability of each nucleotide at each position)

}
\usage{
CalculScores(sequencia, logodds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequencia}{
Sequence to analyze
}
  \item{logodds}{
logodds matrix for a given motif. 
}
}

\value{
score: score of the studied sequence
}
\references{
 Gary D. Stormo. DNA binding sites: representation and discovery. Bioinformatics (2000) 16(1): 16-23 doi:10.1093/bioinformatics/16.1.16
}
\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}


\seealso{
CalculPSSM, CalculInformation
}
\examples{
require("MEET")
#First example, calculating logodds of a Sequence 
data(TFlogodds)
data(Sequence)
CalculScores(sequencia=Sequence, logodds=TFlogodds)
#given a Transfac matrix, calculate first loggods and then scores
data(TranscriptionFactor)
data(BackgroundOrganism)
data(Sequence)
Factortrans<-TranscriptionFactor
suma<-apply(Factortrans,2,function(y){sum(y=="-")})
Factortrans<-Factortrans[, suma==0]
logodds<-CalculInformation(matriu=Factortrans, Prob=Prob)
logodds <- logodds[,2:dim(logodds)[2]]
CalculScores(sequencia=Sequence,logodds=logodds)
  
}

