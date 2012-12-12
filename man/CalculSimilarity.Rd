\name{CalculSimilarity}
\alias{CalculSimilarity}
\title{
Similarity Score between a Sequence and a PSSM model
}
\description{ The Similarity between a Sequence and a PSSM model is calculated using the score of the Sequence, and the minimum and maximum scores that can be obtained with the model. As the score approaches to 1, the sequence is more likely to belong to the modeled motif 
}
\usage{
CalculSimilarity(current, minim, maxim)
}
\arguments{
  \item{current}{
Score of the studied sequence
}
  \item{minim}{
minimum score obtained using the PSSM model
}
  \item{maxim}{
maximum score obtained using the PSSM model
}
}

\value{
similarity: between 0 and 1, similarity of the given sequence to the model
}
\references{
A.E. Kel , E. Gossling , I. Reuter , E. Cheremushkin , O.V. Kel-Margoulis , and E. Wingender
MATCHTM: a tool for searching transcription factor binding sites in DNA sequences 
Nucl. Acids Res. 31: 3576-3579.
}
\author{
Erola Pairo <epeiro@ibec.pcb.ub.es>
}
\examples{
require("MEET")
data(TFlogodds)
data(BackgroundOrganism)
data(Sequence)
current<-CalculScores(Sequence, TFlogodds)

maxim<-0
minim<-0

for(j in 1:dim(TFlogodds)[1]){
	  minim <- min(TFlogodds[j,])+minim
	  maxim <- max(TFlogodds[j,])+maxim
	  }
Similarity<-CalculSimilarity(current, minim, maxim)
  }


