\name{entropy.corrected}
\alias{entropy.corrected}
\title{
Correction of the Finite Sample Size Effect
}
\description{
Each training matrix is formed by a finite number of samples. The probability estimation error using the nucleotide frequency causes a bias on the uncertainty measurement 
}
\usage{entropy.corrected(H, ErrorHX, HXmax)}

\arguments{
  \item{H}{Entropy}
  \item{ErrorHX}{Error Entropy}
  \item{HXmax}{Maximum Entropy}}
\details{
This function uses the results of correction.entropy and correctionaprox functions.
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
correctionaprox, correction.entropy
}
\examples{
require("MEET")
data(BackgroundOrganism)
data(iicc)
Factortrans<-iicc$Transcriptionfactor
correction<-correction.entropy(q=iicc$q,p=nrow(Factortrans),long=1,iicc)
ErrorHX<-slot(correction,"sderror")[nrow(Factortrans)]
HXmax<-iicc$HXmax
prob<-probability(Factortrans,Prob)
H<-entropy.Shannon(prob)
test<-entropy.corrected(H,ErrorHX,HXmax)
}

