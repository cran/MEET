\name{correctionaprox}
\alias{correctionaprox}
\title{Correction Entropy Approximate from the Finite Sample Size Effect.}
\description{
Each training matrix is formed by a finite number of samples. The probability estimation error using the nucleotide frequency causes a bias on the uncertainty measurement The approximation method is used for sample size is more than or equal  to 50.}

\usage{
correctionaprox(x, matriu, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Maximum entropy}
  \item{matriu}{A set of aligned nucleotide sequences.}
  \item{s}{Number of symbols}
}
\details{Sample size has to be > 50.}

\author{
Joan Maynou <joan.maynouatupc.edu>
}




\seealso{
entropy, redundancy, correction.redundancy, correction.entropy
}
\examples{
require("MEET")
data(RenyiOrder)
data(iicc)
s<-4
correctionaprox(x=iicc$HXmax,matriu=iicc$TranscriptionFactor,s) 
}

