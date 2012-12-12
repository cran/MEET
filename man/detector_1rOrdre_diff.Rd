\name{detector_1rOrdre_diff}
\alias{detector_1rOrdre_diff}
\title{Detection of Transcription Factor Binding Sites Through Differential Renyi Entropy}
\description{
This detection algorithm is based on Information Theory.  Specifically, entropy algorithm uses a parametric uncertainty measurement called Renyi entropy.  This computational method evaluates the variation on the total Renyi entropy of a set of sequences when a candidate sequence is assumed to be a a true binding site belonging to the set.The measurement of the variation of the total redundancy when the candidate sequence is added to the set has been computed by using the difference between the redundancy profile. This technique assumes independecy between positions in the binding sequence. 
}
\usage{
detector_1rOrdre_diff(training.set, val.set, iicc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.set}{A set of aligned nucleotide sequences}
  \item{val.set}{A candidate sequence}
  \item{iicc}{A set of inicial conditions for the MEET-package: mode, method,background,alignment,threshold,parameters,Transcriptionfactor,
  		 nummotif,lenmotif,sentit,position,missing,vector,gapopen,maxiters,gapextend}
}
\details{Options parameter has to contain the next arguments: maximum entropy (HXmax), correction entropy and redundancy from the Finite 
		 Sample Size Effect(correction, Redundancia_corregida,
		 Herror and ErrorHX)}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
detector_2Ordre, MEME, MDscan, Q-residuals and MATCH
}
\examples{
require("MEET")
data(iicc)
test<-detector_1rOrdre_diff(training.set=iicc$Transcriptionfactor,val.set=NULL,iicc) }

