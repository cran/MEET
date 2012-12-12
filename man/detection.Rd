\name{detection}
\alias{detection}
\title{
Detection: A set of functions for detection of TFBS}
\description{
This function contains a set of functions for detection of Transcription Factor Binding Sites. Specifically, motif discovery algorithms are ITEME, MEME, MDscan, MATCH and Qresiduals.
}
\usage{
detection(iicc)
}
\arguments{
  \item{iicc}{Set of initial conditions for the MEET-package: mode, method, 
                 background,alignment,pvalor,parameters,Transcriptionfactor,
		 nummotif,lenmotif,sentit,position,missing,vector,gapopen,
		 maxiters,gapextend}
}
\details{
	The main parameter is method ( a set of motif discovery algorithms: 
	ITEME(Entropy, divergence), Q-residuals, MEME, Match and MDscan)}
\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>}

\seealso{
MEET, ConstructModel}
\examples{
require("MEET")
data(iicc)
detection(iicc)}
	