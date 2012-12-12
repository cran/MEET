\name{ConstructModel}
\alias{ConstructModel}
\title{
A set of functions for training of motif discovery algorithms.
}
\description{
This function contains a set of functions for training of motif discovery algorithms. Specifically, the motif discovery algorithms are ITEME, MEME, MDscan, MATCH and Qresiduals.}
\usage{
ConstructModel(iicc, TF)
}

\arguments{
  \item{iicc}{Set of initial conditions for the MEET-package: mode, method, background,alignment,threshold,parameters,Transcriptionfactor,nummotif,lenmotif,sentit,position,missing,vector,gapopen,maxiters,gapextend}
  \item{TF}{A set of nucleotide sequences}
}

\details{This function has two output: validation scores and ROC curve.}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>}
\seealso{
MEET, detection}
