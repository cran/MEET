\name{kfold.transMEME}
\alias{kfold.transMEME}
\title{
Leave-one-out cross-validation for MEME/MAST through training.matrix aligned with MUSCLE or CLUSTALW.
}
\description{
This function does leave-one-out cross-validation for MEME/MAST.  In this case, a set of nucleotide sequences is lined up MUSCLE and CLUSTALW.  This is main difference between transMEME and MEME.
}
\usage{
kfold.transMEME(iicc, TF)
}

\arguments{
  \item{iicc}{
Set of initial conditions for the MEET-package: mode, method, background,alignment,threshold,parameters,Transcriptionfactor,nummotif,lenmotif,sentit,position,missing,vector,gapopen,maxiters,gapextend}
  \item{TF}{A set of nucleotide sequence}
}
\details{
This function needs MEME/MAST software.
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>
}

\seealso{MEET, kfold.Entropy, kfold.transMEME,kfoldMEME, kfold.Divergence, kfold.PCA}



