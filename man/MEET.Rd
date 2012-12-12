\name{MEET}
\alias{MEET}
\title{
MEET: Motif Elements Estimation Toolkit}
\description{
MEET (Motif Elements Estimation Toolkit) is a R-package that integrates a set of computational algorithms for the detection of Transcription Factor Binding Sites (TFBS)}
\usage{
MEET(TF,nameTF, seqin, alg, method, mode, org, vector, num_motif, len_motif, direction, threshold, order, model,position, mv, gapopen, maxiters, gapextend,optionsFile)
}

\arguments{
  \item{TF}{ A set of nucleotide sequences}
  \item{nameTF}{Transcription Factor name}
  \item{seqin}{Candidate sequence}
  \item{alg}{A set of Multiple Sequence Alignment}
  \item{method}{A set of motif discovery algorithms}
  \item{mode}{Training or detection}
  \item{org}{Background organism}
  \item{vector}{Characteristic parameters methodfor leave-one-out cross-validation.}
  \item{num_motif}{Number of motif}
  \item{len_motif}{Length of motif}
  \item{direction}{Forward, reverse or both}
  \item{threshold}{P-value}
  \item{order}{Renyi Order}
  \item{model}{model}
  \item{position}{Binding site position for training}
  \item{mv}{missing values}
  \item{gapopen}{The gap open score}
  \item{maxiters}{Maximum number of iterations}
  \item{gapextend}{The gap extend score}
  \item{optionsFile}{Path for MEME/MAST, clustalx and muscle}
}
 


\details{
This function need the next packages: seqinr, fields, pcaMethods, Matrix, ROCR, Hmisc, KernSmooth. Moreover, it needs the next software: MEME/MAST Version 4.4.0, MATCH Version 1.0, MDscan, ClustalW and Muscle Version 3.8.}

