\name{run.read.MDscan}
\alias{run.read.MDscan}
\title{
Run and read MDscan on validation 
}
\description{
This function runs and reads MDscan on validation mode.}
\usage{
run.read.MDscan(input, k, len_motif, num_motif,call.MDscan)
}
\arguments{
  \item{input}{DNA sequence}
  \item{k}{Length DNA candidate seqeunce}
  \item{len_motif}{Length motif}
  \item{num_motif}{Number motif}
  \item{call.MDscan}{Path MDscan}
}
\details{
This function needs seqinr-Package and MDscan software.}


\author{Erola Pairo <epeiroatibec.pcb.ub.es>  and Joan Maynou <joan.maynouatupc.edu>}


\seealso{MEET, kfold.PredictInformationTheory, MEME, kfold.PredictDivergence, kfold.PCA}


