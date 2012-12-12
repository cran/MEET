\name{PredicttransMEME}
\alias{PredicttransMEME}
\title{MAST algorithm to detect TFBS in a sequence}
\description{ PWM in transfac format is used to find TFBS into a query sequence. }
\usage{PredicttransMEME(iicc,TF)}
\arguments{
  \item{iicc}{options of the MEET program}
 \item{TF}{Transcription Factor}}
\details{Uses the transMEME constructed motif}
\value{output: the detected sequences and its position within the sequence.}
\author{Erola Pairo <epeiroatibec.pcb.ub.es>}
\seealso{MEET, kfold.transMEME, ModeltransMEME}
