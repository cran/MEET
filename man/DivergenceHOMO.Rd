\name{DivergenceHOMO}
\alias{DivergenceHOMO}
\title{
DivergenceHOMO: Given a Transcription factor chooses the model for a specific organism and method.}
\description{
This function uses transcription factor name to choose the model for the Homo organism and Divergence algorithm.
}
\usage{
DivergenceHOMO(nameTF)
}
\arguments{
   \item{nameTF}{Transcription Factor name}
}

\details{
	The needed iicc is the vector of parameters to study }
\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>}

\seealso{
MEET,EntropyHOMO,EntropyDROSOPHILA,EntropyMUS,EntropyRATTUS,DivergenceHOMO,DivergenceMUS,DivergenceRATTUS,DivergenceMUS}
\examples{
require("MEET")
data(iicc)
DivergenceHOMO(nameTF=iicc$nameTF)}
