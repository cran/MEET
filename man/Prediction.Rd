\name{Prediction}
\alias{Prediction}
\title{To detect Transcription Factor Binding sites by means of a model}
\description{This function, Prediction, detects Transcription Factor Binding Sites, TFBS, from the model. There are different kinds of model: Entropy, Divergence, PCA, MEME, MDscan, MATCH and transMEME.  This model is obtained from ConstructModel function.}
\usage{
Prediction(iicc)
}
\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package}
}

\author{
Joan Maynou <joan.maynouatupc.edu> and Erola Pairo  <epeiroatibec.pcb.ub.es>
}
\examples{
require("MEET")
data(iicc)
Prediction(iicc)
}

