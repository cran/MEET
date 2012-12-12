\name{Alignment}
\alias{Alignment}
\title{To line up Transcription Factor Binding sites through Multiple Sequence Alignment (MSA)}
\description{This function, Alignment, lines up Transcription Factor Binding Sites, TFBS, through MSA. There are different kinds of MSA: MEME, Clustalw and Muscle}

\usage{Alignment(TF,iicc)}
\arguments{
\item{TF}{Transcription Factor}
\item{iicc}{A set of initial conditions for the MEET-package}
	    }

\author{
Joan Maynou <joan.maynouatupc.edu> and Erola Pairo  <epeiroatibec.pcb.ub.es>
}
\examples{
require("MEET")
data(iicc)
pathMEET <- system.file("sequences", package = "MEET")
Alignment(TF=paste(pathMEET, "AP1.fa", sep = "/"),iicc)
}

