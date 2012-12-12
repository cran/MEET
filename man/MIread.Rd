\name{MIread}
\alias{MIread}
\title{To read PredictDivergence values
}
\description{ This function reads divergence values saved in memory.  From the divergence values, MIread calculates the variation of the total divergence
when the candidate sequence is added to the set.}
\usage{
MIread(training.set,val.set,iicc)
}

\arguments{
\item{training.set}{ A set of nucleotide sequences}
\item{val.set}{A candidate sequence}
\item{iicc}{A set of inicial conditions for the MEET-package}
 
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
MImemory}

\examples{
require("MEET")
data(DrosophilaDivergence)
model<-list()
model$D<-iicc[["a1"]]$model$parameterModel$D
model$HXmax<-iicc[["a1"]]$model$parameterModel$HXmax
model$correctioc_1rOrdre<-iicc[["a1"]]$model$parameterModel$correction_1rOrdre
model$Entropy<-iicc[["a1"]]$model$parameterModel$HX
model$Mperfil<-iicc[["a1"]]$model$parameterModel$Mperfil
model$interA<-iicc[["a1"]]$model$parameterModel$interA
model$interB<-iicc[["a1"]]$model$parameterModel$interB
model$Divergence<-iicc[["a1"]]$model$model

test<-MIread(training.set=iicc[["a1"]]$Transcriptionfactor, val.set=iicc[["a1"]]$Transcriptionfactor[1,],iicc=model)
}

