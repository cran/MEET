\name{ModelMATCH}
\alias{ModelMATCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Match algorithm to detect TFBS in a sequence
}
\description{ Match algorithm is used to construct a PWM model of the intput TFBS sequences. The returned values are the PWM model and the parameters needed to predict using MATCH and the PWM, as the position of the Core, the minimum and maximum scores for a query sequence and its core, the Corecut parameter and the dimensions of teh TFBS matrix. 

}
\usage{
ModelMATCH(iicc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iicc}{options of the MEET program}
}
\details{
The specific parameters for this detection are: the transcription factor to model and the background probabilities for each nucleotide}
\value{
    output: a list with the model and the parameters
  \item{model }{logodds matrix using information per site}
  \item{parametersModel }{posCore=position of the core, minim_core=minimum score for the core, maxim_core= maximum score for the core, minim=minimu score, maxim=maximum score, Corecut= percentage of core score that we want to consider, core=logodds for the core, ncolTFBS=dimension of the TFBS matrix}
  
}
\references{
E. Kel , E. Gossling , I. Reuter , E. Cheremushkin , O.V. Kel-Margoulis , and E. Wingender MATCHTM: a tool for searching transcription factor binding sites in DNA sequences Nucl. Acids Res. 31: 3576-3579.}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es>}

\seealso{
MEET, kfold.MATCH, ModelPCA}
\examples{
require("MEET")
data(iicc)
data(TranscriptionFactor)
iicc$method<-"MATCH"

iicc$pvalor<-0.8
iicc$parameter<-60
ModelMATCH(iicc)
}
