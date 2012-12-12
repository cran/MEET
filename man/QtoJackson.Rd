\name{QtoJackson}
\alias{QtoJackson}
\title{
Q to Jackson: transform a Q-residual into a confidence interval
}
\description{
Given Q-residual, this function uses the parameters calculated in JAcksonParameters to convert a Q-residual into a confidence interval}
\usage{
QtoJackson(Q,h0, x1, x2, x3)

}
\arguments{
  \item{Q}{Q-residual}
  \item{h0}{h0 parameter}
  \item{x1}{x1 parameter}
  \item{x2}{x2 parameter}
  \item{x3}{x3 parameter}
 
  
}
\value{
 output:confidence interval of the given residual.
}

\author{
Erola Pairo <epairoatibec.pcb.ub.es>}
\seealso{
kfold.PCA, Predict.PCA, JacksonParameters
}
\examples{
require("MEET")
data(TranscriptionFactor)
data(iicc)
nPCs<-1
TFBS<-iicc$Transcriptionfactor
Prob<-iicc$background
missing<-iicc$missing
NumericalMatrix<-numericalDNA(Prob)
suma<-apply(TFBS,2,function(y){sum(y=="-")})
threshold<-floor(nrow(TFBS)*missing/100)
TFBS<-TFBS[, suma<=threshold]
ncolTFBS<-ncol(TFBS)
TFBSnum<-apply(TFBS,1,function(x){as.vector(t(NumericalMatrix[x,]))})
TFBSnum<-t(TFBSnum)
model<-pca(TFBSnum, nPcs=2, method="svd", center=TRUE)
JacksonPars<-JacksonParameters(nPCs,TFBSnum)
Qres<-3.45
confidence<-QtoJackson(Qres,h0=1,x1=1,x2=1,x3=1)
}
