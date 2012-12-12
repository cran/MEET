\name{JacksonParameters}
\alias{JacksonParameters}
\title{
JacksonParameters: To calculates the parameters needed to transform a Q-residual to a confidence interval
}
\description{
Given PCA model, this function calculates the parameters to perform the Jackson statistics of the Q-residuals.}
\usage{
JacksonParameters(nPCs, TFBS)
}
\arguments{
  \item{nPCs}{ number of principal components of the model}
  \item{TFBS}{numerical TFBS matrix}
 
  
}
\value{
 output:list including the parameters h0 x1, x2 and x3 needed to calculate the Q-residuals statistics
}

\author{
Erola Pairo <epairoatibec.pcb.ub.es>}
\seealso{
kfold.PCA, Predict.PCA, QtoJackson
}
\examples{
require("MEET")
data(TranscriptionFactor)
data(iicc)
nPCs<-1
TFBS<-iicc$Transcriptionfactor
Prob<-iicc$background
missing<-50
NumericalMatrix<-numericalDNA(Prob)
suma<-apply(TFBS,2,function(y){sum(y=="-")})
threshold<-floor(nrow(TFBS)*missing/100)
TFBS<-TFBS[, suma<=threshold]
ncolTFBS<-ncol(TFBS)
TFBSnum<-apply(TFBS,1,function(x){as.vector(t(NumericalMatrix[x,]))})
TFBSnum<-t(TFBSnum)
model<-pca(TFBSnum, nPcs=nPCs, method="svd", center=TRUE)
JacksonPars<-JacksonParameters(nPCs,TFBSnum)
}
