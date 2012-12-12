\name{kfold.PCA}
\alias{kfold.PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA}
\description{For a vector of principal components the Q-residuals of a DNA symbolical sequence are calculated using a PCA model of the numerical DNA motif. First sequences are converted to numerical DNA sequences and then the model is applied. A matrix is returned for each number of PC with the value of the Q-residuals in each position and a label indicating if the sequence belong to a binding site or not}
\usage{
kfold.PCA(iicc, TF)
}
\arguments{
  \item{iicc}{ Options described  in the MEET function}
  \item{TF}{DNA sequences used to construct the motif model}
}
\details{ Alignment method has to be installed in your computer. 

}
\value{
As a list, for each nPCs
  \item{matriuROC[[nPCs]] }{matrix with two colums, in the first one the Q-residuals for each studied sequence, and the second one indicates if the sequence belong to a TFBS}

}
\references{
Jolliffe I.T. Principal Component Analysis, Series: Springer Series in Statistics, 2nd ed., Springer, NY, 2002, XXIX, 487 p. 28 illus. ISBN 978-0-387-95442-4
Stacklies, Wolfram, Redestig, Henning, Scholz, Matthias, Walther, Dirk, and Selbig, Joachim: pcaMethods a bioconductor package providing PCA methods for incomplete data, Bioinformatics 23(9), volume 23, 1164-1167, 2007
}
\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}


\seealso{
MEET, kfold.Entropy, kfold.MATCH, kfold.Divergence, PCanalysis, llegir_DNA, convertDNA, numericalDNA. 
}
\examples{
require("MEET")
require("seqinr")
pathMEET <- system.file("sequences", package = "MEET")
data(iicc)
data(TranscriptionFactor)
iicc$method<-"PCA"
#Define the numner of principal components
iicc$vector<-c(1,3,5)
kfold.PCA(iicc, TF = paste(pathMEET, "AP1.fa", sep = "/"))
  

}
