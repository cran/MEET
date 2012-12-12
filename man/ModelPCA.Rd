\name{ModelPCA}
\alias{ModelPCA}
\title{ PCA model for a set of TFBS}
\description{Performs a principal components analysis of the input DNA aligned sequences. It can be used to construct a model with the number of components chosen using the validation method, or the number of components entered by the user. The PCA model for the aligned sequences, and the parameters to calculate the Q-residuals statistics are returned.}
\usage{ModelPCA(iicc)}

\arguments{
  \item{iicc}{options described in the MEET program}
}

\details{
The specific options for this program are: Order (Number of Principal Components used),  and Missing (Percentage threshold of unknown nucleotides in a given position to take into account this position. Default 50)}

\value{
 output:list with the PCA model (same output than in pcaMEthods package), and the parameters of the model: the parameters needed to calculate the Jackon statistics, the numerical TFBS matrix and the dimensions of the TFBS matrix used once the missing values are estimated.
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es>
}

\seealso{
PredictPCA, kfold.PCA, MEET, detection}

\examples{
require("MEET")
data(iicc)
data(TranscriptionFactor)
iicc$parameterIdeal<-2
iicc$pvalor<-0.1
iicc$parameter<-2
TFBSmodel<-ModelPCA(iicc)
}