\name{kfold.MEME}
\alias{kfold.MEME}
\title{
Leave-one-out cross-validation for MEME
}
\description{
Given a training sequence set, the optimal length and number motif has been estimated by means of leave-one-out cross-validation from length and number motif set.  From each value, the ROC curve has been calculated.  From this results, the optimal value has been considered according to the area under conver surface maximum
}
\usage{
kfold.MEME(iicc, TF)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
  \item{TF}{ A set of nucleotide sequences}
}
\details{
This function needs MEME/MAST software.
}
\references{
[5] T. Bailey and C. Elkan, Fitting a mixture model by expectation maximization to discover motifs in biopolymers, in Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology.AAAI Press, August 1994, pp. 28-36
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
kfold.PredictDivergence, kfold.PredictInformationTheory, MDscan, kfold.MATCH and kfold.PCA}


