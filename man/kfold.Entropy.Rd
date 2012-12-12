\name{kfold.Entropy}
\alias{kfold.Entropy}

\title{
Leave-one-out cross-validation for Renyi entropy (ITEME)
}
\description{
Given a training sequence set, the optimal value for Renyi entropy has been estimated by means of leave-one-out cross-validation from q-value set. For each q-value, the ROC curve has been calculated. From this results, the optimal q-value has been considered according to the area under convex surface maximum.
}
\usage{
kfold.Entropy(iicc, TF)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
  \item{TF}{ A set of nucleotide sequences}
}
\details{
 This function integrates the Shannon entropy for Renyi Order equal 1. Moreover, it contains a set of function for the detection of transcription factor binding sites:correction.entropy.R, correction.redundancy.R, entropy.Shannon.R,entropy.Renyi.R, entropy.corrected.R, probability.R, CalculRedundancy.R, diff.instructions.R, redundancy.R, ROCmodel.R, detector_1rOrdre_diff.R, pvalue.R.
}


\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
kfold.Divergence, kfold.MEME, kfold.MDscan, kfold.MATCH and kfold.PCA
}


