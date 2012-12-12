\name{kfold.Divergence}
\alias{kfold.Divergence}
\title{
Leave-one-out cross-validation for parametric divergence (ITEME).
}
\description{
Given a training sequence set, the optimal value for parametric divergence has been estimated by means of leave-one-out cross-validation from q-value set. For each q-value, the ROC curve has been calculated. From this results, the optimal q-value has been considered according to the area under convex surface maximum.
}
\usage{
kfold.Divergence(iicc, TF)
}

\arguments{
  \item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
  \item{TF}{ A set of nucleotide sequences}
}
\details{
  This function integrates the Mutual information (Renyi Order equal 1) and parametric divergence (Renyi Order different 1). Moreover, it contains a set of function for the detection of transcription factor binding sites:correction.entropy.R, correction.redundancy.R, entropy.shannon.R, entropy.renyi.R, entropy.corrected.R, probability.R, CalculRedundancy.R, diff.instructions.R, redundancy.R, ROCmodel.R, detector_2nOrdre.R, pvalue.R.
}
\author{
Joan Maynou <joan.maynouatupc.edu>
}
\seealso{
kfold.Entropy, kfold.MEME, kfold.MDscan, kfold.MATCH and kfold.PCA
}
