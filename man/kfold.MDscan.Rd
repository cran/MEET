\name{kfold.MDscan}
\alias{kfold.MDscan}
\title{
Leave-one-out cross-validation for MDscan.
}
\description{
Given a training sequence set, the optimal length and number motif has been estimated by means of leave-one-out cross-validation from length and number motif set.  From each value, the ROC curve has been calculated.  From this results, the optimal value has been considered according to the area under conver surface maximum
}
\usage{
kfold.MDscan(iicc, TF)
}

\arguments{
\item{iicc}{A set of inicial conditions for the MEET-package (mode, method, background, alignment, threshold, parameters, Transcriptionfactor, nummotif, lenmotif, sentit, position, missing, vector, gapopen, maxiters, gapextend)}
  \item{TF}{ A set of nucleotide sequences}
}

\author{
Erola Pairo <epeiroatibec.pcb.ub.es> and Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
kfold.Divergence, kfold.Entropy, kfold.MEME, kfold.MATCH and kfold.PCA
}


