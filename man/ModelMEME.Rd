\name{ModelMEME}
\alias{ModelMEME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MEME algortihm to detect TFBS within a sequence
}
\description{ MEME algorithm is used to construct a PWM model of the intput TFBS sequences. The returned values is the file with the motif

}
\usage{
ModelMEME(iicc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iicc}{options of the MEET program}
}
\details{
The output is the name of the file where the results are writen}
\value{
    output: 
  \item{nameMEME}{name of the MEME file}
  
  
}


\author{
Erola Pairo <epeiroatibec.pcb.ub.es>}

\seealso{
MEET, kfold.MEME}
