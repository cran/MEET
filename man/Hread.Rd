\name{Hread}
\alias{Hread}
\title{To read Entropy values
}
\description{ This function reads entropy values saved in memory.  From the entropy values, Hread calculates the variation of the total entropy
when the candidate sequence is added to the set.}
\usage{
Hread(training.set.mes.rand,val.set,iicc)
}

\arguments{
\item{training.set.mes.rand}{ A set of nucleotide sequences}
\item{val.set}{A candidate sequence}
\item{iicc}{A set of inicial conditions for the MEET-package}
 
}

\author{
Joan Maynou <joan.maynouatupc.edu>
}

\seealso{
Hmemory}

\examples{
require("MEET")
require("seqinr")
write.fasta <- get("write.fasta",pos="package:seqinr")
read.fasta <- get("read.fasta",pos="package:seqinr")
data(iicc)
training.set.mes.rand<-iicc$Transcriptionfactor
val.set<-c("A","T","C","C","G","T","A","T","T","A","C","G")
test<-Hread(training.set.mes.rand, val.set, iicc)
}

