\name{avs.preprocess}
\alias{avs.preprocess}

\title{A preprocessing function for objects of class AVS.}

\description{
This is a generic function.
}

\usage{
avs.preprocess(object, nrand=1000, mergeColinked=TRUE, reldata="RTNdata.LDHapMap.rel27", snpop="all", verbose=TRUE)
}

\arguments{
\item{object}{
an object. When this function is implemented as the S4 method of class 
\code{\link[RTN:AVS-class]{AVS-class}}, this argument is an object of class 'AVS'.
}
\item{nrand}{
a single integer value specifying the size to the random variant set.
}
\item{mergeColinked}{
a single logical value specifying to merge co-linked markers eventually 
present in more than one cluster in the AVS (when mergeColinked=TRUE) or 
not (when mergeColinked=FALSE).
}
\item{reldata}{
a single character value specifying the RTNdata package used to build the associated variant set. Current options are "RTNdata.LDHapMap.rel27" and "RTNdata.LD1000g.rel20130502" (for additional details on the Linkage Disequilibrium data, see the RTNdata package documentation).
}
\item{snpop}{
a single character value specifying the universe size used to build the null distributions for the VSE/EVSE analysis (i.e. random AVSs). Options: 'all' or 'ld'. The 1st option represents all the SNPs genotyped in the Linkage Disequilibrium (LD) dataset, while the 2nd represents all the SNPs observed after the LD threshold is applyed. Alternativelly, 'snpop' can be a customized 'BED file' data frame with rs# identifiers.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:AVS-class]{AVS-class}}
}

\examples{

\dontrun{
#This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDHapMap.rel27)
data(bcarisk)
avs <- new("AVS", markers=bcarisk)
avs <- avs.preprocess(avs, nrand=100)
avs.get(avs)
}
}

\keyword{methods}

