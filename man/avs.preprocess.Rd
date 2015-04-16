\name{avs.preprocess}
\alias{avs.preprocess}

\title{A preprocessing function for objects of class AVS.}

\description{
This is a generic function.
}

\usage{
avs.preprocess(object, nrand=1000, mergeColinked=TRUE, reldata="rel27CEU-NCBIB36", 
ldfilter="DprimeLOD", snpop=NULL, verbose=TRUE)
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
a single character value specifying the the HapMap Linkage Disequilibrium data used to build the associated variant set. Current option is release "rel27CEU-NCBIB36" (for additional details see the RTNdata package).
}
\item{ldfilter}{
a single character value specifying the parameter used to filter the HapMap Linkage Disequilibrium data. Current options are "DprimeLOD" and "Rsquare" (for additional details see the RTNdata package).
}
\item{snpop}{
an optional vector with rs# identifiers used to select random markers from the Linkage Disequilibrium dataset. This option is intended to set the universe size used to build the null distribution for the EVSE analysis (i.e. the random AVS).
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
# This example requires the RTNdata package!
# RTNdata is currently available under request <mauro.a.castro at gmail.com>
library(RTNdata)
data(bcarisk)
avs <- new("AVS", markers=bcarisk$rsid)
avs <- avs.preprocess(avs, nrand=100)
avs.get(avs)
}
}

\keyword{methods}

