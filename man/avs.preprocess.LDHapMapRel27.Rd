\name{avs.preprocess.LDHapMapRel27}
\alias{avs.preprocess.LDHapMapRel27}

\title{A preprocessing function for objects of class AVS.}

\description{
This function builds an associated variant set from the 'RTNdata.LDHapMapRel27' package.
}

\usage{
avs.preprocess.LDHapMapRel27(markers, nrand=1000, mergeColinked=TRUE, 
snpop="HapMapRel27", verbose=TRUE)
}

\arguments{
\item{markers}{
a data frame,  'BED file' format with rs# markers mapped to the same genome build of the LD source in the RTNdata package.
}
\item{nrand}{
a single integer value specifying the size to the random variant set.
}
\item{mergeColinked}{
a single logical value specifying to merge co-linked markers eventually 
present in more than one cluster in the AVS (when mergeColinked=TRUE) or 
not (when mergeColinked=FALSE).
}
\item{snpop}{
a single character value specifying the universe size used to build the null distributions for the random AVSs. Options: 'dbSNP' and 'HapMapRel27'. The 1st option represents random SNPs listed in the dbSNP database release 130, while the 2nd random SNPs listed in the HapMap release 27 project. Alternativelly, 'snpop' can be a customized 'BED file' data frame with rs# identifiers (for further details on the LD data, please see 'RTNdata.LDHapMapRel27' package).
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
library(RTNdata.LDHapMapRel27)
data(bcarisk)
avs <- avs.preprocess.LDHapMapRel27(bcarisk, nrand=100)
avs.get(avs)
}
}

\keyword{methods}

