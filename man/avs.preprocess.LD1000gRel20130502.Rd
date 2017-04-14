\name{avs.preprocess.LD1000gRel20130502}
\alias{avs.preprocess.LD1000gRel20130502}

\title{A preprocessing function for objects of class AVS.}

\description{
This function builds an associated variant set from the 'RTNdata.LD1000gRel20130502' package.
}

\usage{
avs.preprocess.LD1000gRel20130502(markers, nrand=1000, mergeColinked=TRUE, 
snpop="1000g", ld_window_kb=200, ld_threshold_pval=1e-7, 
pAdjustMethod="bonferroni", ld_threshold_r2=NULL, verbose=TRUE)
}

\arguments{
\item{markers}{
a data frame, 'BED file' format with rs# markers mapped to the same genome build of the Linkage Disequilibrium (LD) source in the RTNdata package.
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
a single character value specifying the universe size used to build the null distributions for the random AVSs. Options: 'dbSNP' and '1000g'. The 1st option represents random SNPs listed in the dbSNP database release 147, while the 2nd random SNPs listed in the 1000 Genome project, release 20130502. Alternativelly, 'snpop' can be a customized 'BED file' data frame with rs# identifiers (for further details on the LD data, please see 'RTNdata.LD1000gRel20130502' package).
}
\item{ld_window_kb}{
a single integer value (<= 1000). This means only assess LD between markers within a specified window (e.g. within 1000 kb).
}
\item{ld_threshold_pval}{
a single numeric value specifying the cutoff for p-values considered significant. The approach derives Chi-square statistics based on R2 estimates and returns the associated p-values, which are then used to decide on the linkage disequilibrium between genetic markers.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be used (see 'p.adjust' for details).
}
\item{ld_threshold_r2}{
an optinal single numeric value (>= 0.2 and <=1.0) specifying the pairwise LD threshold between markers. 
(NOTE: this option overrides the 'ld_threshold_pval' statistics).
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
library(RTNdata.LD1000gRel20130502)
data(bcarisk)
avs <- avs.preprocess.LD1000gRel20130502(bcarisk, nrand=100)
avs.get(avs)
}
}

\keyword{methods}

