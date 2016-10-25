\name{tni.gsea2}
\alias{tni.gsea2}

\title{
Two-tailed Gene Set Enrichment Analysis (GSEA) over a list of regulons.
}

\description{
This function takes a TNI object and returns a CMAP-like analysis obtained by two-tailed GSEA over a list of regulons in a transcriptional network.
}

\usage{
tni.gsea2(object, minRegulonSize=15, doSizeFilter=FALSE, exponent=1, tnet="dpi",
  tfs=NULL, samples=NULL, features=NULL, refsamp=NULL, log=FALSE, verbose=TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from the analysis.
}
\item{doSizeFilter}{
a logical value. If TRUE, regulons are be filtered by size and eliminated if not 
larger than minRegulonSize.
}
\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see 'gseaScores' function at HTSanalyzeR).
}
\item{tnet}{
can take values of 'refnet', 'dpi' or 'cdt'. It refers to the version of the regulatory
network that will be used for GSEA analysis.
}
\item{tfs}{
an optional vector with transcription factor identifiers (this option overrides the 'stepFilter' argument).
}
\item{samples}{
an optional vector with sample numbers which will be exclusively used to compute the
GSEA.
}
\item{features}{
a string vector containing features for feature selection.
}
\item{refsamp}{
an optional string, containing the reference sample for differential expression calculations. If not given, the expression average will be used as a reference value.
}
\item{log}{
a logical value. If TRUE, differential expression calculations will be computed in log
space.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
a list with enrichment scores for all samples in the TNI.
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link[RTN:tna.gsea2]{tna.gsea2}}
\code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}}
}

\examples{

data(dt4rtn)

tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3")
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

rtni <- tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)

#run GSEA2 analysis pipeline
EScores <- tni.gsea2(rtni)

}
}

\keyword{GSEA}