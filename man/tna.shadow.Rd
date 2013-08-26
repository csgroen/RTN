\name{tna.shadow}
\alias{tna.shadow}

\title{
shadow analysis over a list of regulons.
}

\description{
This function takes a TNA object and returns the results of the shadow analysis 
over a list of regulons in a transcriptional network
(with multiple hypothesis testing corrections).
}

\usage{
tna.shadow(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, minIntersectSize=1,
                  nPermutations=1000, exponent=1, tnet="ref", orderAbsValue=TRUE, stepFilter=TRUE,
                  tfs=NULL, verbose=TRUE)
}
\arguments{

\item{object}{
a preprocessed object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details).
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from the analysis.
}
\item{minIntersectSize}{
a single integer or numeric value specifying the minimum number of 
elements in the intersect between any two regulons in the shadow analysis 
(as percentage value).
}
\item{nPermutations}{
a single integer or numeric value specifying the number of permutations 
for deriving p-values in GSEA.
}
\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see 'gseaScores' function at HTSanalyzeR).
}
\item{tnet}{
a single character value specifying which transcriptional network should to used to 
compute the shadow and shadow analyses. Options: "dpi" and "ref".
}
\item{orderAbsValue}{
a single logical value indicating whether the values should be converted to absolute values and then 
ordered (if TRUE), or ordered as they are (if FALSE).
}
\item{stepFilter}{
a single logical value specifying to use a step-filter algorithm removing non-significant 
regulons derived from \code{\link[RTN:tna.gsea1]{tna.gsea1}} (when stepFilter=TRUE) or not 
(when stepFilter=FALSE). It may have a substantial impact on the overall processing time.
}
\item{tfs}{
an optional vector with transcription factor identifiers (this option overrides the 'stepFilter' argument).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
a data frame in the slot "results", see 'shadow' in \code{\link[RTN:tna.get]{tna.get}}.
}

\author{Mauro Castro}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}
\code{\link[RTN:tna.shadow]{tna.shadow}}
}

\examples{

data(dt4rtn)

tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

rtni <- tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
rtni<-tni.dpi.filter(rtni)
rtna<-tni2tna.preprocess(rtni, phenotype=dt4rtn$pheno, hits=dt4rtn$hits, phenoIDs=dt4rtn$phenoIDs)

#run overlap analysis pipeline
rtna <- tna.overlap(rtna)

#run shadow analysis pipeline
rtna <- tna.shadow(rtna,stepFilter=FALSE)

#get results
tna.get(rtna,what="shadow")

# run parallel version with SNOW package!
options(cluster=makeCluster(4, "SOCK"))
rtna <- tna.shadow(rtna)
stopCluster(getOption("cluster"))
}
}

\keyword{shadow}
