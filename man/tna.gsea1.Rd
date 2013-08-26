\name{tna.gsea1}
\alias{tna.gsea1}

\title{
One-tailed Gene Set Enrichment Analysis (GSEA) over a list of regulons.
}

\description{
This function takes a TNA object and returns the results of the GSEA analysis over a list of regulons in a transcriptional network (with multiple hypothesis testing corrections).
}

\usage{
tna.gsea1(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, nPermutations=1000, 
        exponent=1, tnet="dpi", orderAbsValue=TRUE, stepFilter=TRUE, tfs=NULL, verbose=TRUE)
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
compute the GSEA analysis. Options: "dpi" and "ref".
}
\item{orderAbsValue}{
a single logical value indicating whether the values should be converted to absolute values and then 
ordered (if TRUE), or ordered as they are (if FALSE).
}
\item{stepFilter}{
a single logical value specifying to use a step-filter algorithm removing 
non-significant regulons derived from \code{\link[RTN:tna.mra]{tna.mra}}
(when stepFilter=TRUE) or not (when stepFilter=FALSE). It may have a substantial impact 
on the overall processing time.
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
a data frame in the slot "results", see 'gsea1' option in \code{\link[RTN:tna.get]{tna.get}}.
}

\author{
Mauro Castro, Xin Wang
}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}
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

#run GSEA1 analysis pipeline
rtna <- tna.gsea1(rtna,stepFilter=FALSE)

#get results
tna.get(rtna,what="gsea1")

# run parallel version with SNOW package!
options(cluster=makeCluster(4, "SOCK"))
rtna <- tna.gsea1(rtna)
stopCluster(getOption("cluster"))
}
}

\keyword{GSEA}
