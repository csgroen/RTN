\name{tna.overlap}
\alias{tna.overlap}

\title{
Overlap analysis over a list of regulons.
}

\description{
This function takes a TNA object and returns the results of the overlap analysis among 
regulons in a transcriptional network (with multiple hypothesis testing corrections).
}

\usage{
tna.overlap(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, tnet="ref", 
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
\item{tnet}{
a single character value specifying which transcriptional network should to used to 
compute the overlap analysis. Options: "dpi" and "ref".
}
\item{tfs}{
an optional vector with transcription factor identifiers.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
a data frame in the slot "results", see 'overlap' option in \code{\link[RTN:tna.get]{tna.get}}.
}

\author{Mauro Castro}

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

#run overlap analysis pipeline
rtna <- tna.overlap(rtna)

#get results
tna.get(rtna,what="overlap")
}
}

\keyword{overlap}
