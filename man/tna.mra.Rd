\name{tna.mra}
\alias{tna.mra}

\title{
Master Regulator Analysis (MRA) over a list of regulons.
}

\description{
This function takes a TNA object and returns the results of the RMA analysis over a 
list of regulons from a transcriptional network (with multiple hypothesis testing corrections).
}

\usage{
tna.mra(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, tnet="dpi", verbose=TRUE)
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
compute the MRA analysis. Options: "dpi" and "ref".
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
a data frame in the slot "results", see 'rma' option in \code{\link[RTN:tna.get]{tna.get}}.
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

#run MRA analysis pipeline
rtna <- tna.mra(rtna)

#get results
tna.get(rtna,what="mra")

}
}

\keyword{RMA}
