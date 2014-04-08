\name{tna.get}
\alias{tna.get}

\title{
Get information from individual slots in a TNA object.
}

\description{
Get information from individual slots in a TNA object. Available results from a previous analysis can be selected either by pvalue cutoff (default) or top significance.
}

\usage{
tna.get(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, 
idkey=NULL)
}

\arguments{
\item{object}{
an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. 
Options: 'tnet', 'refnet', 'tfs', 'pheno', 'regulons','refregulons', 'para', 'mra', 'gsea1', 
'gsea2','overlap', 'synergy', 'shadow', 'summary' and 'status'. Regulons can also be retrieved 
mapped to the available phenotype vector ('regulons.and.pheno' or 'refregulons.and.pheno') or
mapped to the assiged mode of action ('regulons.and.mode' or 'refregulons.and.mode').
}
\item{order}{
a single logical value specifying whether or not the output data should be ordered by significance. 
Valid only for 'gsea1', 'gsea2', 'overlap', 'synergy' or 'shadow' options.
}
\item{ntop}{
a single integer value specifying to select how many results of top significance from 'gsea', 'overlap', 
'synergy' or 'shadow' options.
}
\item{reportNames}{
a single logical value specifying to report regulons with 'names' (when reportNames=TRUE) or not (when reportNames=FALSE). This option is effective only if transcription factors were named with alternative identifiers in the pre-processing analysis. It takes effect on 'mra', 'gsea', 'overlap', 'synergy' and 'shadow' options.
}
\item{idkey}{
an optional single character value specifying an ID name from the available 'TNA' annotation to be used as alias for data query outputs (obs. it has no effect on consolidated tables).
}
}

\value{
get the slot content from an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}

\author{
Mauro Castro
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

# run MRA analysis pipeline
rtna <- tna.mra(rtna)

# check summary
tna.get(rtna,what="summary")

# get results, e.g., MRA analysis
tna.get(rtna,what="mra")
}
}


\keyword{GSEA}
