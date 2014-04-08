\name{tni2tna.preprocess}
\alias{tni2tna.preprocess}

\title{A preprocessing function for objects of class TNI.}

\description{
This is a generic function.
}

\usage{
tni2tna.preprocess(object, phenotype=NULL, hits=NULL, phenoIDs=NULL, 
duplicateRemoverMethod="max", verbose=TRUE)
}

\arguments{
\item{object}{
a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
\code{\link[RTN:tni.permutation]{tni.permutation}}, \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} 
and \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}}.
}
\item{phenotype}{
a numeric or integer vector of phenotypes named by gene identifiers. 
Required for gsea, synergy and shadow methods (see \code{\link[RTN:tna.gsea1]{tna.gsea1}}).
}
\item{hits}{
a character vector of gene identifiers for those considered as hits.
Required for \code{\link[RTN:tna.mra]{tna.mra}} and \code{\link[RTN:tna.overlap]{tna.overlap}} methods.
}
\item{phenoIDs}{
an optional 2cols-matrix used to aggregate genes in the 'phenotype' (e.g. probe-to-gene ids; 
in this case, col 1 should correspond to probe ids). 
}
\item{duplicateRemoverMethod}{
a single character value specifying the method to remove the duplicates. The current version provides 
"min" (minimum), "max" (maximum), "average". Further details in 'duplicateRemover' function at the 
HTSanalyzeR package.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) or not 
(when verbose=FALSE).
}
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:TNI-class]{TNI-class}}
  \code{\link[RTN:TNA-class]{TNA-class}}
}

\examples{

data(dt4rtn)

tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
rtni<-tni.dpi.filter(rtni)
rtna<-tni2tna.preprocess(rtni, phenotype=dt4rtn$pheno, hits=dt4rtn$hits, phenoIDs=dt4rtn$phenoIDs)
}
}


\keyword{methods}
