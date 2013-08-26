\name{tni.get}
\alias{tni.get}

\title{
Get information from individual slots in a TNI object.
}

\description{
Get information from individual slots in a TNI object and any available results from a previous analysis.
}

\usage{
tni.get(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE)
}

\arguments{
\item{object}{
an object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. Options: "gexp", "tfs", "para", 
"refnet", "tnet","refregulons","regulons","cdt", "summary" and "status". Regulons can also be retrieved mapped
to the assiged mode of action ("regulons.and.mode" or "refregulons.and.mode").
}
\item{order}{
a single logical value specifying whether or not the output data should be ordered by significance. 
Valid only for "cdt" option.
}
\item{ntop}{
a single integer value specifying to select how many results of top significance from "cdt" option.
}
\item{reportNames}{
a single logical value specifying to reporte regulators with 'names'
(when reportNames=TRUE) or not (when reportNames=FALSE). This option is effective 
only if regulators were named with alternative identifiers. It takes effect for 
the option "cdt".
}
}

\value{
get the slot content from an object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}

\author{
Mauro Castro
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

# check summary
tni.get(rtni,what="summary")

# get reference/relevance network
refnet<-tni.get(rtni,what="refnet")

# get transcriptional network
tnet<-tni.get(rtni,what="tnet")

# get status of the pipeline
tni.get(rtni,what="status")
}
}

\keyword{methods}
