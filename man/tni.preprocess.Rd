\name{tni.preprocess}
\alias{tni.preprocess}

\title{A preprocessing function for objects of class TNI.}

\description{
This is a generic function.
}

\usage{
tni.preprocess(object, gexpIDs=NULL, cvfilter=TRUE, verbose=TRUE)
}

\arguments{
\item{object}{
an object. When this function is implemented as the S4 method of class 
\code{\link[RTN:TNI-class]{TNI-class}}, this argument is an object of class 'TNI'.
}
\item{gexpIDs}{
an optional data frame or matrix with probe-to-gene annotation. Column 1 must provide all probe ids listed in the 'gexp' matrix. Ideally, col1 = PROBEID, col2 = GENEID, and col3 = SYMBOL.
}
\item{cvfilter}{
a single logical value specifying to remove duplicated genes in the gene expression matrix using the probe-to-gene annotation. In this case, 'gexpIDs' must be provided, with col1 = PROBEID and col2 = GENEID. The decision is made based on the maximum dinamic range (i.e. keeping the probes with max coefficient of variation across all samples).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(dt4rtn)

rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs)
rtni <- tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
}

\keyword{methods}

