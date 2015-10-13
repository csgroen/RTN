\name{avs.get}
\alias{avs.get}

\title{
Get information from individual slots in an AVS object.
}

\description{
Get information from individual slots in an AVS object.
}

\usage{
avs.get(object, what="summary", pValueCutoff=NULL)
}

\arguments{
\item{object}{
an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. 
Options: 'markers', 'validatedMarkers', 'variantSet', 'randomSet', 'linkedMarkers', 'randomMarkers', 'vse', 'evse', 'summary' and 'status'.
}
\item{pValueCutoff}{
an optional single numeric value specifying the cutoff to retrive results for p-values considered significant.
}

}

\value{
get the slot content from an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}

\author{
Mauro Castro
}

\examples{

\dontrun{
#This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDrel27)
data(bcarisk)
avs <- new("AVS", markers=bcarisk)
avs <- avs.preprocess(avs, nrand=100)
avs.get(avs)
}
}

\keyword{methods}

