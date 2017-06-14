\name{TNI.data}
\alias{stni}
\docType{data}

\title{A pre-processed TNI for demonstration purposes only.}

\description{
A minimum TNI object that can be used to demonstrate RTN functionalities.
}

\usage{data(stni)}

\format{
  \code{stni}{A TNI-class with a subset of samples and genes from 
  the Fletcher2013b package.}
}

\details{
The TNI consists of a TNI-class with a subsetted gene expression matrix and
reduced list of transcription factors. It should be regarded as a toy example 
for demonstration purposes only, despite being extracted, pre-processed and 
size-reduced from Fletcher et al. (2013) and Curtis et al. (2012).
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Curtis C. et al., \emph{The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.} Nature 486, 7403. 2012.
}

\value{
a TNI-class.
}

\examples{
    data(stni)
}

\keyword{dataset}

