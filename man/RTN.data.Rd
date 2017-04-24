\name{RTN.data}
\alias{dt4rtn}
\docType{data}

\title{A pre-processed dataset for the RTN package.}

\description{
A minimum dataset used to demonstrate RTN main functions.
}

\usage{data(dt4rtn)}

\format{
  \code{dt4rtn}{ List containing 6 R objects: 'gexp', 'gexpIDs', 'pheno', 'phenoIDs', 'hits' and 'tfs'.}
}

\details{

The dataset consists of 6 R objects used in the RTN vignettes. It should be regarded as a toy example for demonstration purposes only, despite being extracted, pre-processed and size-reduced from Fletcher et al. (2013) and Curtis et al. (2012).

\describe{
  \item{dt4rtn$gexp}{a named gene expression matrix with 250 samples (a subset from the Fletcher2013b package for demonstration purposes only).}
  \item{dt4rtn$gexpIDs}{a data.frame of characters with Probe-to-ENTREZ annotation.}
  \item{dt4rtn$pheno}{a named numeric vector with differential gene expression data.}
  \item{dt4rtn$phenoIDs}{a data.frame of characters with probe ids matching a secundary annotation source (e.g. Probe-to-ENTREZ).}
  \item{dt4rtn$hits}{a character vector with genes differentially expressed.}
  \item{dt4rtn$tfs}{a named vector with transcriptions factors.}
}
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Curtis C. et al., \emph{The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.} Nature 486, 7403. 2012.
}

\examples{
    data(dt4rtn)
}

\keyword{dataset}

