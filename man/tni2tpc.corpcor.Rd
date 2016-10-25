\name{tni2tpc.corpcor}
\alias{tni2tpc.corpcor}

\title{
Inference of transcriptional network using partial correlations
}

\description{
The function uses the 'pcor.shrink' method in the \code{corpcor} R package to compute shrinkage estimates of partial correlations between TFs and target genes. The function 'fdrtool' from \code{fdrtool} package is used to estimate parameters of the null
distribution adaptively from the data. Returns for each inferred edge in the transcriptional network the p-value, the density-based FDR (q-value) and the posterior probability (1 - local FDRs).
}

\usage{
tni2tpc.corpcor(object, verbose=TRUE)
}

\arguments{

  \item{object}{
  a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
  }
  
  \item{verbose}{
  a single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
  }

}

\value{
An object of the class  \code{\link[RTN:TPC-class]{TPC-class}} with the transcriptional network in the slot 'results.pcor'. The slot 'results.pcor' contains a list of items describing the edges in the transcriptional network:
\itemize{
  \item{\code{tn} shows the partial correlation scores}
  \item{\code{tfmode} shows the direction of the regulation, either positive (1) or negative (-1)}
  \item{\code{prob} shows the posterior probabilities}
  \item{\code{qval} shows the FDR adjusted p-values}
  \item{\code{pval} shows the p-values} 
  \item{\code{z.score} shows the z-scores inferred directly from the p-values.}
}
}



\author{Ines de Santiago}

\seealso{
\code{\link[RTN:TPC-class]{TPC-class}}
}

\examples{

data(dt4rtn)

# create a new TNI object
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs)

# run
rtni <- tni.preprocess(rtni)
rtni <- tni2tpc.corpcor(rtni)
}

\keyword{methods}

