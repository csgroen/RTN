\name{tpc.fuseNetwork}
\alias{tpc.fuseNetwork}

\title{
Fuse partial correlations networks
}

\description{
Given a series of networks, \code{tpc.fuseNetwork} will form a new 'fused network' by combining Z-scores edgewise using Stouffer's weighted Z-score method. Weights are taken from the number of samples used to infer each individual network. The function \code{fdrtool} from \code{fdrtool} package is then used to estimate the density-based FDR (q-value) and the posterior probability (1 - local FDRs) of each edge of the fused network.}

\usage{
tpc.fuseNetwork(objectlist, verbose=TRUE)
}

\arguments{

  \item{objectlist}{
  A \link[base:list]{list} of networks. Each network is a processed object of class \code{\link[RTN:TPC-class]{TPC-class}}
  }
  
  \item{verbose}{
  a single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
  }

}

\value{
An object of the class  \code{\link[RTN:TPC-class]{TPC-class}} with the fused transcriptional network in the slot 'results.pcor'. The slot 'results.pcor' contains a list of items describing the edges in the transcriptional network:
\itemize{
  \item{\code{tfmode} shows the direction of the regulation, either positive (1) or negative (-1)}
  \item{\code{prob} shows the posterior probabilities}
  \item{\code{qval} shows the FDR adjusted p-values}
  \item{\code{pval} shows the p-values inferred directly from the z-scores} 
  \item{\code{z.score} shows the z-scores inferred using Stouffer's weighted Z-score method.}
}
}

\author{Ines de Santiago}

\seealso{
\code{\link[RTN:TPC-class]{TPC-class}}
}

\examples{

data(dt4rtn)

# TODO!

}

\keyword{methods}

