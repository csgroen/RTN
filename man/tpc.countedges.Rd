\name{tpc.countedges}
\alias{tpc.countedges}

\title{
Number of significant edges in a partial correlation network
}

\description{
The function counts the number of significant edges in a partial correlation network
}

\usage{
tpc.countedges(object, cutoff.method="prob", cutoff.tr = 0.9, verbose=TRUE)
}

\arguments{

  \item{object}{
  a processed object of class 'TPC' \code{\link[RTN:TPC-class]{TPC-class}} evaluated by the methods 
  }
  
  \item{cutoff.method}{
  one of 'qval' or 'prob' (default is 'prob')
  }

  \item{cutoff.tr}{
  treshold to determine significant edges. If cutoff.method is 'prob', significant edges will be those with prob > 
  cutoff.tr. If cutoff.method is 'qval', significant edges will be those with qval <= cutoff.tr
  }

  \item{verbose}{
  a single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE)
  }
}

\value{
Returns the number of significant edges in the network
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
tpc.countedges(rtni)
}

\keyword{methods}

