\name{TPC-class}
\Rdversion{1.1}
\docType{class}
\alias{TPC-class}
\alias{tpc.countedges,TPC-method}
\alias{tpc.fuseNetwork,TPC-method}

\title{Class \code{"TPC"}: an S4 class for Partical Correlation Transcriptional Network Inference.}

\description{
This S4 class includes a series of methods to do transcriptional network inference 
for high-throughput gene expression.
}

\section{Objects from the Class}{
Objects can be created directly from objects of the class \code{\link[RTN:TNI-class]{TNI-class}} by calls of the form \code{new("TPC", tniclassObj=object)}.
}


\section{Slots}{
  \describe{
    \item{\code{gexp}:}{Object of class \code{"matrix"},
    a gene expression matrix.
    }
    \item{\code{transcriptionFactors}:}{Object of class \code{"char_Or_NULL"},
    a vector with transcription factor identifiers.
    }
    \item{\code{results.pcor}:}{Object of class \code{"list"}, a list of results (see the returned values in the function \code{\link[RTN:tni2tpc.corpcor]{tni2tpc.corpcor}}).} 
		\item{\code{status}:}{Object of class \code{"character"}, a character value specifying the status of the TNI object based on the available methods.
    }
  }
}

\section{Methods}{
  \describe{
    \item{tpc.countedges}{\code{signature(object = "TPC")}: see  \code{\link[RTN:tpc.countedges]{tpc.countedges}}}
    \item{tpc.fuseNetwork}{\code{signature(object = "TPC")}: see  \code{\link[RTN:tpc.fuseNetwork]{tpc.fuseNetwork}}}
	  }
}

\author{
Ines de Santiago
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

## see tni2tpc.corpcor
}

\keyword{classes}
