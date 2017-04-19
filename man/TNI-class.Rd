\name{TNI-class}
\Rdversion{1.1}
\docType{class}
\alias{TNI-class}
\alias{tni.preprocess,TNI-method}
\alias{tni.permutation,TNI-method}
\alias{tni.bootstrap,TNI-method}
\alias{tni.dpi.filter,TNI-method}
\alias{tni.get,TNI-method}
\alias{tni.graph,TNI-method}
\alias{tni.conditional,TNI-method}
\alias{tni.gsea2,TNI-method}
\alias{tni2tna.preprocess,TNI-method}

\title{Class \code{"TNI"}: an S4 class for Transcriptional Network Inference.}

\description{
This S4 class includes a series of methods to do transcriptional network inference 
for high-throughput gene expression.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TNI", gexp, transcriptionFactors)}.
}

\section{Slots}{
  \describe{
    \item{\code{gexp}:}{Object of class \code{"matrix"},
    a gene expression matrix.
    }
    \item{\code{transcriptionFactors}:}{Object of class \code{"char_Or_NULL"},
    a vector with transcription factor identifiers.
    }
    \item{\code{modulators}:}{Object of class \code{"char_Or_NULL"},
    a vector with modulator identifiers.
    }
    \item{\code{annotation}:}{Object of class \code{"data.frame"},
    a data frame with probe-to-gene information.
    }
    \item{\code{para}:}{Object of class \code{"list"},
  	a list of parameters for transcriptional network inference. These parameters are 
    those listed in the functions \code{\link[RTN:tni.permutation]{tni.permutation}}, 
    \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} and \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}}.
    }
    \item{\code{results}:}{Object of class \code{"list"},
  	a list of results (see the returned values in the functions 
		\code{\link[RTN:tni.permutation]{tni.permutation}}).   
    }
    \item{\code{summary}:}{Object of class \code{"list"},
  	a list of summary information for \code{gexp}, \code{transcriptionFactors}, 
    \code{para}, and \code{results}.    
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the TNI object 
    based on the available methods.
    }
  }
}

\section{Methods}{
  \describe{
    \item{tni.preprocess}{\code{signature(object = "TNI")}: see  \code{\link[RTN:tni.preprocess]{tni.preprocess}}}
    \item{tni.permutation}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.permutation]{tni.permutation}} }
    \item{tni.bootstrap}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} }
    \item{tni.dpi.filter}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} }    
    \item{tni.conditional}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.conditional]{tni.conditional}} }   
    \item{tni.get}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.get]{tni.get}} }        
    \item{tni.graph}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.graph]{tni.graph}} } 
    \item{tni.gsea2}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni.gsea2]{tni.gsea2}} }
    \item{tni2tna.preprocess}{\code{signature(object = "TNI")}: see \code{\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}} }
	 }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}
}

\examples{

data(dt4rtn)

rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs)
}

\keyword{classes}
