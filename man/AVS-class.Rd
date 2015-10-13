\name{AVS-class}
\Rdversion{1.1}
\docType{class}
\alias{AVS-class}
\alias{avs.preprocess,AVS-method}
\alias{avs.vse,AVS-method}
\alias{avs.evse,AVS-method}
\alias{avs.get,AVS-method}

\title{Class \code{"AVS"}: an S4 class for variant set enrichment analysis.}

\description{
This S4 class includes a series of methods to do enrichment analyses in Associated Variant Sets (AVSs).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AVS", markers)}.
}

\section{Slots}{
  \describe{
    \item{\code{markers}:}{Object of class \code{"character"}, a character vector with rs# markers.
    }
    \item{\code{validatedMarkers}:}{Object of class \code{"data.frame"}, 
    a data frame with genome positions of the validated markers.
    }  
    \item{\code{variantSet}:}{Object of class \code{"list"}, 
    an associated variant set.
    }
    \item{\code{randomSet}:}{Object of class \code{"list"}, 
    a random associated variant set. 
    }
    \item{\code{para}:}{Object of class \code{"list"},
    a list of parameters for variant set enrichment analysis.
    }
    \item{\code{results}:}{Object of class \code{"list"},
    a list of results (see return values in the AVS methods).
    }
    \item{\code{summary}:}{Object of class \code{"list"},
    a list of summary information for \code{markers}, \code{para}, 
    and \code{results}.    
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the AVS object 
    based on the available methods.
    }    
  }
}

\section{Methods}{
  \describe{
    \item{avs.preprocess}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.preprocess]{avs.preprocess}} } 
    \item{avs.vse}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.vse]{avs.vse}} } 
    \item{avs.evse}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.evse]{avs.evse}} }
    \item{avs.get}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.get]{avs.get}} }    
	 }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}.
}

\examples{

\dontrun{
#This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDrel27)
data(bcarisk)
avs <- new("AVS", markers=bcarisk$rsid)
}
}

\keyword{classes}
