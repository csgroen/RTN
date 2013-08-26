\name{TNA-class}
\Rdversion{1.1}
\docType{class}
\alias{TNA-class}
\alias{tna.mra,TNA-method}
\alias{tna.overlap,TNA-method}
\alias{tna.gsea1,TNA-method}
\alias{tna.gsea2,TNA-method}
\alias{tna.synergy,TNA-method}
\alias{tna.shadow,TNA-method}
\alias{tna.get,TNA-method}
\alias{tna.graph,TNA-method}

\title{Class \code{"TNA"}: an S4 class for Transcriptional Network Analysis.}

\description{
This S4 class includes a series of methods to do enrichment, synergy, shadow and overlap 
analyses in transcriptional networks.
}

\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{new("TNA", referenceNetwork, transcriptionalNetwork, transcriptionFactors, phenotype, hits)}.
}

\section{Slots}{
  \describe{
    \item{\code{referenceNetwork}:}{Object of class \code{"matrix"}, 
    an optional partial co-expression matrix.
    }
    \item{\code{transcriptionalNetwork}:}{Object of class \code{"matrix"}, 
    a partial co-expression matrix.
    }
    \item{\code{transcriptionFactors}:}{Object of class \code{"char_Or_NULL"}, 
    a vector of transcription factors.
    }
    \item{\code{phenotype}:}{Object of class \code{"num_Or_int"}, a numeric or integer 
    vector of phenotypes named by gene identifiers. 
    }
    \item{\code{hits}:}{Object of class \code{"character"}, a character vector of 
    gene identifiers for those considered as hits. 
    }
    \item{\code{annotation}:}{Object of class \code{"data.frame"},
    a data frame with transcriptional network annotation.
    }
    \item{\code{listOfReferenceRegulons}:}{Object of class \code{"list"}, a list of regulons 
    derived from the referenceNetwork.
    }
    \item{\code{listOfRegulons}:}{Object of class \code{"list"}, a list of regulons 
    derived from the transcriptionalNetwork (a 'regulon' is a vector of genes or potential 
    transcription factor targets). 
    }    
    \item{\code{listOfModulators}:}{Object of class \code{"list"}, a list of modulators 
    derived from the \code{\link[RTN:tni.conditional]{tni.conditional}} analysis.
    }
    \item{\code{para}:}{Object of class \code{"list"}, a list of parameters for 
    transcriptional network analysis. These parameters are those listed in the functions 
    \code{\link[RTN:tna.mra]{tna.mra}}, \code{\link[RTN:tna.overlap]{tna.overlap}},
    \code{\link[RTN:tna.gsea1]{tna.gsea1}}, \code{\link[RTN:tna.synergy]{tna.synergy}}, 
    \code{\link[RTN:tna.shadow]{tna.shadow}} and \code{\link[RTN:tna.gsea2]{tna.gsea2}}.
    }
    \item{\code{results}:}{Object of class \code{"list"},
    a list of results (see return values in the functions 
		\code{\link[RTN:tna.mra]{tna.mra}},\code{\link[RTN:tna.gsea1]{tna.gsea1}},
    \code{\link[RTN:tna.overlap]{tna.overlap}}, \code{\link[RTN:tna.synergy]{tna.synergy}},
    \code{\link[RTN:tna.shadow]{tna.shadow}} and \code{\link[RTN:tna.gsea2]{tna.gsea2}})
    }
    \item{\code{summary}:}{Object of class \code{"list"},
  	a list of summary information for \code{transcriptionalNetwork}, 
    \code{transcriptionFactors}, \code{phenotype},\code{listOfRegulons},
    \code{para}, and \code{results}.    
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a single character value specifying the status of the TNI object 
    based on the available methods.
    }    
  }
}

\section{Methods}{
  \describe{
    \item{tna.mra}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.mra]{tna.mra}} }    
    \item{tna.overlap}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.overlap]{tna.overlap}} }
    \item{tna.gsea1}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.gsea1]{tna.gsea1}} }
    \item{tna.gsea2}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.gsea2]{tna.gsea2}} }
    \item{tna.synergy}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.synergy]{tna.synergy}} }    
    \item{tna.shadow}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.shadow]{tna.shadow}} }    
    \item{tna.get}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.get]{tna.get}} }
    \item{tna.graph}{\code{signature(object = "TNA")}: see \code{\link[RTN:tna.graph]{tna.graph}} }    
	 }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}.
}

\examples{

data(dt4rtn)

rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs)

}

\keyword{classes}
