\name{avs.plot1}
\alias{avs.plot1}

\title{
Plot results from AVS methods, single plots.
}

\description{
This function takes an AVS object and plots results from the VSE and EVSE methods.
}

\usage{
avs.plot1(object, what="vse", fname=what, ylab="genomic annotation", 
xlab="Number of clusters mapping to genomic annotation", breaks="Sturges", 
maxy=200, pValueCutoff=1e-2, width=8, height=3)
}

\arguments{
\item{object}{
an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}
\item{what}{
a character value specifying which analysis should be used. Options: "vse" and "evse".
}
\item{fname}{
a character value specifying the name of output file.
}
\item{ylab}{
a character value specifying the y-axis label.
}
\item{xlab}{
a character value specifying the x-axis label.
}
\item{breaks}{
breaks in the histogram, see \code{\link[graphics:hist]{hist}} function.
}
\item{maxy}{
a numeric value specifying the max y-limit.
}
\item{pValueCutoff}{
a numeric value specifying the cutoff for p-values considered significant.
}
\item{width}{
a numeric value specifying the width of the graphics region in inches.
}
\item{height}{
a numeric value specifying the height of the graphics region in inches.
}

}

\author{
Mauro Castro
}

\examples{

###

}

\keyword{VSE}
