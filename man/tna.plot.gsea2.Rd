\name{tna.plot.gsea2}
\alias{tna.plot.gsea2}
\title{
Plot enrichment analyses from TNA objects.
}

\description{
This function takes a TNA object and plots the two-tailed GSEA results for individual regulons.
}

\usage{
tna.plot.gsea2(object, labPheno="tna", file=labPheno, filepath=".", regulon.order="size", 
ntop=NULL, tfs=NULL, ylimPanels=c(-1.0,3.0,-0.5,0.5), heightPanels=c(2.0,0.8,5.0), width=3.5, 
height=3.8, ylabPanels=c("Phenotype","Regulon","Enrichment score"), 
xlab="Position in the ranked list of genes", alpha=1.0, 
sparsity=10, autoformat=TRUE, ...)
}

\arguments{
\item{object}{
an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{file}{
a character string naming a file.
}
\item{filepath}{
a single character value specifying where to store GSEA2 figures.
}
\item{regulon.order}{
a single character value specifying whether regulons should be ordered by 'size', 
'score', 'pvalue', 'adj.pvalue' and 'name' (or 'none' to keep the input ordering).
}
\item{ntop}{
a single integer value specifying how many regulons of top significance will be plotted.
}
\item{tfs}{
an optional vector with transcription factor identifiers (this option overrides the 'ntop' argument).
}
\item{ylimPanels}{
a numeric vector of length=4 specifying y coordinates ranges of the 1st and 3th plots 
(i.e. ylim for 'Phenotypes' and 'Running enrichment score').
}
\item{heightPanels}{
a numeric vector of length=3 specifying the relative height of each panel in the plot.
}
\item{width}{
a single numeric value specifying the width of the graphics region in inches.
}
\item{height}{
a single numeric value specifying the height of the graphics region in inches.
}
\item{ylabPanels}{
a character vector of length=3 specifying the the title for the y axes.
}
\item{xlab}{
a single character specifying the the title for the x axis.
}
\item{labPheno}{
a single character specifying a label for the phenotype (will also be used as the name of output file).
}
\item{alpha}{
a single numeric value in [0,1] specifying the transparency of the hits in the ranked list.
}
\item{sparsity}{
a single integer value (>1) specifying the density of the dots representing the running score.
}
\item{autoformat}{
a single logical value specifying to set the graph format using predefined themes. This option 
overrides the "ylimPanels" argument.
}
\item{\dots}{
other arguments used by the function pdf.
}
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:tna.gsea2]{tna.gsea2}}
}

\examples{

data(dt4rtn)

tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

rtni <- tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
rtni<-tni.dpi.filter(rtni)
rtna<-tni2tna.preprocess(rtni, phenotype=dt4rtn$pheno, hits=dt4rtn$hits, phenoIDs=dt4rtn$phenoIDs)

# run GSEA2 analysis pipeline
rtna <- tna.gsea2(rtna, stepFilter=FALSE)

# plot available GSEA2 results
tna.plot.gsea2(rtna, labPheno="test")
}
}

\keyword{GSEA2}
