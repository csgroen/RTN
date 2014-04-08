\name{tni.graph}
\alias{tni.graph}

\title{
Compute a graph from TNI objects.
}

\description{
Extract results from a TNI object and compute a graph.
}

\usage{
tni.graph(object, tnet = "dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", 
amapCutoff=NULL, mask=FALSE, ntop=NULL)
}

\arguments{
\item{object}{
an object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{tnet}{
a single character value specifying which network information should be used to compute the graph. 
Options: "ref" and "dpi".
}
\item{gtype}{
a single character value specifying the graph type.
Options: "rmap", "amap", "mmap" and "mmapDetailed". The "rmap" option returns regulatory maps represented by 
TFs and targets (regulons); "amap" computes association maps among regulons (estimates the overlap using the Jaccard Coefficient); "mmap" and "mmapDetailed" return modulated maps derived from the \code{\link[RTN:tni.conditional]{tni.conditional}} function.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon. Regulons with fewer than this number are 
removed from the graph.
}
\item{tfs}{
a vector with transcription factor identifiers.
}
\item{amapFilter}{
a single character value specifying which method should be used to filter association maps 
(only when gtype="amap"). Options: "phyper","quantile" and "custom".
}
\item{amapCutoff}{
a single numeric value (>=0 and <=1) specifying the cutoff for the association map filter. When 
amapFilter="phyper", amapCutoff corresponds to a pvalue cutoff; when amapFilter="quantile", 
amapCutoff corresponds to a quantile threshold; and when amapFilter="custom", amapCutoff is 
a JC threshold.
}
\item{mask}{
a logical value specifying to apply a mask on the 'amapFilter', keeping at least the 
best weighted edge (when verbose=TRUE) or not (when verbose=FALSE).
}
\item{ntop}{
when gtype="mmapDetailed", ntop is an optional single integer value (>=1) specifying the number of TF's targets that should be used to compute the modulated map. The n targets is derived from the top ranked TF-target interations, as defined in the mutual information analysis used to construct the regulon set.
}
}

\value{
a graph object.
}

\author{
Mauro Castro
}

\examples{

data(dt4rtn)

tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
rtni<-tni.dpi.filter(rtni, eps=0.05)

# compute regulatory maps
g<-tni.graph(rtni, tnet="dpi", gtype="rmap", tfs=tfs4test)

# option: plot the igraph object using RedeR
library(RedeR)
rdp<-RedPort()
calld(rdp)
addGraph(rdp,g)
addLegend.shape(rdp,g)
addLegend.color(rdp,g,type="edge")
relax(rdp,p1=50,p5=20)

# compute association maps
resetd(rdp)
g<-tni.graph(rtni, tnet="ref", gtype="amap", tfs=tfs4test)
addGraph(rdp,g)
addLegend.size(rdp,g)
addLegend.size(rdp,g,type="edge")

}
}

\keyword{methods}
