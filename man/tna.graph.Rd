\name{tna.graph}
\alias{tna.graph}

\title{
Compute a graph from TNA objects.
}

\description{
Extract results from a TNA object and compute a graph.
}

\usage{
tna.graph(object, tnet = "dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL)
}

\arguments{
\item{object}{
an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{tnet}{
a single character value specifying which network information should be used to compute the graph. 
Options: "ref" and "dpi".
}
\item{gtype}{
a single character value specifying the graph type.
Options: "rmap" and "amap". The "rmap" option returns regulatory maps represented by TFs and targets 
(regulons) and "amap" computes association maps among regulons (estimates the overlap using the Jaccard Coefficient).
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
a single numeric value (>=0 and <=1) specifying the cutoff for association map filter. When 
amapFilter="phyper", amapCutoff corresponds to a pvalue cutoff; when amapFilter="quantile", 
amapCutoff corresponds to a quantile threshold; and when amapFilter="custom", amapCutoff is 
a JC threshold.
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

#run MRA analysis pipeline
rtna<-tni2tna.preprocess(rtni, phenotype=dt4rtn$pheno, hits=dt4rtn$hits, phenoIDs=dt4rtn$phenoIDs)
rtna <- tna.mra(rtna)

# compute regulatory maps
g<-tna.graph(rtna, tnet="dpi", gtype="rmap", tfs=tfs4test)

# option: plot the igraph object using RedeR
#library(RedeR)
#rdp<-RedPort()
#calld(rdp)
#addGraph(rdp,g)
#relax(rdp,p1=50,p5=20)

# compute association maps
g<-tna.graph(rtna, tnet="ref", gtype="amap", tfs=tfs4test)

}
}

\keyword{methods}
