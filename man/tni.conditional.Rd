\name{tni.conditional}
\alias{tni.conditional}

\title{
Modulators of transcription factor (TF) activity assessed by conditional mutual information analysis.
}

\description{
This function takes a TNI object and a list of candidate modulators, and computes the conditional mutual 
information over the TF-target interactions in a transcriptional network (with multiple hypothesis 
testing corrections). For each TF, the method measures the change in the mutual information between 
the TF and its targets conditioned to the gene expression of a modulator.
}

\usage{
tni.conditional(object, modulators=NULL, tfs=NULL, sampling=35, pValueCutoff=0.01, 
pAdjustMethod="bonferroni", minRegulonSize=15, minIntersectSize=5, 
miThreshold="md", prob=0.95, pwtransform=FALSE, verbose=TRUE)
}

\arguments{

\item{object}{
a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
\code{\link[RTN:tni.permutation]{tni.permutation}}, \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} 
and \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}}.
}
\item{modulators}{
a vector with gene identifiers for those considered as candidate modulators. If NULL, the function 
will assess all TFs in the network, provided that the candidate passes all filtering steps.
}
\item{tfs}{
a vector with TF identifiers. If NULL, the function will assess all TFs in the network.
}
\item{sampling}{
a single integer value specifying the percentage of the available samples that should be included in 
the analysis. For example, for each TF-target interaction of a given hub, 'sampling = 35' means that 
the conditional mutual information will be computed from the top and bottom 35\% of the samples ranked 
by the gene expression of a given candidate modulator.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be used 
(see 'p.adjust' for details).
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon. Gene sets with fewer than this number are 
removed from the analysis.
}
\item{minIntersectSize}{
a single integer or numeric value specifying the minimum number of observed modulated 
elements in a regulon (as percentage value).
}
\item{miThreshold}{
a single character value specifying the underlying distribution used to estimate the mutual information threshold. Options: 'md', 'md.tf' and 'md.tf.tar'. In the 1st case, 'miThreshold' is estimated from a pooled null distribution representing random modulators. In the 2nd case, specific mutual information thresholds are estimated for each TF conditioned on random modulators, while in the last case the 'miThreshold' is estimated for each individual TF-target interaction conditioned on random modulators. In all three options the 'miThreshold' is estimated by permutation analysis (see 'prob'). Alternatively, users can provide a custom mutual information threshold (numeric value >= 0).
}
\item{prob}{
a probability value in [0,1] used to estimate the 'miThreshold' based on the underlying quantile distribution.
}
\item{pwtransform}{
a single logical value specifying whether a Box-Cox power transformation should be used to improve the data symmetry (this procedure might be carefully used as a complementary test to check the validity of the inferred associations).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}

}

\value{
a data frame in the slot "results", see 'cdt' option in \code{\link[RTN:tni.get]{tni.get}}.
}

\author{Mauro Castro}

\references{
Wang, K. et al. \emph{Genome-wide identification of post-translational modulators of transcription factor activity in human B cells}. Nat Biotechnol, 27(9):829-39, 2009.

Castro, M.A.A. et al. \emph{RTN: Reconstruction and Analysis of Transcriptional Networks}. Journal Paper (in preparation), 2012.
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(dt4rtn)

# a few TFs for quick demonstration!
tfs4test<-dt4rtn$tfs[c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")]

# create a new TNI object
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=tfs4test)

\dontrun{

# preprocessing
rtni <- tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)

# permutation analysis (infers the reference/relevance network)
rtni<-tni.permutation(rtni)

# dpi filter (infers the transcriptional network)
rtni<-tni.dpi.filter(rtni)

# ..and a few candidate modulators for demonstration!
mod4test<-rownames(rtni@gexp)[sample(1:nrow(rtni@gexp),200)]

# conditional analysis
rtni<-tni.conditional(rtni, modulators=mod4test, pValueCutoff=1e-3)

#get results
cdt<-tni.get(rtni,what="cdt")

#get summary on a graph object
g<-tni.graph(rtni,gtype="mmap")

###---------------------------------------------
### optional: plot the igraph object using RedeR
library(RedeR)

#--load reder interface
rdp<-RedPort()
calld(rdp)

#---add graph and legends
addGraph(rdp,g)
addLegend.shape(rdp,g)
addLegend.size(rdp,g)
addLegend.color(rdp,g,type="edge")
relax(rdp,p1=50,p5=20)

}
}

\keyword{methods}

