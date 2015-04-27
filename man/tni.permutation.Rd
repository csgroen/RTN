\name{tni.permutation}
\alias{tni.permutation}

\title{
Inference of transcriptional networks.
}

\description{
This function takes a TNI object and returns a transcriptional network inferred by mutual information
(with multiple hypothesis testing corrections).
}

\usage{
tni.permutation(object, pValueCutoff=0.01, pAdjustMethod="BH", globalAdjustment=TRUE,
        estimator="pearson",nPermutations=1000, pooledNullDistribution=TRUE, 
        parChunks=50, verbose=TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details).
}
\item{globalAdjustment}{
a single logical value specifying to run global p.value adjustments (when 
globalAdjustment=TRUE) or not (when globalAdjustment=FALSE).
}
\item{estimator}{
a character string indicating which estimator to be used for mutual information 
computation. One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
}
\item{nPermutations}{
a single integer value specifying the number of permutations for deriving TF-target 
p-values in the mutual information analysis. If running in parallel, nPermutations 
should be greater and multiple of parChunks.
}
\item{pooledNullDistribution}{
a single logical value specifying to run the permutation analysis with pooled regulons 
(when pooledNullDistribution=TRUE) or not (when pooledNullDistribution=FALSE).
}
\item{parChunks}{
an optional single integer value specifying the number of permutation chunks to 
be used in the parallel analysis (effective only for "pooledNullDistribution = TRUE").
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE)
}

}

\value{
a mutual information matrix in the slot "results" containing a reference transcriptional network, 
see 'tn.ref' option in \code{\link[RTN:tni.get]{tni.get}}.
}

\author{Mauro Castro}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link[snow:makeCluster]{makeCluster}}
}

\examples{

data(dt4rtn)

# just a few TFs for quick demonstration!
tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")

# create a new TNI object
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[tfs4test])

\dontrun{

# preprocessing
rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)

# linear version!
rtni<-tni.permutation(rtni)

# parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(3, "SOCK"))
rtni<-tni.permutation(rtni)
stopCluster(getOption("cluster"))
}
}


\keyword{methods}


