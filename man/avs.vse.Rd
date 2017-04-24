\name{avs.vse}
\alias{avs.vse}

\title{Variant set enrichment (VSE) analysis.}

\description{
The VSE method tests the enrichment of an AVS for a particular trait in a genomic annotation.
}

\usage{
avs.vse(object, annotation, maxgap=0, pValueCutoff=0.05, boxcox=TRUE, 
lab="annotation", glist=NULL, minSize=100, verbose=TRUE)
}

\arguments{
\item{object}{
an object. When this function is implemented as the S4 method of class 
\code{\link[RTN:AVS-class]{AVS-class}}, this argument is an object of class 'AVS'.
}
\item{annotation}{
a data frame with genomic annotations listing chromosome coordinates to which a particular property or function has been attributed. It should include the following columns: <CHROM>, <START>, <END> and <ID>. The <ID> column can be any genomic identifier, while values in <CHROM> should be listed in ['chr1', 'chr2', 'chr3' ..., 'chrX']. Both <START> and <END> columns correspond to chromosome positions mapped to the human genome assembly used to build the AVS object (see \code{\link[RTN:avs.preprocess.LDHapMapRel27]{avs.preprocess.LDHapMapRel27}}).
}
\item{maxgap}{
a single integer value specifying the max distant (kb) between the AVS and the annotation used to compute the enrichment analysis.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered significant.
}
\item{boxcox}{
a single logical value specifying to use Box-Cox procedure to find a transformation of the null that approaches normality (when boxcox=TRUE) or not (when boxcox=FALSE). See \code{\link[car:powerTransform]{powerTransform}} and \code{\link[car:bcPower]{bcPower}}.
}
\item{lab}{
a single character value specifying a name for the annotation dataset (this option is overrided if 'glist' is used).
}
\item{glist}{
an optional list with character vectors mapped to the 'annotation' data via <ID> column. This option can be used to run a batch mode for gene sets and regulons.
}
\item{minSize}{
if 'glist' is provided, this argument is a single integer or numeric value specifying the minimum number of elements for each gene set in the 'glist'. Gene sets with fewer than this number are removed from the analysis.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:AVS-class]{AVS-class}}
}

\examples{

\dontrun{
# This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDHapMapRel27)
library(Fletcher2013b)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

##################################################
### Build AVS and random AVSs (mapped to hg18)
##################################################

#--- step 1: load 'risk SNPs' data (e.g. BCa risk SNPs from the GWAS catalog)
data(bcarisk, package="RTNdata.LDHapMapRel27")

#--- step 2: build an AVS and 1000 matched random AVSs for the input 'risk SNPs'
bcavs <- avs.preprocess.LDHapMapRel27(bcarisk, nrand=1000)

##################################################
### Example of VSE analysis for ERa and FOXA1 
### cistromes (one genomic annotation each time)
##################################################

#--- step 1: load a precomputed AVS (same 'bcavs' object as above!)
data(bcavs, package="RTNdata.LDHapMapRel27")

#--- step 2: load cistrome data from the Fletcher2013b package
#NOTE: Fletcher2013b is a large data package, but only two 'bed files' 
#are used to illustrate this analysis (ESR1bdsites and FOXA1bdsites). 
#these bed files provide ERa and FOXA1 binding sites mapped by 
#ChIP-seq experiments
data(miscellaneous)

#--- step 3: run the avs.vse pipeline
bcavs <- avs.vse(bcavs, annotation=ESR1bdsites$bdsites, pValueCutoff=0.001, lab="ERa")
bcavs <- avs.vse(bcavs, annotation=FOXA1bdsites$bdsites, pValueCutoff=0.001, lab="FOXA1")

#--- step 4: generate the VSE plots
avs.plot2(bcavs,"vse",height=2.2)

##################################################
### Example of VSE analysis for sets of genomic 
### annotations (e.g. regulons, gene sets, etc.)
##################################################

#--- step 1: load the precomputed AVS (same 'bcavs' object as above!)
data("bcavs", package="RTNdata.LDHapMapRel27")

#--- step 2: load genomic annotation for all genes
genemap <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg18.knownGene))
genemap <- genemap[,c("seqnames","start","end","gene_id")]
colnames(genemap) <- c("CHROM","START","END","ID")

#--- step 3: load a TNI object, or any other source of regulons (e.g. gene sets)
#--- and prepare a gene set list (gene ids should be the same as in the 'genemap' object)
data("rtni1st")
glist <- tni.get(rtni1st,what="refregulons",idkey="ENTREZ")
glist <- glist[ c("FOXA1","GATA3","ESR1") ] #reduce the list just for demonstration!

#--- step 4: run the avs.vse pipeline
bcavs<-avs.vse(bcavs, annotation=genemap, glist=glist, pValueCutoff=0.05)

#--- step 5: generate the VSE plots
avs.plot2(bcavs,"vse",height=2.5)

### NOTE REGARDING THIS EXAMPLE ####
#- This example is for demonstration purposes only;
#- we recommend using the EVSE/eQTL approach when analysing genes/regulons.
#- Also, the AVS object here is not the same as the one used in the study that  
#- extended the method (doi:10.1038/ng.3458), so the results are not comparable;
#- (here fewer risk SNPs are considered, and without the eQTL step).
####################################

}
}

\keyword{methods}

