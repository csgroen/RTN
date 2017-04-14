\name{avs.evse}
\alias{avs.evse}

\title{An eQTL/VSE pipeline for variant set enrichment analysis.}

\description{
The VSE method (\code{\link[RTN:avs.vse]{avs.vse}}) provides a robust framework to cope with the heterogeneous structure of haplotype blocks, and has been designed to test enrichment in cistromes and epigenomes. In order to extend the variant set enrichment to genes this pipeline implements an additional step using expression quantitative trait loci (eQTLs).
}

\usage{
avs.evse(object, annotation, gxdata, snpdata, maxgap=250, pValueCutoff=0.05, 
boxcox=TRUE, lab="annotation", glist=NULL, minSize=100,  fineMapping=TRUE, 
verbose=TRUE)
}

\arguments{
\item{object}{
an object. When this function is implemented as the S4 method of class 
\code{\link[RTN:AVS-class]{AVS-class}}, this argument is an object of class 'AVS'.
}
\item{annotation}{
a data frame with genomic annotations listing chromosome coordinates to which a particular 
property or function has been attributed. It should include the following columns: <CHROM>, <START>, <END> and <ID>. The <ID> column can be any genomic identifier, while values in <CHROM> should be listed in ['chr1', 'chr2', 'chr3' ..., 'chrX']. Both <START> and <END> columns correspond to chromosome positions mapped to the human genome assembly used to build the AVS object (see \code{\link[RTN:avs.preprocess.LDHapMapRel27]{avs.preprocess.LDHapMapRel27}}).
}
\item{gxdata}{
object of class "matrix", a gene expression matrix.
}
\item{snpdata}{
either an object of class "matrix" or "ff", a single nucleotide polymorphism (SNP) matrix.
}
\item{maxgap}{
a single integer value specifying the max distant (kb) between the AVS and the annotation used to compute the eQTL analysis.
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
if 'glist' is provided, this argument is a single integer or numeric value specifying the minimum number of elements for each gene set in the 'glist'. Gene sets with fewer than this number are removed from the analysis. if 'fineMapping=FALSE', an alternative min size value can be provided as a vector of the form c(minSize1, minSize2) used to space the null distributions (see 'fineMapping').
}
\item{fineMapping}{
if 'glist' is provided, this argument is a single logical value specifying to compute individual null distributions, sized for each gene set (when fineMapping=TRUE). This option has a significant impact on the running time required to perform the computational analysis, especially for large gene set lists. When fineMapping=FALSE, a low resolution analysis is performed by pre-computing a fewer number of null distributions of different sizes (spaced by 'minSize'), and then used as a proxy of the nulls.
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
### Example of EVSE analysis for sets of genomic 
### annotations (e.g. regulons, gene sets, etc.)
##################################################

#--- step 1: load a precomputed AVS (same 'bcavs' object as above!)
data(bcavs, package="RTNdata.LDHapMapRel27")

#--- step 2: load genomic annotation for all genes
genemap <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg18.knownGene))
genemap <- genemap[,c("seqnames","start","end","gene_id")]
colnames(genemap) <- c("CHROM","START","END","ID")

#--- step 3: load a TNI object, or any other source of regulons (e.g. gene sets)
#--- and prepare a gene set list (gene ids should be the same as in the 'genemap' object)
data("rtni1st")
glist <- tni.get(rtni1st,what="refregulons",idkey="ENTREZ")
glist <- glist[ c("FOXA1","GATA3","ESR1") ] #reduce the list just for demonstration!

#--- step 4: input matched variation and gene expression datasets!
#--- here we use two "toy" datasets for demonstration purposes only.
data(toy_snpdata, package="RTNdata.LDHapMapRel27")
data(toy_gxdata, package="RTNdata.LDHapMapRel27")

#--- step 5: run the avs.evse pipeline
bcavs<-avs.evse(bcavs, annotation=genemap, gxdata=toy_gxdata, snpdata=toy_snpdata, 
                glist=glist, pValueCutoff=0.01)

#--- step 6: generate the EVSE plots
avs.plot2(bcavs,"evse",height=2.5)

### NOTE REGARDING THIS EXAMPLE ####
#- This example is for demonstration purposes only. Despite the toy datasets, 
#- both the AVS and regulons are derived from true observations. So, any 
#- eventual positive/negative associations derived from these datasets are 
#- not comparable with the original studies that described the method 
#- (doi: 10.1038/ng.3458; 10.1038/ncomms3464).
####################################

}
}

\keyword{methods}

