\name{RTN-package}
\alias{RTN-package}
\alias{RTN}
\docType{package}

\title{
Reconstruction and Analysis of Transcriptional Networks.
}

\description{
This package provides classes and methods for transcriptional network inference and analysis. Modulators of transcription factor activity are assessed by conditional mutual information, and master regulators are mapped to phenotypes using different strategies, e.g., gene set enrichment, shadow and synergy analyses. Additionally, master regulators can be linked to genetic markers using eQTL/VSE analysis, taking advantage of the haplotype block structure mapped to the human genome in order to explore risk-associated SNPs identified in GWAS studies.
}

\details{
\tabular{ll}{
Package: \tab RTN\cr
Type: \tab Package\cr
Depends: \tab R (>= 2.15), methods, igraph\cr
Imports: \tab RedeR, minet, snow, limma, data.table, ff, car, IRanges\cr
Suggests: \tab HTSanalyzeR, RUnit, BiocGenerics\cr
License: \tab Artistic-2.0\cr
biocViews: \tab NetworkInference, GeneRegulation, GeneExpression, GraphsAndNetworks\cr
Collate: \tab ClassUnions.R, AllChecks.R, AllClasses.R, AllGenerics.R, AllSupplements.R, AllPlots.R, TNA-methods.R, TNI-methods.R,AVS-methods.R\cr
LazyLoad: \tab yes\cr
}
}

\section{Index}{

\tabular{ll}{
\link[RTN:TNI-class]{TNI-class}: \tab an S4 class for Transcriptional Network Inference.\cr
\link[RTN:tni.preprocess]{tni.preprocess}: \tab a preprocessing method for objects of class TNI.\cr
\link[RTN:tni.permutation]{tni.permutation}: \tab inference of transcriptional networks.\cr
\link[RTN:tni.bootstrap]{tni.bootstrap}: \tab inference of transcriptional networks.\cr
\link[RTN:tni.dpi.filter]{tni.dpi.filter}: \tab data processing inequality (DPI) filter.\cr
\link[RTN:tni.conditional]{tni.conditional}: \tab conditional mutual information analysis.\cr
\link[RTN:tni.get]{tni.get}: \tab get information from individual slots in a TNI object.\cr
\link[RTN:tni.graph]{tni.graph}: \tab compute a graph from TNI objects.\cr
\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}: \tab a preprocessing method for objects of class TNI.\cr
\link[RTN:TNA-class]{TNA-class}: \tab an S4 class for Transcriptional Network Analysis.\cr
\link[RTN:tna.mra]{tna.mra}: \tab master regulator analysis (MRA) over a list of regulons.\cr
\link[RTN:tna.overlap]{tna.overlap}: \tab overlap analysis over a list of regulons.\cr
\link[RTN:tna.gsea1]{tna.gsea1}: \tab one-tailed gene set enrichment analysis (GSEA) over a list of regulons.\cr
\link[RTN:tna.gsea2]{tna.gsea2}: \tab two-tailed gene set enrichment analysis (GSEA) over a list of regulons.\cr
\link[RTN:tna.synergy]{tna.synergy}: \tab synergy analysis over a list of regulons.\cr
\link[RTN:tna.shadow]{tna.shadow}: \tab shadow analysis over a list of regulons.\cr
\link[RTN:tna.get]{tna.get}: \tab get information from individual slots in a TNA object.\cr
\link[RTN:tna.plot.gsea1]{tna.plot.gsea1}: \tab plot results from the one-tailed GSEA.\cr
\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}: \tab plot results from the two-tailed GSEA.\cr
\link[RTN:AVS-class]{AVS-class}: \tab an S4 class to do enrichment analyses in associated variant sets (AVSs).\cr
\link[RTN:avs.preprocess.LDHapMapRel27]{avs.preprocess.LDHapMapRel27}: \tab constructor for associated variant sets (AVS).\cr
\link[RTN:avs.preprocess.LD1000gRel20130502]{avs.preprocess.LD1000gRel20130502}: \tab constructor for associated variant sets (AVS).\cr
\link[RTN:avs.vse]{avs.vse}: \tab variant set enrichment analysis.\cr
\link[RTN:avs.evse]{avs.evse}: \tab an eQTL/VSE pipeline for variant set enrichment analysis.\cr
\link[RTN:avs.get]{avs.get}: \tab get information from individual slots in an AVS object.\cr
\link[RTN:avs.plot1]{avs.plot1}: \tab plot results from AVS methods, single plots.\cr
\link[RTN:avs.plot2]{avs.plot2}: \tab plot results from AVS methods, multiple plots.\cr
}
                   
Further information is available in the vignettes by typing \code{vignette("RTN")}. Documented 
topics are also available in HTML by typing \code{help.start()} and selecting the RTN package 
from the menu.

}

\author{
Maintainer: Mauro Castro <mauro.a.castro@gmail.com>
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Castro M.A.A. et al., \emph{Regulators of genetic risk of breast cancer identified by integrative network analysis.} Nature Genetics, 48:12-21, 2016.
}

\keyword{ package }
