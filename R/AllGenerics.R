##generic functions
##-------------------------------------------------------------------------
setGeneric("tni.preprocess",
           function(object, gexpIDs=NULL, cvfilter=TRUE, verbose=TRUE)
             standardGeneric("tni.preprocess"), package="RTN")
setGeneric("tni.permutation",
           function(object, pValueCutoff=0.01, pAdjustMethod="BH", globalAdjustment=TRUE, 
                    estimator="pearson",nPermutations=1000, pooledNullDistribution=TRUE, 
                    parChunks=50, verbose=TRUE) 
             standardGeneric("tni.permutation"), package="RTN")
setGeneric("tni.bootstrap",
           function(object,estimator="pearson", nBootstraps=100, consensus=95, 
                    parChunks=10, verbose=TRUE)
             standardGeneric("tni.bootstrap"), package="RTN")
setGeneric("tni.dpi.filter",
           function(object, eps=0, verbose=TRUE)
             standardGeneric("tni.dpi.filter"), package="RTN")
setGeneric("tni.get",
           function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, idkey=NULL) 
             standardGeneric("tni.get"), package="RTN")
setGeneric("tni.graph",
           function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL, ntop=NULL, ...)
             standardGeneric("tni.graph"), package="RTN")
setGeneric("tni.conditional",
           function(object, modulators=NULL, tfs=NULL, sampling=35, pValueCutoff=0.01, 
                    pAdjustMethod="bonferroni", minRegulonSize=15, minIntersectSize=5, 
                    miThreshold="md", prob=0.99, pwtransform=FALSE, medianEffect=FALSE, 
                    verbose=TRUE, ...)
             standardGeneric("tni.conditional"), package="RTN")
setGeneric("tni2tna.preprocess",
           function(object, phenotype=NULL, hits=NULL, phenoIDs=NULL, duplicateRemoverMethod="max", 
                    verbose=TRUE) 
             standardGeneric("tni2tna.preprocess"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("tna.graph",
           function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL, ...)
             standardGeneric("tna.graph"), package="RTN")
setGeneric("tna.mra",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15,
                    tnet="dpi", verbose=TRUE) 
             standardGeneric("tna.mra"), package="RTN")
setGeneric("tna.overlap",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, 
                    tnet="ref", tfs=NULL, verbose=TRUE)
             standardGeneric("tna.overlap"), package="RTN")
setGeneric("tna.gsea1",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, 
                    nPermutations=1000,exponent=1, tnet="dpi", orderAbsValue=TRUE, stepFilter=TRUE, 
                    tfs=NULL, verbose=TRUE) 
             standardGeneric("tna.gsea1"), package="RTN")
setGeneric("tna.gsea2",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, 
                    nPermutations=1000,exponent=1, tnet="dpi", stepFilter=TRUE, tfs=NULL, verbose=TRUE) 
             standardGeneric("tna.gsea2"), package="RTN")
setGeneric("tna.synergy",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, minIntersectSize=1,
                    nPermutations=1000, exponent=1, tnet="ref", orderAbsValue=TRUE, stepFilter=TRUE, 
                    tfs=NULL, verbose=TRUE)
             standardGeneric("tna.synergy"), package="RTN")
setGeneric("tna.shadow",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, minIntersectSize=1,
                    nPermutations=1000, exponent=1, tnet="ref", orderAbsValue=TRUE, stepFilter=TRUE, 
                    tfs=NULL, verbose=TRUE)
             standardGeneric("tna.shadow"), package="RTN")
setGeneric("tna.get",
           function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, idkey=NULL) 
             standardGeneric("tna.get"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("avs.preprocess",
           function(object, nrand=1000, mergeColinked=TRUE, reldata="RTNdata.LDHapMap.rel27", snpop="all", verbose=TRUE)
             standardGeneric("avs.preprocess"), package="RTN")
setGeneric("avs.vse",
           function(object, annotation, maxgap=0, pValueCutoff=0.05, boxcox=TRUE, 
                    lab="annotation", glist=NULL, minSize=100, verbose=TRUE)
             standardGeneric("avs.vse"), package="RTN")
setGeneric("avs.evse",
           function(object, annotation, gxdata, snpdata, maxgap=250000, pValueCutoff=0.05, 
                    boxcox=TRUE, lab="annotation", glist=NULL, minSize=100, 
                    fineMapping=TRUE, verbose=TRUE)
             standardGeneric("avs.evse"), package="RTN")
setGeneric("avs.get",
           function(object, what="summary", report=FALSE, pValueCutoff=NULL) 
             standardGeneric("avs.get"), package="RTN")


