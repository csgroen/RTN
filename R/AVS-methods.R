################################################################################
##########################         AVS Class        ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "AVS",
          function(.Object, markers) {
            ##-----check arguments
            if(missing(markers))stop("NOTE: 'markers' is missing!")            
            tnai.checks(name="markers",markers)
            ##-----initialization
            .Object@markers<-markers
            .Object@validatedMarkers<-data.frame()
            .Object@variantSet<-list()
            .Object@randomSet<-list()
            .Object@results<-list()
            ##-----status matrix
            .Object@status <- rep("[ ]", 1, 3)
            names(.Object@status) <- c("Preprocess", "VSE", "EVSE")
            ##-----summary info
            ##-----markers
            sum.info.markers<-matrix(,1,3)
            rownames(sum.info.markers)<-"Marker"
            colnames(sum.info.markers)<-c("input","valid","colinked.removed")
            ##-----parameters
            sum.info.para <- list()
            sum.info.para$avs<-matrix(,1,3)
            colnames(sum.info.para$avs)<-c("nrand","reldata","ldfilter")
            rownames(sum.info.para$avs)<-"Parameter"
            sum.info.para$vse<-matrix(,1,3)
            colnames(sum.info.para$vse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$vse)<-"Parameter"
            sum.info.para$evse<-matrix(,1,3)
            colnames(sum.info.para$evse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$evse)<-"Parameter"          
            ##-----results
            sum.info.results<-matrix(,2,1)
            colnames(sum.info.results)<-"Annotation"
            rownames(sum.info.results)<-c("VSE","EVSE")
            .Object@summary<-list(markers=sum.info.markers,para=sum.info.para,results=sum.info.results)			
            .Object
          }
)
##------------------------------------------------------------------------------
##pre-processing for an avs
setMethod(
  "avs.preprocess",
  "AVS",
  function(object, nrand=1000, mergeColinked=TRUE, reldata="rel27CEU-NCBIB36", ldfilter="DprimeLOD", 
           snpop=NULL, verbose=TRUE){
    if( ! "package:RTNdata" %in% search() ){
      stop("Please, in order to build an AVS the 'RTNdata' package should be installed!
           This package is currently available under request <mauro.a.castro at gmail.com>.")
    }
    ##-----check input arguments
    tnai.checks(name="nrand",para=nrand)
    tnai.checks(name="mergeColinked",para=mergeColinked)
    tnai.checks(name="reldata",para=reldata)
    tnai.checks(name="ldfilter",para=ldfilter)
    tnai.checks(name="snpop",para=snpop)
    tnai.checks(name="verbose",para=verbose)
    object@para$avs<-list(nrand=nrand)
    object@summary$para$avs[1,]<-c(nrand,reldata,ldfilter)
    object@summary$markers[,"input"]<-length(object@markers)
    
    #---sort markers by chrom position
    object@validatedMarkers<-getmarkers(object@markers)
    object@summary$markers[,"valid"]<-nrow(object@validatedMarkers)
    if(nrow(object@validatedMarkers)<3)stop("not enough valid rs# markers!")
    object@validatedMarkers<-sortPosition(object@validatedMarkers)
    
    #---build AVS
    variantSet<-buildAVS(object@validatedMarkers, reldata=reldata, ldfilter=ldfilter, 
                         verbose=verbose)
    if(verbose)cat("\n")
    
    #---check co-linked AVS
    if(mergeColinked){
      if(verbose)cat("Checking co-linked markers in the AVS...\n")
      variantSet<-check.colinked(variantSet, verbose=verbose)
      mkrs<-getMarkers.vset(variantSet,getlinked=FALSE)
      object@validatedMarkers<-object@validatedMarkers[mkrs,,drop=FALSE]
      object@summary$markers[,"colinked.removed"]<-nrow(object@validatedMarkers)
      if(verbose)cat("\n")
    }
    
    #---build random AVS
    randomSet<-buildRandomAVS(variantSet, nrand=nrand, reldata=reldata, 
                              snpop=snpop, verbose=verbose)
    
    #get IRanges
    if(verbose)cat("-Mapping AVS to range of integer values...\n")
    continued=FALSE
    object@variantSet<-getAvsRanges(variantSet, continued=continued)
    object@randomSet<-getRandomAvsRanges(randomSet, continued=continued, verbose=verbose)
    
    ##-----update status and return results
    object@status["Preprocess"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##pre-processing
setMethod(
  "avs.vse",
  "AVS",
  function(object, annotation, maxgap=0, pValueCutoff=0.05, boxcox=TRUE, 
           lab="annotation", glist=NULL, minSize=100, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!")
    
    #---initial checks
    annotation<-tnai.checks(name="annotation",para=annotation)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="boxcox",para=boxcox)
    tnai.checks(name="lab",para=lab)
    tnai.checks(name="glist",para=glist)
    tnai.checks(name="minSize",para=minSize)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$vse[1,]<-c(maxgap,pValueCutoff,NA)
    object@para$vse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,pAdjustMethod="bonferroni")
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)cat("-Checking agreement between 'glist' and 'annotation' datasets... ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        warning(tp)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        stop(tp)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!")
      }
      #map names to integer values
      annot<-data.table(aid=annotation$ID,ord=1:nrow(annotation))
      setkey(annot,'aid')
      glist<-lapply(glist,function(gl){
        annot[gl,nomatch=0]$ord
      })
    } else {
      glist<-list(annotation$ID)
      names(glist)<-lab
    }
    vSet<-object@variantSet
    rSet<-object@randomSet
    
    #---start vse analysis
    if(isParallel()){
      if(verbose)cat("-Running VSE analysis (parallel version)...\n")
      getTree=FALSE
    } else {
      if(verbose)cat("-Running VSE analysis...\n")
      getTree=TRUE
    }
    junk<-sapply(names(glist),function(lab){
      if(verbose)cat("--For ",lab,"...\n",sep="")
      annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,getTree=getTree,getReduced=TRUE)
      vse<-vsea(vSet,rSet,annot=annot)
      object@results$vse[[lab]]<<-vse
      return(NULL)
    })
    
    #---map avs to annotation
    annot<-getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE,getReduced=FALSE)
    annotdist<-getAnnotOverlap(vSet,annot)
    annotation$OverlapAVS<-FALSE
    idx<-match(names(annotdist),annotation$ID)
    annotation$OverlapAVS[idx]<-annotdist
    object@results$annotation$vse<-annotation
    
    #---compute enrichment stats
    object@results$stats$vse<-vseformat(object@results$vse,pValueCutoff=pValueCutoff,boxcox=boxcox)
    
    #get universe counts (marker and gene counts)
    universeCounts<-getUniverseCounts1(vSet,annotation,maxgap)
    object@results$counts$vse<-universeCounts
    
    ##-----update status and return results
    object@status["VSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##pre-processing
setMethod(
  "avs.evse",
  "AVS",
  function(object, annotation, gxdata, snpdata, maxgap=250000, pValueCutoff=0.05, 
           boxcox=TRUE, lab="annotation", glist=NULL, minSize=100, fineMapping=TRUE,
           verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!")
    
    #---initial checks
    annotation<-tnai.checks(name="annotation",para=annotation)
    tnai.checks(name="gxdata",para=gxdata)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="boxcox",para=boxcox)
    tnai.checks(name="lab",para=lab)
    glist<-tnai.checks(name="glist",para=glist)
    minSize=tnai.checks(name="evse.minSize",para=minSize)
    tnai.checks(name="fineMapping",para=fineMapping)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$evse[1,]<-c(maxgap,pValueCutoff,NA)
    object@para$evse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,pAdjustMethod="bonferroni")
      
    #---check annotation agreement with gxdata
    if(verbose)cat("-Checking agreement between 'annotation' and 'gxdata' datasets... ")
    agreement<-sum(annotation$ID%in%rownames(gxdata))
    agreement<-agreement/nrow(annotation)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
    if(agreement<90){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,"% of the ids in 'annotation' are not represented in the 'gxdata' dataset!",sep="")
      warning(tp)
    } else if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,"% of the ids in 'annotation' are not represented in the 'gxdata' dataset!",sep="")
      stop(tp)
    }
    annotation<-annotation[annotation$ID%in%rownames(gxdata),,drop=FALSE]
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)cat("-Checking agreement between 'glist' and 'annotation' datasets...  ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        warning(tp)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        stop(tp)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize[1]]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!")
      }
      ##if not fine mapping, get a proxy for the nulls with pre-predefined sizes
      if(!fineMapping){
        gsz<-unlist(lapply(glist,length))
        maxSize<-round(max(gsz)/minSize[2])*minSize[2]
        nproxy<-rev(round(seq(minSize[2],maxSize,by=minSize[2])))
        check<-sapply(gsz,function(gz){
          tp<-abs(nproxy-gz)
          nproxy[which(tp==min(tp))[1]]
        })
        nproxy<-nproxy[nproxy%in%unique(check)]
        names(nproxy)<-1:length(nproxy)
        nproxyids<-sapply(gsz,function(gz){
          tp<-abs(nproxy-gz)
          which(tp==min(tp))[1]
        })
        names(nproxyids)<-names(gsz)
      }
    } else {
      glist<-list(annotation$ID)
      names(glist)<-lab
    }
    
    #---check snpdata matrix
    b1<-!is.matrix(snpdata) && !inherits(snpdata, "ff")
    b2<-!is.integer(snpdata[1,])
    if( b1 && b2){
      stop("'snpdata' should be a matrix (or ff matrix) of integer values!",call.=FALSE)
    }
    b1<-is.null(colnames(snpdata)) || is.null(rownames(snpdata))
    b2<-length(unique(rownames(snpdata))) < length(rownames(snpdata))
    b3<-length(unique(colnames(snpdata))) < length(colnames(snpdata))   
    if(  b1 || b2 || b3 ){
      stop("'snpdata' matrix should be named on rows and cols (unique names)!",call.=FALSE)
    }
    
    #---check gxdata/snpdata matching
    b1<-!all(colnames(gxdata)%in%colnames(snpdata))
    b2<-ncol(gxdata)!=ncol(snpdata)
    if(b1 || b2){
      stop("inconsistent 'gxdata' and 'snpdata' colnames!",call.=FALSE)
    }
    idx<-match(colnames(snpdata),colnames(gxdata))
    gxdata<-gxdata[,idx]
    vSet<-object@variantSet
    rSet<-object@randomSet
    
    #if possible coerce IRanges to IntervalTree
    #if(!isParallel()){
    #  if(verbose)cat("-Coercing IRanges to IntervalTree...\n")
    #  vSet<-coerceAvsRanges(object@variantSet)
    #  rSet<-coerceRandomAvsRanges(object@randomSet)
    #}
    
    #---set marker ids to integer in order to improve computational performance  
    if(verbose)cat("-Mapping marker ids to 'snpdata'...\n")
    vSet<-mapvset(vSet,snpnames=rownames(snpdata))
    rSet<-maprset(rSet,snpnames=rownames(snpdata),verbose=verbose)
    cat("\n")
    
    #--- start evse analysis
    if(isParallel()){
      getTree=FALSE
      if(fineMapping){
        if(verbose)cat("-Running EVSE analysis (parallel version)...\n")
      } else {
        if(verbose)cat("-Running EVSE analysis - pooled null (parallel version)...\n")
      }
    } else {
      getTree=TRUE
      if(fineMapping){
        if(verbose)cat("-Running EVSE analysis...\n")
      } else {
        if(verbose)cat("-Running EVSE analysis - pooled null...\n")
      }
    }
    if(fineMapping){
      junk<-sapply(names(glist),function(lab){
        if(verbose)cat("--For ",lab,"...\n",sep="")
        #---run evsea
        annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,getTree=getTree, getReduced=FALSE)
        evse<-evsea(vSet,rSet,annot,gxdata,snpdata,pValueCutoff=pValueCutoff,verbose=verbose)
        #---get individual eqtls
        annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,getTree=FALSE,getReduced=FALSE)
        eqtls<-eqtlExtract(vSet,annot,gxdata,snpdata,pValueCutoff)
        #---check and save
        mtally<-names(evse$mtally[evse$mtally])
        bl<-all(unique(eqtls$RiskSNP)%in%mtally)
        if(!bl){warning("...mismatched 'mtally' counts for ", lab)}
        evse$eqtls<-eqtls
        object@results$evse[[lab]]<<-evse
        return(NULL)
      })
      #---
      #object@results$evsemtx$probs<-getEvseMatrix(object,"probs")
      #object@results$evsemtx$fstat<-getEvseMatrix(object,"fstat")
      #object@results$evsemtx$coef<-getEvseMatrix(object,"coef")
    } else {
      #---run evsea null
      nullproxy<-sapply(1:length(nproxy),function(i){
        if(verbose)cat("-- ",i,"/",length(nproxy),"...\n",sep="")
        annot<-sample(1:nrow(annotation),nproxy[i])
        annot<-getAnnotRanges(annotation[annot,],maxgap=maxgap,getTree=getTree,getReduced=FALSE)
        nullproxy<-evseaproxy(rSet,annot,gxdata,snpdata,pValueCutoff=pValueCutoff,verbose=verbose)
        return(nullproxy)
      })
      #---run evsea
      if(verbose)cat("-- concluding batch processing...\n")
      if(verbose) pb <- txtProgressBar(style=3)
      junk<-sapply(1:length(glist),function(i){
        if(verbose) setTxtProgressBar(pb, i/length(glist))
        lab<-names(glist)[i]
        #---run evsea
        annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,getTree=getTree,getReduced=FALSE)
        #compute evse
        evse<-list()
        evse$mtally<-get.eqtldist(vSet, annot, gxdata, snpdata, pValueCutoff=pValueCutoff)
        evse$nulldist<-nullproxy[,nproxyids[lab]]
        evse$nclusters<-length(evse$mtally)
        #---get individual eqtls (OBS: rever, usa versao antiga do variantSet)
        object@results$evse[[lab]]<<-evse
        return(NULL)
      })
      if(verbose) close(pb)
    }
    
    #---map avs to annotation
    annot<-getAnnotRanges(annotation,maxgap=maxgap, getTree=FALSE, getReduced=FALSE)
    annotdist<-getAnnotOverlap(vSet,annot)
    annotation$OverlapAVS<-FALSE
    idx<-match(names(annotdist),annotation$ID)
    annotation$OverlapAVS[idx]<-annotdist
    object@results$annotation$evse<-annotation
    
    #---compute enrichment stats
    object@results$stats$evse<-vseformat(object@results$evse,pValueCutoff=pValueCutoff,boxcox=boxcox)
    
    #get universe counts (marker and gene counts)
    universeCounts<-getUniverseCounts2(vSet,annotation,maxgap)
    object@results$counts$evse<-universeCounts
    
    ##-----update status and return results
    object@status["EVSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##get slots from AVS 
setMethod(
  "avs.get",
  "AVS",
  function(object, what="summary", pValueCutoff=NULL) {
    ##-----check input arguments
    tnai.checks(name="avs.what",para=what)
    if(!is.null(pValueCutoff))tnai.checks(name="pValueCutoff",para=pValueCutoff)
    ##-----get query
    query<-NULL
    if(what=="markers"){
      query<-object@markers
    } else if(what=="validatedMarkers"){
      query<-object@validatedMarkers      
    } else if(what=="variantSet"){
      query<-object@variantSet      
    } else if(what=="randomSet"){
      query<-object@randomSet
    } else if(what=="randomMarkers"){
      query<-getMarkers.rset(object@randomSet,getlinked=TRUE)
    } else if(what=="linkedMarkers"){
      query<-getMarkers.vset(object@variantSet,getlinked=TRUE)
    } else if(what=="evse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$evse$pValueCutoff
      query<-vseformat(object@results$evse, pValueCutoff=pValueCutoff, boxcox=TRUE)
    } else if(what=="vse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$vse$pValueCutoff
      query<-vseformat(object@results$vse, pValueCutoff=pValueCutoff, boxcox=TRUE)
    } else if(what=="summary"){
      query<-object@summary
    } else if(what=="status"){
      query<-object@status
    }
    return(query)
  }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "AVS",
  function(object) {
    cat("An AVS (Associated Variant Set) object:\n")
    message("--status:")
    print(avs.get(object, what=c("status")), quote=FALSE)
  }
)

