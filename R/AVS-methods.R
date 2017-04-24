################################################################################
##########################         AVS Class        ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "AVS",
          function(.Object, markers) {
            ##-----check arguments
            if(missing(markers))stop("NOTE: 'markers' is missing!",call.=FALSE) 
            markers=avs.checks(name="markers",markers)
            ##-----initialization
            .Object@validatedMarkers<-markers
            .Object@markers<-markers$rsid
            .Object@variantSet<-list()
            .Object@randomSet<-list()
            .Object@results<-list()
            ##-----status matrix
            .Object@status <- rep("[ ]", 1, 3)
            names(.Object@status) <- c("Preprocess", "VSE", "EVSE")
            ##-----summary info
            ##-----markers
            sum.info.markers<-matrix(,1,4)
            rownames(sum.info.markers)<-"Marker"
            colnames(sum.info.markers)<-c("input","valid","universe.removed","colinked.removed")
            ##-----parameters
            sum.info.para <- list()
            sum.info.para$avs<-matrix(,1,3)
            colnames(sum.info.para$avs)<-c("nrand","reldata","snpop")
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
##pre-processing
setMethod(
  "avs.vse",
  "AVS",
  function(object, annotation, maxgap=0, pValueCutoff=0.05, boxcox=TRUE, 
           lab="annotation", glist=NULL, minSize=100, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks
    if(ncol(annotation)<3 && !is.null(glist)){
      stop("'annotation' input should also provide the IDs available the 'glist'! ",call.=FALSE)
    }
    annotation<-tnai.checks(name="annotation.vse",para=annotation)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="boxcox",para=boxcox)
    tnai.checks(name="lab",para=lab)
    glist<-tnai.checks(name="glist",para=glist)
    tnai.checks(name="minSize",para=minSize)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$vse[1,]<-c(maxgap,pValueCutoff,NA)
    object@para$vse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,pAdjustMethod="bonferroni")
    maxgap <- maxgap*1000 #set to bp
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)cat("-Checking agreement between 'glist' and 'annotation' datasets... ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        warning(tp,call.=FALSE)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        stop(tp,call.=FALSE)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
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
    
    #---map avs to annotation REVISAR esse chunk! (vide abaixo)
    annot<-getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE,getReduced=FALSE)
    annotdist<-getAnnotOverlap(vSet,annot)
    annotation$OverlapAVS<-FALSE
    idx<-match(names(annotdist),annotation$ID)
    annotation$OverlapAVS[idx]<-annotdist
    object@results$annotation$vse<-annotation
    
    #---compute enrichment stats
    object@results$stats$vse<-vseformat(object@results$vse,pValueCutoff=pValueCutoff,boxcox=boxcox)
    
    #get universe counts (marker and annotation counts)
    # REVISAR: contagem de anotacao nao relevante p/ VSE, talvez seja desncessaria
    # quando nao entrar com glist... revisar correspondente no EVSE!!!
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
  function(object, annotation, gxdata, snpdata, maxgap=250, pValueCutoff=0.05, 
           boxcox=TRUE, lab="annotation", glist=NULL, minSize=100, fineMapping=TRUE,
           verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks
    annotation<-tnai.checks(name="annotation.evse",para=annotation)
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
    maxgap <- maxgap*1000 #set to bp
    
    #---check gxdata agreement with annotation
    if(verbose)cat("-Checking agreement between 'gxdata' and 'annotation' datasets... ")
    agreement<-sum(rownames(gxdata)%in%annotation$ID)
    agreement<-agreement/nrow(gxdata)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
    if(agreement<90){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,"% of the ids in 'gxdata' are not represented in the 'annotation' dataset!",sep="")
      warning(tp,call.=FALSE)
    } else if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,"% of the ids in 'gxdata' are not represented in the 'annotation' dataset!",sep="")
      stop(tp,call.=FALSE)
    }
    annotation<-annotation[annotation$ID%in%rownames(gxdata),,drop=FALSE]
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)cat("-Checking agreement between 'glist' and 'annotation' datasets...  ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        warning(tp,call.=FALSE)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        stop(tp,call.=FALSE)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize[1]]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
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
    
    #---check avs agreement with snpdata
    if(verbose)cat("-Checking agreement between 'AVS' and 'snpdata' datasets...  ")
    rMarkers <- avs.get(object,what="randomMarkers")
    lMarkers <- avs.get(object,what="linkedMarkers")
    allMarkers <- unique(c(lMarkers,rMarkers))
    agreement<-sum(allMarkers%in%rownames(snpdata))/length(allMarkers)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n\n",sep=""))
    if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp1 <- paste("NOTE: ",idiff,"% of the SNPs in the 'AVS' are not represented in the 'snpdata'!\n",sep="")
      tp2 <- "Although the ideal case would be a perfect matching, it is common\n"
      tp3 <- "to see large GWAS studies interrogating a fraction of the annotated\n"
      tp4 <- "variation. So, given that the 'AVS' object might represent all\n"
      tp5 <- "annotated variation (e.g. the randon sets come from a SNP population),\n"
      tp6 <- "it is expected a certain level of underepresation for the 'snpdata'.\n"
      tp7 <- "Please carefully evaluate whether this number is acceptable for\n"
      tp8 <- "your study."
      warning(tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,call.=FALSE)
    }
    vSet<-object@variantSet
    rSet<-object@randomSet
    
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
        if(!bl){warning("...mismatched 'mtally' counts for ", lab,call.=FALSE)}
        evse$eqtls<-eqtls
        object@results$evse[[lab]]<<-evse
        return(NULL)
      })
      #---
      #object@results$evsemtx$probs<-getEvseMatrix(object,"probs")
      #object@results$evsemtx$fstat<-getEvseMatrix(object,"fstat")
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
  function(object, what="summary", report=FALSE, pValueCutoff=NULL){
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
      if(report){
        query<-report.vset(object@variantSet)
      } else {
        query<-object@variantSet        
      }
    } else if(what=="randomSet"){
      query<-object@randomSet
    } else if(what=="randomMarkers"){
      query<-getMarkers.rset(object@randomSet,getlinked=TRUE)
    } else if(what=="linkedMarkers"){
      query<-getMarkers.vset(object@variantSet,getlinked=TRUE)
    } else if(what=="evse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$evse$pValueCutoff
      if(!is.null(object@results$evse)){
        query<-vseformat(object@results$evse, pValueCutoff=pValueCutoff, boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
    } else if(what=="vse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$vse$pValueCutoff
      if(!is.null(object@results$vse)){
        query<-vseformat(object@results$vse, pValueCutoff=pValueCutoff, boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
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

