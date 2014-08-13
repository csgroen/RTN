################################################################################
##########################         TNI Class        ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "TNI",
          function(.Object, gexp, transcriptionFactors) {
            ##-----check arguments
            if(missing(gexp))stop("NOTE: 'gexp' is missing!",call.=FALSE)    
            if(missing(transcriptionFactors))stop("NOTE: 'transcriptionFactors' is missing!",call.=FALSE)            
            tnai.checks(name="gexp",gexp)
            tnai.checks(name="transcriptionFactors",transcriptionFactors)
            ##-----initialization
            .Object@gexp<-gexp
            .Object@transcriptionFactors<-transcriptionFactors
            .Object@modulators<-character()
            ##-----result slot
            .Object@results<-list()
            ##-----status matrix
            .Object@status <- rep("[ ]", 1, 5)
            names(.Object@status) <- c("Preprocess", "Permutation", "Bootstrap", "DPI.filter", "Conditional")
            ##-----summary info
            ##-----tfs
            sum.info.tfs<-matrix(,1,2)
            rownames(sum.info.tfs)<-"TF"
            colnames(sum.info.tfs)<-c("input","valid")           
            ##-----parameters
            sum.info.para <- list()
            sum.info.para$perm<-matrix(,1,6)
            colnames(sum.info.para$perm)<-c("pValueCutoff","pAdjustMethod", "globalAdjustment",
                                       "estimator", "nPermutations","pooledNullDistribution")
            rownames(sum.info.para$perm)<-"Parameter"
            sum.info.para$boot<-matrix(,1,3)
            colnames(sum.info.para$boot)<-c("estimator", "nBootstraps", "consensus")        
            rownames(sum.info.para$boot)<-"Parameter"
            sum.info.para$dpi<-matrix(,1,1)
            colnames(sum.info.para$dpi)<-c("eps")       
            rownames(sum.info.para$dpi)<-"Parameter"
            sum.info.para$cdt<-matrix(,1,8)
            colnames(sum.info.para$cdt)<-c("sampling","pValueCutoff","pAdjustMethod","minRegulonSize",
                                           "minIntersectSize","miThreshold","prob","pwtransform")
            rownames(sum.info.para$cdt)<-"Parameter"
            ##-----results
            sum.info.results<-list()
            sum.info.results$tnet<-matrix(,2,3)
            colnames(sum.info.results$tnet)<-c("TFs","Targets","Edges")
            rownames(sum.info.results$tnet)<-c("tnet.ref","tnet.dpi")
            .Object@summary<-list(tfs=sum.info.tfs,para=sum.info.para,results=sum.info.results)			
            .Object
          }
)
##------------------------------------------------------------------------------
##pre-processing
setMethod(
  "tni.preprocess",
  "TNI",
  function(object, gexpIDs=NULL, cvfilter=TRUE, verbose=TRUE){
    ##-----check input arguments
    gexpIDs=tnai.checks(name="gexpIDs",para=gexpIDs)
    tnai.checks(name="cvfilter",para=cvfilter)
    tnai.checks(name="verbose",para=verbose)
    ##-----preprocessing
    if(verbose)cat("-Preprocessing for input data...\n")
    ##-----check gexpIDs if available
    if(!is.null(gexpIDs)){
      if(verbose)cat("--Mapping 'gexp' to 'gexpIDs'...\n")
      if( any(!rownames(object@gexp)%in%rownames(gexpIDs)) ){
        stop("NOTE: all rownames in 'gexp' should be available in col1 of 'gexpIDs'!",call.=FALSE)
      }
      if(cvfilter){
        if(verbose)cat("--Removing duplicated genes (keep max coefficient of variation!)...\n")
        #i.e. col1=probe, col2=gene (collapse cv by col2)
        cvres<-cv.filter(object@gexp, gexpIDs)
        object@gexp<-cvres$gexp
        object@annotation<-cvres$ids
      } else {
        #or leave by the user!!!
        object@annotation<-gexpIDs[rownames(object@gexp),]
        tp1<-paste("by setting 'cvfilter=FALSE', please note that both 'gexp' and 'gexpIDs'\n")
        tp2<-paste("should be provided with unique and matched probe-to-gene identifiers!", sep="")          
        if(verbose)warning(tp1,tp2,call.=FALSE)
      }
      #correct 'symbol' column if a valid name 
      #ps. annotation is already check for duplicated col names!
      idx<-toupper(colnames(object@annotation))%in%"SYMBOL"
      if(any(idx)){
        idx<-which(idx)[1]
        colnames(object@annotation)[idx]<-"SYMBOL"
        #..remove any empty space or NA from SYMBOL!!!
        object@annotation$SYMBOL<-as.character(object@annotation$SYMBOL)
        idx<-is.na(object@annotation$SYMBOL)
        object@annotation$SYMBOL[idx]<-rownames(object@annotation)[idx]
        idx<-object@annotation$SYMBOL==""|object@annotation$SYMBOL=="NA"
        object@annotation$SYMBOL[idx]<-rownames(object@annotation)[idx]
      } else {
        tp1<-paste("NOTE: to get better gene summary across the pipelines, 'gexpIDs' annotation\n")
        tp2<-paste("should provide an extra column named SYMBOL!", sep="")       
        if(verbose)warning(tp1,tp2,call.=FALSE)
      }
    }
    #----check sd in gexp
    sd.check<-apply(object@gexp,1,sd)
    sd.check<-sd.check==0
    if(any(sd.check)){
      if(verbose)cat("--Removing inconsistent data: standard deviation is zero for", sum(sd.check),"gene(s)! \n")
      object@gexp<-object@gexp[!sd.check,]
    }
    #-----check TFs in gexp
    object@summary$tfs[,"input"]<-length(object@transcriptionFactors)
    if(verbose) cat("--Checking TFs in the dataset...\n")
    idxtfs<-object@transcriptionFactors%in%rownames(object@gexp)
    object@transcriptionFactors<-object@transcriptionFactors[idxtfs]
    object@summary$tfs[,"valid"]<-length(object@transcriptionFactors)
    if(length(object@transcriptionFactors)==0)stop("NOTE: input 'transcriptionFactors' contains no useful data!\n",call.=FALSE)
    ##-----make sure 'transcriptionFactors' is a named character vector
    if(!is.character(object@transcriptionFactors)){
      nm<-names(object@transcriptionFactors)
      object@transcriptionFactors<-as.character(object@transcriptionFactors)
      names(object@transcriptionFactors)<-nm
    }
    if(is.null(names(object@transcriptionFactors))){
      #..if null names, add available ones
      if(!is.null(object@annotation$SYMBOL)){
        names(object@transcriptionFactors)<-object@annotation[object@transcriptionFactors,"SYMBOL"]
      } else {
        names(object@transcriptionFactors)<-object@transcriptionFactors
      }
    } else {
      #..else remove any empty space or NA from TF names!!!
      tfnames<-names(object@transcriptionFactors)
      idx<-tfnames==""|tfnames=="NA"
      names(object@transcriptionFactors)[idx]<-object@transcriptionFactors[idx]
      #..and check possible incositency between 'gexpIDs' and 'transcriptionFactors' names
      if(!is.null(object@annotation$SYMBOL)){
        tp<-object@annotation[object@transcriptionFactors,"SYMBOL"]
        if(any(tp!=names(object@transcriptionFactors))){
          tp1<-"NOTE: inconsistent symbol(s) found in the named vector 'transcriptionFactors'!\n"
          tp2<-"Please, use symbols consistent with col <SYMBOL> in 'gexpIDs'!"
          warning(tp1,tp2)
        }
      }
    }
    ##-----updade status and return
    object@status["Preprocess"] <- "[x]"
    object@status["Permutation"] <- "[ ]"
    object@status["Bootstrap"] <- "[ ]"
    object@status["DPI.filter"] <- "[ ]"
    if(verbose)cat("-Preprocessing complete!\n\n")
    return(object)
  }
)
##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.permutation",
  "TNI",
  function(object, pValueCutoff=0.01, pAdjustMethod="BH", globalAdjustment=TRUE, estimator="pearson",
           nPermutations=1000, pooledNullDistribution=TRUE, parChunks=50, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="globalAdjustment",para=globalAdjustment)
    tnai.checks(name="estimator",para=estimator)  
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="pooledNullDistribution",para=pooledNullDistribution)
    tnai.checks(name="parChunks",para=parChunks)
    tnai.checks(name="verbose",para=verbose)
    object@para$perm<-list(pValueCutoff=pValueCutoff,pAdjustMethod=pAdjustMethod,globalAdjustment=globalAdjustment,
                      estimator=estimator,nPermutations=nPermutations,pooledNullDistribution=pooledNullDistribution)
    object@summary$para$perm[1,]<-unlist(object@para$perm)
    ###compute reference network###
    ##---permutation analysis
    if(object@para$perm$pooledNullDistribution){
      res<-tni.perm.pooled(object, parChunks, verbose)
    } else {
      res<-tni.perm(object,verbose)
    }
    #object@results$adjpv<-res$adjpv
    object@results$tn.ref<- res$tn.ref * tni.cor(object@gexp,res$tn.ref)
    object@status["Permutation"] <- "[x]"
    if(verbose)cat("-Permutation analysis complete! \n\n")
    ##update summary and return results
    bin<-object@results$tn.ref
    bin[bin!=0]<-1
    object@summary$results$tnet[1,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    return(object)
  }
)

##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.bootstrap",
  "TNI",
  function(object, estimator="pearson", nBootstraps=100, consensus=95, 
           parChunks=10, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing and permutation analysis!")
    if(object@status["Permutation"]!="[x]")stop("NOTE: input data need permutation analysis!")
    ##-----check and assign parameters
    tnai.checks(name="estimator",para=estimator)  
    tnai.checks(name="nBootstraps",para=nBootstraps)    
    tnai.checks(name="consensus",para=consensus)
    tnai.checks(name="parChunks",para=parChunks)
    tnai.checks(name="verbose",para=verbose)
    object@para$boot<-list(estimator=estimator,nBootstraps=nBootstraps,consensus=consensus)
    object@summary$para$boot[1,]<-unlist(object@para$boot)
    ##---bootstrap analysis
    object@results$tn.ref<-tni.boot(object,parChunks,verbose)
    object@status["Bootstrap"] <- "[x]"
    if(verbose)cat("-Bootstrap analysis complete! \n\n")
    ##update summary and return results
    bin<-object@results$tn.ref
    bin[bin!=0]<-1
    object@summary$results$tnet[1,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    return(object)
  }
)

##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.dpi.filter",
  "TNI",
  function(object, eps=0, verbose=TRUE){
    if(object@status["Permutation"]!="[x]")stop("NOTE: input data need permutation/bootstrep analysis!")
    ##-----check and assign parameters
    tnai.checks(name="eps",para=eps)
    tnai.checks(name="verbose",para=verbose)
    object@para$dpi<-list(eps=eps)
    object@summary$para$dpi[1,]<-unlist(object@para$dpi)
    ##---apply dpi filter
    if(verbose)cat("-Applying dpi filter...\n")
    object@results$tn.dpi<-tni.dpi(abs(object@results$tn.ref), eps=object@para$dpi$eps)
    object@results$tn.dpi<-object@results$tn.dpi * tni.cor(object@gexp,object@results$tn.dpi)
    if(verbose)cat("-DPI filter complete! \n\n")
    object@status["DPI.filter"] <- "[x]"
    ##update summary and return results
    bin<-object@results$tn.dpi
    bin[bin!=0]<-1
    object@summary$results$tnet[2,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    return(object)
  }
)


##------------------------------------------------------------------------------
##pre-processing
setMethod(
  "tni2tna.preprocess",
  "TNI",
  function(object, phenotype=NULL, hits=NULL, phenoIDs=NULL, duplicateRemoverMethod="max", verbose=TRUE) {
    if(object@status["Preprocess"]!="[x]")stop("NOTE: TNI object is not compleate: requires preprocessing!")
    if(object@status["Permutation"]!="[x]")stop("NOTE: TNI object is not compleate: requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")stop("NOTE: TNI object is not compleate: requires DPI filter!")
    ##-----check input arguments
    tnai.checks(name="TNI",para=object)
    tnai.checks(name="phenotype",para=phenotype)
    tnai.checks(name="hits",para=hits)
    phenoIDs<-tnai.checks(name="phenoIDs",para=phenoIDs)
    tnai.checks(name="duplicateRemoverMethod",para=duplicateRemoverMethod)
    tnai.checks(name="verbose",para=verbose)
    ##-----generate a new object of class TNA
    .object <- new("TNA",
                   referenceNetwork=object@results$tn.ref,
                   transcriptionalNetwork=object@results$tn.dpi, 
                   transcriptionFactors=object@transcriptionFactors, 
                   phenotype=phenotype,
                   hits=hits)
    if(nrow(object@annotation)>0).object@annotation<-object@annotation
    if(!is.null(object@results$conditional) && length(object@results$conditional)>0){
      cdt<-tni.get(object,what="cdt")
      lmod<-lapply(cdt,function(reg){
        if(nrow(reg)>0){
          tp<-reg$Mode
          names(tp)<-rownames(reg)
        } else {
          tp=character()
        }
        tp
      })
      .object@listOfModulators<-lmod
    }
    .object <- tna.preprocess(.object,
                              phenoIDs=phenoIDs,
                              duplicateRemoverMethod=duplicateRemoverMethod,
                              verbose=verbose)
    return(.object)
  }
)
##------------------------------------------------------------------------------
##get slots from TNI 
setMethod(
  "tni.get",
  "TNI",
  function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, idkey=NULL) {
    ##-----reset compatibility with old args
    if(what=="tn.dpi")what="tnet"
    if(what=="tn.ref")what="refnet"
    ##-----check input arguments
    tnai.checks(name="tni.what",para=what)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="idkey",para=idkey)
    tnai.checks(name="reportNames",para=reportNames)
    ##-----get query
    query<-NULL
    if(what=="gexp"){
      query<-object@gexp
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"matrixAndNames",reportNames)
    } else if(what=="tfs"){
      query<-object@transcriptionFactors
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,"vecAndContent",reportNames)
    } else if(what=="para"){
      query<-object@para
    } else if(what=="refnet"){
      query<-object@results$tn.ref
      if(is.null(query))stop("NOTE: empty slot!",call.=FALSE)
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,"matrixAndNames",reportNames)
    } else if(what=="tnet"){
      query<-object@results$tn.dpi
      if(is.null(query))stop("NOTE: empty slot!",call.=FALSE)
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,"matrixAndNames",reportNames)
    } else if(what=="refregulons" || what=="refregulons.and.mode"){
      query<-list()
      for(i in object@transcriptionFactors){
        idx<-object@results$tn.ref[,i]!=0
        query[[i]]<-rownames(object@results$tn.ref)[idx]
      }
      if(what=="refregulons.and.mode"){
        for(i in names(query)){
          tp<-object@results$tn.ref[query[[i]],i]
          names(tp)<-query[[i]]
          query[[i]]<-tp
        }
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
      } else {
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,"listAndContent",reportNames)
      }
    } else if(what=="regulons" || what=="regulons.and.mode"){
      query<-list()
      for(i in object@transcriptionFactors){
        idx<-object@results$tn.dpi[,i]!=0
        query[[i]]<-rownames(object@results$tn.dpi)[idx]
      }
      if(what=="regulons.and.mode"){
        for(i in names(query)){
          tp<-object@results$tn.dpi[query[[i]],i]
          names(tp)<-query[[i]]
          query[[i]]<-tp
        }
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
      } else {
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,"listAndContent",reportNames)
      }
    } else if(what=="cdt" || what=="cdtrev"){
      query<-object@results$conditional$count
      for(nm in names(query)){
        qry<-query[[nm]]
        qry<-qry[qry[,2 ]=="FALSE" | qry[,2 ]=="0",-2,drop=FALSE]
        qry<-qry[qry[,2 ]=="FALSE" | qry[,2 ]=="0",-2,drop=FALSE]
        query[[nm]]<-qry
      }
      #obs. daqui resultados saem ordenados, segundo stat disponivel!
      if(what=="cdtrev"){
        query<-cdt.getReverse(query,object@para$cdt$pAdjustMethod)
      } else {
        query<-cdt.get(query,object@para$cdt$pAdjustMethod)
      }
      if(is.null(ntop)){
        for(nm in names(query)){
          qry<-query[[nm]]
          b1 <- qry[,"AdjPvFET"] <= object@para$cdt$pValueCutoff
          b2 <- qry[,"AdjPvKS"]  <= object@para$cdt$pValueCutoff
          if(!is.null(qry$AdjPvSNR)){
            b3 <- qry[,"AdjPvSNR"] <= object@para$cdt$pValueCutoff
            qry<-qry[b1 & b2 & b3,,drop=FALSE]
          } else {
            qry<-qry[b1 & b2,,drop=FALSE]
          }
          query[[nm]]<-qry
        }
      } else {
        for(nm in names(query)){
          qry<-query[[nm]]
          qryntop<-ntop
          if(nrow(qry)>1){
            if(qryntop>nrow(qry) || qryntop<0)qryntop=nrow(qry)
            qry<-qry[1:qryntop,,drop=FALSE]
            query[[nm]]<-qry
          }
        }
      }
      query<-query[unlist(lapply(query,nrow))>0]
      if(reportNames){
        for(nm in names(query)){
          if(nrow(query[[nm]])>0){
            idx<-match(query[[nm]][,"Modulator"],object@modulators)
            query[[nm]][,"Modulator"]<-names(object@modulators)[idx]
            idx<-match(query[[nm]][,"TF"],object@transcriptionFactors)
            query[[nm]][,"TF"]<-names(object@transcriptionFactors)[idx]
          }
        }
      }
      #sort 1st list hierarchy by nrow
      if(length(query)>0)query<-sortblock.cdt(query)
      if(!is.null(idkey))warning("'idkey' argument has no effect on consolidated tables!")
    } else if(what=="summary"){
      query<-object@summary
    } else if(what=="status"){
      query<-object@status
    } else if(what=="gsea2"){
      getqs<-function(query,order=TRUE,reportNames=TRUE,ntop=NULL){
        if(is.data.frame(query) && nrow(query)>0 ){
          if(is.null(ntop)){
            query<-query[query[,"Adjusted.Pvalue"] <= object@para$gsea2$pValueCutoff,,drop=FALSE]
          } else {
            if(ntop>nrow(query)|| ntop<0)ntop=nrow(query)
            if(nrow(query)>1){
              idx<-sort.list(query[,"Pvalue"]) 
              query<-query[idx[1:ntop],,drop=FALSE]
            }
          }
          if(order){
            if(nrow(query)>1) query<-query[order(query[,"Observed.Score"]),,drop=FALSE]
          }
          if(reportNames){
            idx<-match(query[,1],object@transcriptionFactors)
            query[,1]<-names(object@transcriptionFactors)[idx]
          }
        }
        query
      }
      query<-list()
      if(is.null(ntop)){
        tp<-rownames(getqs(object@results$GSEA2.results$differential))
        tp<-intersect(tp,rownames(getqs(object@results$GSEA2.results$positive)))
        tp<-intersect(tp,rownames(getqs(object@results$GSEA2.results$negative)))
        dft<-getqs(object@results$GSEA2.results$differential,order,reportNames)
        dft<-dft[rownames(dft)%in%tp,]
        query$differential<-dft
        query$positive<-object@results$GSEA2.results$positive[rownames(dft),,drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[rownames(dft),,drop=FALSE]
      } else {
        query$differential<-getqs(object@results$GSEA2.results$differential,order,reportNames,ntop)
        query$positive<-object@results$GSEA2.results$positive[rownames(query$differential),,drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[rownames(query$differential),,drop=FALSE]
      }
    }
    
    return(query)
  }
)

##------------------------------------------------------------------------------
##get graph from TNI
## experimental args:
## mask: a logical value specifying to apply a mask on the 'amapFilter', keeping at least the 
## ......best weighted edge (when verbose=TRUE) or not (when verbose=FALSE).
## hcl: an hclust object with TF's IDs
## overlap: overlapping nodes used for the Jaccard (options: 'all', 'pos', 'neg')

setMethod(
  "tni.graph",
  "TNI",
  function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL, 
           ntop=NULL, mask=FALSE, hcl=NULL, overlap="all", xlim=c(30,80,5), breaks=NULL, nquant=NULL){
    # chech igraph compatibility
    b1<-"package:igraph0" %in% search()
    b2<- "igraph0" %in%  loadedNamespaces()
    if( b1 || b2) {
      stop("\n\n ...conflict with 'igraph0': please use the new 'igraph' package!")
    }
    ##-----check input arguments
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status["DPI.filter"]!="[x]")stop("NOTE: input data need dpi analysis!")
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="tni.gtype",para=gtype)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="amapFilter",para=amapFilter)
    tnai.checks(name="amapCutoff",para=amapCutoff)
    tnai.checks(name="mask",para=mask)
    if(!is.null(hcl))gtype="amapDend"
    if(gtype=="mmap" || gtype=="mmapDetailed")tnet="dpi"
    if(tnet=="ref"){
      tnet<-object@results$tn.ref
    } else {
      tnet<-object@results$tn.dpi
    }
    if(is.null(tfs)){
      tfs<-object@transcriptionFactors
      minsz<-colnames(tnet)[colSums(tnet!=0)>=minRegulonSize]
      tfs<-tfs[tfs%in%minsz]
    } else {
      tfs<-as.character(tfs)
      idx<-which(names(object@transcriptionFactors)%in%tfs | object@transcriptionFactors%in%tfs)
      if(length(idx)==0)stop("NOTE: input 'tfs' contains no useful data!\n")
      tfs<-object@transcriptionFactors[idx]
    }
    if(!is.null(hcl)){
      if(!all(hcl$labels%in%tfs))
        stop("all labels in the 'hclust' object should be listed as 'transcriptionFactors'!")
      tfs<-tfs[tfs%in%hcl$labels]
    }
    
    #-----------------------------------------
    #-----------------------------------------
    
    if(gtype=="mmap" || gtype=="mmapDetailed"){ #get modulatory maps
      
      ##-----check input arguments
      if(object@status["Conditional"]!="[x]")stop("NOTE: input need conditional analysis!")
      #get tfs and modulators
      cdt<-tni.get(object,what="cdt")
      if(length(cdt)==0)stop("NOTE: input conditional analysis is empty")
      #cdt<-tni.get(object,what="cdt",ntop=5)
      testedtfs<-names(cdt)
      testedtfs<-object@transcriptionFactors[object@transcriptionFactors%in%testedtfs]
      testedtfs<-testedtfs[testedtfs%in%tfs]
      if(length(testedtfs)==0)stop("NOTE: input 'tfs' contains no useful data!\n")
      modulators<-sapply(testedtfs,function(tf){
        rownames(cdt[[tf]])
      })
      modulators<-unlist(modulators)
      modulators<-object@modulators[object@modulators%in%modulators]
      othertfs<-object@transcriptionFactors
      othertfs<-othertfs[!othertfs%in%testedtfs]
      othertfs<-othertfs[othertfs%in%modulators]
      #get adjmt
      tnet<-tnet[unique(c(testedtfs,setdiff(modulators,testedtfs))),testedtfs,drop=FALSE]
      mnet<-tnet;mnet[,]=0
      junk<-sapply(colnames(mnet),function(i){
        tp<-cdt[[i]]
        mnet[rownames(tp),i]<<-tp$Mode
        NULL
      })
      pvnet<-tnet;pvnet[,]=1
      junk<-sapply(colnames(mnet),function(i){
        tp<-cdt[[i]]
        pvnet[rownames(tp),i]<<-tp$PvKS
        NULL
      })
      #---
      if(gtype=="mmapDetailed"){
        #---experimental!!!
        #return a lista with:
        #1st level: a TF
        #2nd level: all MDs of a TF
        #3rd level: a graph
        g<-tni.mmap.detailed(object,mnet,testedtfs,ntop=ntop)
      } else {
        #get mmap
        #tnet[,]<-0
        g<-tni.mmap(object,mnet,tnet,pvnet,othertfs,testedtfs,modulators)
      }
      return(g)
      
    } else if(gtype=="rmap"){
      
      tnet<-tnet[,tfs,drop=FALSE]
      g<-tni.rmap(tnet)
      #add annotation
      if(nrow(object@annotation)>0)g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      #set target names if available
      if(!is.null(V(g)$SYMBOL)){
        g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
      } else {
        V(g)$nodeAlias<-V(g)$name
      }
      #set TF names
      V(g)$tfs<-as.numeric(V(g)$name%in%tfs)
      idx<-match(tfs,V(g)$name)
      V(g)$nodeAlias[idx]<-names(tfs)
      V(g)$nodeColor<-"black"
      V(g)$nodeLineColor<-"black"
      g<-att.setv(g=g, from="tfs", to='nodeShape',title="")
      g$legNodeShape$legend<-c("nTF","TF")
      g<-att.setv(g=g, from="tfs", to='nodeSize', xlim=c(20,50,1))
      g<-att.setv(g=g, from="tfs", to='nodeFontSize',xlim=c(10,32,1))
      #remove non-usefull legends
      g<-remove.graph.attribute(g,"legNodeSize")
      g<-remove.graph.attribute(g,"legNodeFontSize")
      if(ecount(g)>0){
        #set edge attr
        g<-att.sete(g=g, from="modeOfAction", to='edgeColor',cols=c("#96D1FF","grey80","#FF8E91"), 
                    title="ModeOfAction",categvec=-1:1)
        g$legEdgeColor$legend<-c("Down","NA","Up")
        E(g)$edgeWidth<-1.5
        #map modeOfAction to node attribute (compute the average of the interactions)
        el<-data.frame(get.edgelist(g),E(g)$modeOfAction,stringsAsFactors=FALSE)
        nid<-V(g)$name
        mdmode<-sapply(nid,function(id){
          idx<-el[,2]==id
          median(el[idx,3])
        })
        mdmode[V(g)$tfs==1]=NA
        V(g)$medianModeOfAction<-as.integer(mdmode)
        #assign mode to targets
        g<-att.setv(g=g, from="medianModeOfAction", to='nodeColor',cols=c("#96D1FF","grey80","#FF8E91"), 
                    title="ModeOfAction",categvec=-1:1,pal=1,na.col="grey80")
        V(g)$nodeLineColor<-V(g)$nodeColor
        g<-remove.graph.attribute(g,"legNodeColor")
      }
      return(g)
      
    } else if(gtype=="amap"){
      
      tnet<-tnet[,tfs,drop=FALSE]
      adjmt<-tni.amap(tnet,overlap)
      #-------------------filter J.C.
      if(mask){
        #set a mask to keep at least the best weighted edge
        mask<-sapply(1:ncol(adjmt),function(i){
          tp<-adjmt[,i]
          tp==max(tp)
        })
        nc<-ncol(mask);nr<-nrow(mask)
        mask<-mask+mask[rev(nr:1),rev(nc:1)]>0
      } else {
        mask<-array(0,dim=dim(adjmt))
      }
      if(amapFilter=="phyper"){
        #filter based phyper distribution (remove non-significant overlaps)
        if(is.null(amapCutoff))amapCutoff=0.01
        pvalue<-amapCutoff
        pmat<-tni.phyper(tnet)
        adjmt[pmat>pvalue & mask==0]=0
      } else if(amapFilter=="quantile"){
        #filter based on quantile distribution
        if(is.null(amapCutoff))amapCutoff=0.75
        jc<-as.integer(amapCutoff*100)+1
        tp<-as.numeric(adjmt)
        jc<-quantile(tp[tp>0],probs = seq(0, 1, 0.01), na.rm=TRUE)[jc]
        adjmt[adjmt<jc & mask==0]=0
      } else {
        #custom filter
        if(is.null(amapCutoff))amapCutoff=0
        adjmt[adjmt<amapCutoff & mask==0]=0
      }
      #-------------------
      g<-igraph::graph.adjacency(adjmt, diag=FALSE, mode="undirected", weighted=TRUE)
      if(nrow(object@annotation)>0)g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      sz<-apply(tnet!=0, 2, sum)
      idx<-match(V(g)$name,tfs)
      V(g)$nodeAlias<-names(tfs)[idx]
      V(g)$degree<-sz[idx]
      #---set main attribs
      if(ecount(g)>0)g<-att.sete(g=g, from="weight", to='edgeWidth', nquant=5, xlim=c(1,15,1),roundleg=2)
      g<-att.setv(g=g, from="degree", to='nodeSize', xlim=c(20,100,1), nquant=5, roundleg=1,title="Regulon size")
      V(g)$nodeFontSize<-20
      return(g)
      
    } else if(gtype=="amapDend"){
      
      if(!is.null(hcl)){
        gg<-hclust2igraph(hcl)
      } else {
        x<-tni.amap(tnet[,tfs], overlap)
        diag(x)=1
        hcl <- hclust(as.dist(1-cor(x)), method='complete')
        gg<-hclust2igraph(hcl)
      }
      gg$hcl<-hcl
      sz<-apply(tnet!=0, 2, sum)
      idx<-match(V(gg$g)$name,tfs)
      V(gg$g)$nodeAlias<-names(tfs)[idx]
      V(gg$g)$nodeAlias[is.na(idx)]<-"$hcnode"
      idx<-match(V(gg$g)$name,names(sz))
      V(gg$g)$degree<-sz[idx]
      #---set main attribs
      gg$g<-att.setv(g=gg$g, from="degree", to='nodeSize', xlim=xlim, breaks=breaks, 
                     nquant=nquant, roundleg=0, title="Regulon size")
      V(gg$g)$nodeFontSize<-20
      V(gg$g)$nodeFontSize[V(gg$g)$nodeAlias=="$hcnode"]<-1
      V(gg$g)$nodeColor<-"black"
      V(gg$g)$nodeLineColor<-"black"
      E(gg$g)$edgeColor<-"black"
      return(gg)
    }
    
  }
  
)
##------------------------------------------------------------------------------
##run conditional mutual information analysis
setMethod(
  "tni.conditional",
  "TNI",
  function(object, modulators=NULL, tfs=NULL, sampling=35, pValueCutoff=0.01, 
           pAdjustMethod="bonferroni", minRegulonSize=15, minIntersectSize=5, 
           miThreshold="md", prob=0.99, pwtransform=FALSE, medianEffect=FALSE, 
           verbose=TRUE, mdStability=FALSE){
    ##-----check input arguments
    if(object@status["DPI.filter"]!="[x]")stop("NOTE: input data need dpi analysis!")
    tnai.checks(name="modulators",para=modulators)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="sampling",para=sampling)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="minIntersectSize",para=minIntersectSize)
    tnai.checks(name="miThreshold",para=miThreshold)
    tnai.checks(name="pwtransform",para=pwtransform)
    tnai.checks(name="medianEffect",para=medianEffect)
    tnai.checks(name="prob",para=prob)
    tnai.checks(name="verbose",para=verbose)
    #check additional args (experimental)
    if(is.logical(mdStability)){
      mrkboot<-NULL
    } else {
      mrkboot<-tnai.checks(name="mdStability.custom",para=mdStability)
      mdStability<-TRUE
    }
    ##-----par info
    object@para$cdt<-list(sampling=sampling, pValueCutoff=pValueCutoff,
                          pAdjustMethod=pAdjustMethod, minRegulonSize=minRegulonSize, 
                          minIntersectSize=minIntersectSize, miThreshold=NA, prob=prob,
                          pwtransform=pwtransform)
    ##-----summary info
    cdt<-unlist(object@para$cdt)
    object@summary$para$cdt<-matrix(cdt,nrow=1,ncol=8)
    rownames(object@summary$para$cdt)<-"Parameter"
    colnames(object@summary$para$cdt)<-names(cdt)

    if(verbose)cat("-Preprocessing for input data...\n")
    ##-----make sure all 'tfs' are valid
    if(is.null(tfs)){
      tfs<-object@transcriptionFactors
    } else {
      if(verbose)cat("--Checking TFs in the dataset...\n")
      tfs<-as.character(tfs)
      idx<-which(!tfs%in%object@transcriptionFactors & !tfs%in%names(object@transcriptionFactors))
      if(length(idx)>0){
        message(paste("Note: input 'tfs' contains", length(idx)," element(s) not listed in the network!\n"))
      }
      idx<-which(object@transcriptionFactors%in%tfs | names(object@transcriptionFactors)%in%tfs)
      if(length(idx)<1){
        stop(paste("NOTE: input 'tfs' contains no useful data!\n"))
      }      
      tfs<-object@transcriptionFactors[idx]
    }
    ##-----make sure 'modulators' are set to character
    if(is.null(modulators)){
      modulators<-object@transcriptionFactors
    } else {
      if(verbose)cat("--Checking modulators in the dataset...\n")
      if(!is.character(modulators)){
        nm<-names(modulators)
        modulators<-as.character(modulators)
        names(modulators)<-nm
      }
      ##-----make sure all 'modulators' are valid
      idx<-which(!modulators%in%rownames(object@gexp))
      if(length(idx)>0){
        message(paste("Note: input 'modulators' contains", length(idx)," element(s) not listed in the network!\n"))
      }
      idx<-which(modulators%in%rownames(object@gexp))
      if(length(idx)<1){
        stop(paste("NOTE: input 'modulators' contains no useful data!\n"))
      }
      modulators<-modulators[idx]
      ##-----make sure 'modulators' is a named vector
      if(is.null(names(modulators))){
        if(!is.null(object@annotation$SYMBOL)){
          names(modulators)<-object@annotation[modulators,"SYMBOL"]
        } else {
          names(modulators)<-modulators
        }
      } else {
        #check possible inconsitency between 'annotation' and 'modulators' 
        if(!is.null(object@annotation$SYMBOL)){
          tp<-object@annotation[modulators,"SYMBOL"]
          if(any(tp!=names(modulators))){
            warning("one or more symbols in the named vector 'modulators' seem to differ from 'annotation' slot!")
          }
        }
      }
      if(length(modulators)==0)stop("NOTE: incorrect number of dimensions: range constraint step requires 1 or more valid modulators!")
      #final check: remove any empty space or NA from md names!!!
      mdnames<-names(modulators)
      idx<-mdnames==""|mdnames=="NA"
      names(modulators)[idx]<-modulators[idx] 
    }
    object@modulators<-modulators
    
    ##-----get TF-targets from tnet
    if(verbose)cat("--Extracting TF-targets...\n")
    tfTargets<-list()
    tfAllTargets<-list()
    for(tf in tfs){
      idx<-object@results$tn.dpi[,tf]!=0
      tfTargets[[tf]]<-rownames(object@results$tn.dpi)[idx]
      idx<-object@results$tn.ref[,tf]!=0
      tfAllTargets[[tf]]<-rownames(object@results$tn.ref)[idx]
    }
    ##-----check regulon size
    gs.size <- unlist(
      lapply(tfTargets, length)
    )
    tfs<-tfs[tfs%in%names(gs.size[gs.size>minRegulonSize])]
    tfTargets<-tfTargets[tfs]
    tfAllTargets<-tfAllTargets[tfs]
    tnetAllTargets<-rownames(object@results$tn.dpi)[rowSums(object@results$tn.dpi!=0)>0]
    ##-----Checking independence of modulators and TFs
    if(verbose)cat("--Applying modulator independence constraint...\n")
    IConstraintList<-list()
    for(tf in tfs){
      idx<-object@results$tn.ref[,tf]!=0
      IConstraintList[[tf]]<-c(tf,rownames(object@results$tn.ref)[idx])
    }
    ##-----set sub-sample idx
    spsz<-round(ncol(object@gexp)*sampling/100,0)
    idxLow<-1:spsz
    idxHigh<-(ncol(object@gexp)-spsz+1):ncol(object@gexp)
    ##-----start filtering
    if(verbose)cat("--Applying modulator range constraint...\n")
    gxtemp<-object@gexp
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)
    gxtemp<-t(apply(gxtemp[modulators,,drop=FALSE],1,sort))[,c(idxLow,idxHigh),drop=FALSE]
    if(length(modulators)==1){
      gxtemp<-rbind(gxtemp,gxtemp)
    }
    ##--run limma
    t <- factor(c(rep("low",spsz),rep("high",spsz)))
    design <- model.matrix(~0+t)
    fit <- lmFit(gxtemp,design)
    thigh=tlow=NULL
    contrasts <- makeContrasts(thigh-tlow, levels=design)
    ct.fit <- eBayes(contrasts.fit(fit, contrasts))
    res.fit<-unclass(decideTests(ct.fit, adjust.method=object@para$cdt$pAdjustMethod, p.value=object@para$cdt$pValueCutoff))
    RConstraintList<-rownames(res.fit)[res.fit<=0]
    ##--get samples (sorted index) for each pre-selected modulator
    if(verbose)cat("--Selecting subsamples...\n")
    gxtemp<-object@gexp
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)    
    idx<-t(apply(gxtemp[modulators,,drop=FALSE],1,sort.list))
    idxLow<-idx[,idxLow,drop=FALSE]
    idxHigh<-idx[,idxHigh,drop=FALSE]
    #----power transformation
    if(pwtransform){
      if(verbose)cat("--Applying power transformation...\n")
      junk<-sapply(tfs,function(tf){
        x<-gxtemp[tf,]
        if(shapiro.test(x)$p.value<0.05){
          if(any(x<=0))x<-x+1-min(x)
          l<-coef(powerTransform(x),round=TRUE)
          x<-bcPower(x,l, jacobian.adjusted=TRUE)
          gxtemp[tf,]<<-x
        }
        NULL
      }) 
    }
    
    ##-----estimate mutual information threshold
    if(is.character(miThreshold)){
      if(verbose)cat("\n")
      if(verbose)cat("-Estimating mutual information threshold...\n")
      if(miThreshold=="md.tf"){
        mimark<-miThresholdMdTf(gxtemp,tfs=tfs,nsamples=spsz,prob=prob,nPermutations=object@para$perm$nPermutations, 
                                estimator=object@para$perm$estimator,verbose=verbose)
      } else {
        mimark<-miThresholdMd(gxtemp,nsamples=spsz,prob=prob,nPermutations=object@para$perm$nPermutations, 
                              estimator=object@para$perm$estimator,verbose=verbose)
      }
    } else {
      mimark<-sort(miThreshold)
      miThreshold<-"md"
      if(length(mimark)==1){
        mimark<-abs(mimark)
        mimark<-c(-mimark,mimark)
      } else {
        if(sum(mimark>0)!=1)stop("'miThreshold' upper and lower bounds should have different signals!")
      }
      object@summary$para$cdt[,"prob"]<-"custom"
      object@para$cdt$prob<-NA
    }
    ##-----update miThreshold
    object@summary$para$cdt[,"miThreshold"]<-miThreshold
    object@para$cdt$miThreshold<-mimark
    
    ##-----set data object to save results
    reseffect<-lapply(tfTargets,function(tar){
      data.frame(targets=tar,stringsAsFactors=FALSE)
    })
    rescount<-lapply(tfTargets,function(tar){
      res<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,stringsAsFactors=FALSE)
      colnames(res)<-c("Modulator","irConstraint","nConstraint","TF","UniverseSize","EffectSize","RegulonSize","Expected",
                       "Observed","Negative","Positive","Mode","PvFET","AdjPvFET","KS","PvKS","AdjPvKS")
      res
    })
    
    ##-----start conditional mutual information analysis
    
    modregulons<-list()
    glstat<-list()
    if(verbose)cat("\n")
    if(verbose)cat("-Performing conditional mutual information analysis...\n")
    if(verbose)cat("--For", length(tfs), "tfs and" , length(modulators), "candidate modulators \n")
    if(verbose && !mdStability) pb <- txtProgressBar(style=3)
    
    for(i in 1:length(modulators)){
      
      md<-modulators[i]
      #get sample ordering
      lw<-idxLow[md,]
      hg<-idxHigh[md,]
      #compute mi on both tails
      milow<-tni.pmim(gxtemp[,lw],tfs,estimator=object@para$perm$estimator)
      mihigh<-tni.pmim(gxtemp[,hg],tfs,estimator=object@para$perm$estimator)
      #get mi delta
      miDelta<-mihigh-milow
      #identify modulations above mi threshold
      if(miThreshold=="md.tf" && length(tfs)>1){
        sigDelta<-t(apply(miDelta,1,"<",mimark[,1])) | t(apply(miDelta,1,">",mimark[,2]))
      } else {
        sigDelta<- miDelta<mimark[1] | miDelta>mimark[2]
      }
      miDelta<-miDelta/(mihigh+milow)
      
      #check modulator stability (experimental)
      #computational cost forbids default use or advanced customization
      #better only for final verification, and one modulator each time!
      if(mdStability){
        if(is.null(mrkboot)){
          if(verbose)cat("\n-Estimating stability threshold...\n")
          mrkboot<-miThresholdMd(gxtemp,nsamples=spsz,prob=c(0.05, 0.95),
                                 nPermutations=1000,estimator=object@para$perm$estimator,
                                 verbose=verbose)
        }
        if(verbose)cat("-Checking modulation stability for", names(md), "\n")
        stabt<-cdt.stability(gxtemp,object@para$perm$estimator,mrkboot,miThreshold,
                             spsz,md,tfs,sigDelta,nboot=100,consensus=75,verbose=verbose)
        sigDelta[!stabt]<-FALSE
      } else {
        if(verbose) setTxtProgressBar(pb, i/length(modulators))
      }
      
      #decision
      miDelta[!sigDelta]<-0
      
      #run main analysis
      sapply(tfs,function(tf){
        
        dtvec<-miDelta[,tf]
        tftar<-tfTargets[[tf]]
        tfalltar<-tfAllTargets[[tf]]
        #---get SZs
        #all modulated
        EffectSZ<-sum(dtvec[tnetAllTargets]!=0)
        #all tested
        UniSZ<-length(tnetAllTargets)
        #others
        RegSZ<-length(tftar)
        ExpOV<-(EffectSZ*RegSZ)/UniSZ
        ObsOV<-sum(dtvec[tftar]!=0)
        
        #---run stats
        #check exclusion list (independence/range constraint)
        irconst<-md%in%IConstraintList[[tf]] || md%in%RConstraintList
        #check minimum number of modulated targets for testing (n constraint)
        nconst<-(ObsOV/RegSZ*100)<minIntersectSize || ExpOV<1
        if(irconst || nconst ){
          ObsPos<-NA;ObsNeg<-NA;Mode<-NA;pvfet<-NA; 
          dks<-NA;pvks<-NA;dtvec[]<-0
        } else {
          #---get mode of actions
          ObsPos<-sum(dtvec[tftar]>0)
          ObsNeg<-sum(dtvec[tftar]<0)
          Mode<-if(ObsNeg>ObsPos) -1 else if(ObsNeg<ObsPos) 1 else 0
          #---set obs to predicted mode
          Obs<-if(Mode==-1) abs(ObsNeg) else if(Mode==1) ObsPos else ObsOV
          #---run fet with phyper (obs-1)
          pvfet <- phyper(Obs-1, RegSZ, UniSZ-RegSZ, EffectSZ, lower.tail=FALSE)
          #---run ks test
          #pheno
          pheno<-abs(object@results$tn.ref[,tf])
          pheno<-pheno[pheno!=0]
          #hits
          hits<-dtvec[tftar]
          hits<-hits[hits!=0]
          hits<-which(names(pheno)%in%names(hits))
          if(length(hits)>length(pheno)/2){
            dks<-1
            pvks<-0
          } else {
            kst<-ks.test(pheno[-hits],pheno[hits],alternative="greater")
            dks<-kst$statistic
            pvks<-kst$p.value
          }
          #count (+) and (-) tf-targets in the modulated set
          #expressed by the ratio of (+) or (-) targets, respectively
          if(Mode>=0){
            mdtf.tar<-dtvec[tftar]
            mdtf.tar<-names(mdtf.tar)[mdtf.tar>0]
            tf.tar<-object@results$tn.dpi[tftar,tf]
            mdtf.tar<-tf.tar[mdtf.tar]
            p1<-sum(mdtf.tar>0)/sum(tf.tar>0)
            p2<-sum(mdtf.tar<0)/sum(tf.tar<0)
            bl<-!is.nan(p1) && !is.nan(p2)
          } else {
            mdtf.tar<-dtvec[tftar]
            mdtf.tar<-names(mdtf.tar)[mdtf.tar<0]
            tf.tar<-object@results$tn.dpi[tftar,tf]
            mdtf.tar<-tf.tar[mdtf.tar]
            p1<-sum(mdtf.tar>0)/sum(tf.tar>0)
            p2<-sum(mdtf.tar<0)/sum(tf.tar<0)
            bl<-!is.nan(p1) && !is.nan(p2)
          }
        }
        
        #---add results to a list
        reseffect[[tf]][[md]]<<-dtvec[tftar]
        rescount[[tf]][md,]<<-c(NA,NA,NA,NA,UniSZ,EffectSZ,RegSZ,ExpOV,ObsOV,ObsNeg,ObsPos,Mode,pvfet,NA,dks,pvks,NA)
        rescount[[tf]][md,c(1,2,3,4)]<<-c(md,irconst,nconst,tf)
        
        #---retain modulated targets
        mdtftar<-tftar[ dtvec[tftar]!=0]
        if(length(mdtftar)>1){
          modregulons[[md]][[tf]]<<-mdtftar
        } else {
          modregulons[[md]][[tf]]<<-c(NA,NA)
          modregulons[[md]][[tf]]<<-mdtftar
        }
        NULL
        
      })
      
      #compute mi differential score for each regulon (signal-to-noise ratio)
      #this is a global stats, only used to assess the median effect 
      #in the selected regulons
      if(medianEffect){
        sig2noise<-sapply(names(modregulons[[md]]),function(tf){
          tftar<-modregulons[[md]][[tf]]
          h<-mihigh[tftar,tf]
          l<-milow[tftar,tf]
          (median(h)-median(l))/(sd(h)+sd(l)) 
        })
        sig2noise[is.na(sig2noise)]<-0
        glstat$observed[[md]]$sig2noise<-sig2noise
      }
      
    }
    
    if(verbose && !mdStability) close(pb)
    
    #set data format
    junk<-sapply(names(rescount),function(tf){
      results<-rescount[[tf]][-1,,drop=FALSE]
      if(nrow(results)>0){
        results[,"Expected"]<-round(results[,"Expected"],2)
        results[,"KS"]<-round(results[,"KS"],2)
      }
      rescount[[tf]]<<-results
      NULL
    })
    
    ##global p.adjustment
    #rescount<-p.adjust.cdt(cdt=rescount,pAdjustMethod=pAdjustMethod, p.name="PvKS",adjp.name="AdjPvKS",sort.name="PvKS",roundpv=FALSE)
    #rescount<-p.adjust.cdt(cdt=rescount,pAdjustMethod=pAdjustMethod, p.name="PvFET",adjp.name="AdjPvFET",roundpv=FALSE)
    rescount<-sortblock.cdt(cdt=rescount,coln="PvFET")
    
    ##update summary
    object@results$conditional$count<-rescount
    object@results$conditional$effect<-reseffect
    
    #compute null based on each regulon's distribution
    #this is a global stats, only used to assess the median effect on regulons
    #...not use to infer the modulated targets
    if(medianEffect){
      if(verbose)cat("\n")
      if(verbose)cat("-Checking median modulation effect...\n") 
      modulatedTFs<-tni.get(object,what="cdt",ntop=-1)
      modulatedTFs<-unlist(lapply(modulatedTFs,nrow))
      modulatedTFs<-names(modulatedTFs)[modulatedTFs>0]      
      if(length(modulatedTFs)>0){
        if(verbose)cat("--For", length(modulators), "candidate modulators \n")
        res<-checkModuationEffect(gxtemp,tfs,modregulons,modulatedTFs,glstat,spsz,
                                  minRegulonSize,pValueCutoff,
                                  nPermutations=object@para$perm$nPermutations,
                                  estimator=object@para$perm$estimator,
                                  pAdjustMethod=pAdjustMethod,
                                  count=object@results$conditional$count,
                                  verbose)
        res$md2tf$count<-p.adjust.cdt(cdt=res$md2tf$count,pAdjustMethod=pAdjustMethod,p.name="PvSNR",
                                      adjp.name="AdjPvSNR",roundpv=FALSE, global=FALSE)       
        object@results$conditional$mdeffect$md2tf$null<-res$md2tf$null
        object@results$conditional$mdeffect$md2tf$observed<-res$md2tf$observed
        object@results$conditional$mdeffect$tf2md$null<-res$tf2md$null
        object@results$conditional$mdeffect$tf2md$observed<-res$tf2md$observed       
        object@results$conditional$count<-res$md2tf$count
        if(verbose)cat("\n")
      }
  }
  object@status["Conditional"] <- "[x]"
  if(verbose)cat("-Conditional analysis complete! \n\n")
  return(object)
  }
)
#supplementary information: get simple correlation between tfs and modulator candidates
#rnet<-tni.tfmdcor(object@gexp,tfs, modulators)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "TNI",
  function(object) {
    cat("A TNI (Transcriptional Network Inference) object:\n")
    message("--status:")
    print(tni.get(object, what=c("status")), quote=FALSE)
  }
)

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------------------------TNI INTERNAL FUNCTIONS-------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
##This function returns alternative annotations for get.tni/get.tna methods
translateQuery<-function(query,idkey,object,annottype,reportNames){
  annotation<-object@annotation
  if(is.null(query))return(query)
  cnames<-colnames(annotation)
  if(!idkey%in%cnames){
    tp1<-"'NOTE: <idkey> not available! please use one of: "
    tp2<-paste(cnames,collapse=", ")
    stop(tp1,tp2,call.=FALSE)
  }
  # get TF's lab
  tfs<-object@transcriptionFactors
  if(!reportNames){
    idx<-tfs%in%rownames(annotation)
    names(tfs)[idx]<-annotation[tfs[idx],idkey]
  }
  if(annottype=="matrixAndNames"){
    idx<-colnames(query)%in%tfs
    colnames(query)[idx]<-names(tfs)[idx]
    idx<-rownames(query)%in%rownames(annotation)
    rownames(query)[idx]<-annotation[rownames(query)[idx],idkey]
  } else if(annottype=="listAndNames"){
    idx<-names(query)%in%tfs
    names(query)[idx]<-names(tfs)[idx]
    query<-lapply(query,function(qry){
      idx<-names(qry)%in%rownames(annotation)
      names(qry)[idx]<-annotation[names(qry)[idx],idkey]
      qry
    })
  } else if(annottype=="listAndContent"){
    idx<-names(query)%in%tfs
    names(query)[idx]<-names(tfs)[idx]
    query<-lapply(query,function(qry){
      nms<-names(qry)
      idx<-qry%in%rownames(annotation)
      qry[idx]<-annotation[qry[idx],idkey]
      names(qry)<-nms
      qry<-qry[!is.na(qry)]
      unique(qry)
    })
  } else if(annottype=="vecAndContent"){
    nms<-names(query)
    idx<-query%in%rownames(annotation)
    query[idx]<-annotation[query[idx],idkey]
    names(query)<-nms
    query<-unique(query[!is.na(query)])
  } else if(annottype=="vecAndNames"){
    idx<-names(query)%in%rownames(annotation)
    names(query)[idx]<-annotation[names(query)[idx],idkey]
  }
  return(query)
}


# #---chisq test
# if(Obs>ExpOV){
#   yt <- min(0.5,Obs-ExpOV)
#   chist<-((Obs-ExpOV)-yt)^2/ExpOV
# } else {
#   chist<-0
# }
# pvalue<-pchisq(chist, df=1, lower.tail=FALSE)
