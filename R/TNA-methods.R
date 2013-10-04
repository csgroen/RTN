################################################################################
##########################         TNA Class        ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
		"TNA",
		function(.Object, transcriptionalNetwork, referenceNetwork, 
             transcriptionFactors, phenotype=NULL, hits=NULL) {
			##-----check arguments
			if(missing(transcriptionalNetwork))stop("NOTE: 'transcriptionalNetwork' is missing!")
			if(missing(referenceNetwork))referenceNetwork=transcriptionalNetwork
			if(missing(transcriptionFactors))stop("NOTE: 'transcriptionFactors' is missing!")
			tnai.checks(name="transcriptionalNetwork",transcriptionalNetwork)
			tnai.checks(name="referenceNetwork",referenceNetwork)
			tnai.checks(name="transcriptionFactors",transcriptionFactors)      
			tnai.checks(name="phenotype",phenotype)
			tnai.checks(name="hits",hits)
      if(is.null(phenotype) && is.null(hits)){
       stop("NOTE: either 'phenotype' or 'hits' should be available!")
      }
      b1<-sum(!colnames(referenceNetwork)==colnames(transcriptionalNetwork)) > 0
			b2<-sum(!rownames(referenceNetwork)==rownames(transcriptionalNetwork)) > 0
			if(b1 || b2) stop("NOTE: col and row names in 'referenceNetwork' should match 'transcriptionalNetwork'!")
			if(sum(!transcriptionFactors%in%colnames(transcriptionalNetwork))>0)
			  stop("NOTE: one or more 'transcriptionFactors' missing in the 'transcriptionalNetwork'!")
			#if(length(transcriptionFactors)<2){
			#  stop("NOTE: require at least two 'transcriptionFactors' for overlap, synergy and shadow analyses!")
			#}
      if(is.null(names(transcriptionFactors)))names(transcriptionFactors)<-transcriptionFactors      
			##-----initialization
			.Object@transcriptionalNetwork<-transcriptionalNetwork
			.Object@referenceNetwork<-referenceNetwork
			.Object@transcriptionFactors<-transcriptionFactors
			.Object@phenotype<-phenotype
			.Object@hits<-hits
			.Object@annotation<-data.frame()
			.Object@listOfRegulons<-list()
			.Object@listOfReferenceRegulons<-list()
			.Object@listOfModulators<-list()
			.Object@para<-list()
			.Object@results<-list()
			#######summary info######
			##-----tnet targets
			sum.info.tar<-matrix(,1,3)
			colnames(sum.info.tar)<-c("input","valid","duplicate.removed")
			rownames(sum.info.tar)<-"Tagets"
			##-----regulon collection
			#sum.info.rgc<-matrix(,1,4)
			#colnames(sum.info.rgc)<-c("input","valid","empty.removed","above.min.size")
			sum.info.rgc<-matrix(,1,3)
			colnames(sum.info.rgc)<-c("input","valid","above.min.size")
			rownames(sum.info.rgc)<-"Regulons"
			##-----gene hits
			sum.info.hits<-matrix(,1,3)
			colnames(sum.info.hits)<-c("input","valid","duplicate.removed")
			rownames(sum.info.hits)<-"Hits"
			##-----gene list
			sum.info.gl<-matrix(,1,3)
			colnames(sum.info.gl)<-c("input","valid","duplicate.removed")
			rownames(sum.info.gl)<-"Genes"
			##-----parameters
			sum.info.para <- list()
			sum.info.para$mra <- matrix(,1,4)
			colnames(sum.info.para$mra) <- c("pValueCutoff", "pAdjustMethod","minRegulonSize", "tnet")
			rownames(sum.info.para$mra) <- "Parameters"
			sum.info.para$overlap <- matrix(,1,4)
			colnames(sum.info.para$overlap) <- c("pValueCutoff", "pAdjustMethod", "minRegulonSize", "tnet")
			rownames(sum.info.para$overlap) <- "Parameters"
			sum.info.para$gsea1 <- matrix(,1,7)
			colnames(sum.info.para$gsea1) <- c("pValueCutoff", "pAdjustMethod", "minRegulonSize", "nPermutations", 
                                        "exponent", "tnet","orderAbsValue") 
			rownames(sum.info.para$gsea1) <- "Parameters"
			sum.info.para$synergy <- matrix(,1,8)
			colnames(sum.info.para$synergy) <- c("pValueCutoff", "pAdjustMethod", "minRegulonSize", "minIntersectSize",
                                           "nPermutations","exponent", "tnet","orderAbsValue")
			rownames(sum.info.para$synergy) <- "Parameters"
			sum.info.para$shadow <- matrix(,1,8)
			colnames(sum.info.para$shadow) <- c("pValueCutoff", "pAdjustMethod", "minRegulonSize", "minIntersectSize",
                                          "nPermutations", "exponent", "tnet","orderAbsValue")
			rownames(sum.info.para$shadow) <- "Parameters"
			sum.info.para$gsea2 <- matrix(,1,6)
			colnames(sum.info.para$gsea2) <- c("pValueCutoff", "pAdjustMethod", "minRegulonSize", "nPermutations", 
			                                  "exponent", "tnet")
			rownames(sum.info.para$gsea2) <- "Parameters"
			##-----results
			sum.info.results<-matrix(,6,1)
			colnames(sum.info.results)<-"TNA"
			rownames(sum.info.results)<-c("MRA","Overlap", "GSEA1", "Synergy","Shadow","GSEA2")
			##-----summary info--initialize
			.Object@summary<-list(tar=sum.info.tar, rgc=sum.info.rgc, hts=sum.info.hits, gl=sum.info.gl,
                            para=sum.info.para, results=sum.info.results)
			#######status info######   
			.Object@status<-list()
			.Object@status$preprocess <- rep("[ ]", 1, 3)
			names(.Object@status$preprocess) <- c("integration","phenotype", "hits")
			.Object@status$analysis <- rep("[ ]", 1, 6)
			names(.Object@status$analysis) <- c("MRA", "Overlap", "GSEA1", "Synergy","Shadow","GSEA2")
			##-----
      
			.Object
		}
)
##------------------------------------------------------------------------------
##get slots from TNA
setMethod(
  "tna.get",
  "TNA",
  function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE) {
    ##-----check input arguments
    tnai.checks(name="tna.what",para=what)
    tnai.checks(name="order",para=order)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="report",para=report)
    ##-----get query
    query<-NULL
    if(what=="tnet"){
      query<-object@transcriptionalNetwork
    } else if(what=="refnet"){
      query<-object@referenceNetwork      
    } else if(what=="tfs"){
      query<-object@transcriptionFactors
    } else if(what=="pheno"){
      query<-object@phenotype
    } else if(what=="regulons" || what=="regulons.and.pheno"){
      query<-object@listOfRegulons
      if( what=="regulons.and.pheno" && !is.null(object@phenotype) ){
        pheno<-object@phenotype
        jk<-lapply(names(query),function(rg){
          tp<-query[[rg]]
          idx<-match(tp,names(pheno))
          idx<-idx[!is.na(idx)]
          tpp<-pheno[idx]
          query[[rg]]<<-tpp
        })
      }
    } else if(what=="refregulons" || what=="refregulons.and.pheno"){
      query<-object@listOfReferenceRegulons
      if( what=="refregulons.and.pheno" && !is.null(object@phenotype) ){
        pheno<-object@phenotype
        jk<-lapply(names(query),function(rg){
          tp<-query[[rg]]
          idx<-match(tp,names(pheno))
          idx<-idx[!is.na(idx)]
          tpp<-pheno[idx]
          query[[rg]]<<-tpp
        })
      }
    } else if(what=="para"){
      query<-object@para
    } else if(what=="mra"){
      query<-object@results$MRA.results
      if(is.data.frame(query) && nrow(query)>0){
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$mra$pValueCutoff,,drop=FALSE]
        } else {
          if(ntop>nrow(query) || ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(order){
          if(nrow(query)>1) query<-query[order(query[,"Pvalue"]),,drop=FALSE]
        }
        if(reportNames){
          idx<-match(query[,1],object@transcriptionFactors)
          query[,1]<-names(object@transcriptionFactors)[idx]
        }
      }
    } 
    else if(what=="gsea1"){
      query<-object@results$GSEA1.results
      if(is.data.frame(query) && nrow(query)>0 ){
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$gsea1$pValueCutoff,,drop=FALSE]
        } else {
          if(ntop>nrow(query)|| ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(order){
          if(nrow(query)>1) query<-query[order(query[,"Pvalue"]),,drop=FALSE]
        }
        if(reportNames){
          idx<-match(query[,1],object@transcriptionFactors)
          query[,1]<-names(object@transcriptionFactors)[idx]
        }
      }
    } else if(what=="gsea2"){
      getqs<-function(query,ntop){
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
            if(nrow(query)>1) query<-query[order(query[,"Pvalue"]),,drop=FALSE]
          }
          if(reportNames){
            idx<-match(query[,1],object@transcriptionFactors)
            query[,1]<-names(object@transcriptionFactors)[idx]
          }
        }
        query
      }
      query<-list()
      query$differential<-getqs(object@results$GSEA2.results$differential,ntop=ntop)
      query$positive<-getqs(object@results$GSEA2.results$positive,ntop=nrow(query$differential))
      query$negative<-getqs(object@results$GSEA2.results$negative,ntop=nrow(query$differential))
    } else if(what=="overlap"){
      query<-object@results$overlap.results
      if(is.data.frame(query) && nrow(query)>0){
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$overlap$pValueCutoff,,drop=FALSE]
        } else {
          if(ntop>nrow(query)|| ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(order){
          if(nrow(query)>1) query<-query[order(query[,"Pvalue"],decreasing=FALSE),,drop=FALSE]
        }
        if(reportNames){
          idx<-match(query[,1],object@transcriptionFactors)
          query[,1]<-names(object@transcriptionFactors)[idx]
          idx<-match(query[,2],object@transcriptionFactors)
          query[,2]<-names(object@transcriptionFactors)[idx]        
        }
      }
    } else if(what=="synergy"){ 
      query<-object@results$synergy.results
      if(is.data.frame(query) && nrow(query)>0 ){
        idx<-!is.na(query[,"Adjusted.Pvalue"])
        query<-query[idx,]
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$synergy$pValueCutoff,,drop=FALSE]
        } else {
          if(ntop>nrow(query) || ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(order){
          if(nrow(query)>1) query<-query[order(query[,"EffectSize"],decreasing=TRUE),,drop=FALSE]
        }
        if(reportNames){
          idx<-match(query[,1],object@transcriptionFactors)
          query[,1]<-names(object@transcriptionFactors)[idx]
          idx<-match(query[,2],object@transcriptionFactors)
          query[,2]<-names(object@transcriptionFactors)[idx]        
        }
      }
    } else if(what=="shadow"){
      query<-object@results$shadow.results$results
      if(is.data.frame(query) && nrow(query)>0){
        idx<-!is.na(query[,"Adjusted.Pvalue"])
        query<-query[idx,]
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$shadow$pValueCutoff,,drop=FALSE]
        } else {
          if(ntop>nrow(query) || ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(order){
          if(nrow(query)>1) query<-query[order(query[,"EffectSize"]),,drop=FALSE]
        }
        if(reportNames){
          idx<-match(query[,1],object@transcriptionFactors)
          query[,1]<-names(object@transcriptionFactors)[idx]
          idx<-match(query[,2],object@transcriptionFactors)
          query[,2]<-names(object@transcriptionFactors)[idx]        
        }
      }
    } else if(what=="regulons.and.mode"){
      query<-list()
      for(i in names(object@listOfRegulons)){
        tp<-object@transcriptionalNetwork[object@listOfRegulons[[i]],i]
        names(tp)<-object@listOfRegulons[[i]]
        query[[i]]<-tp
      }
    } else if(what=="refregulons.and.mode"){
      query<-list()
      for(i in names(object@listOfReferenceRegulons)){
        tp<-object@referenceNetwork[object@listOfReferenceRegulons[[i]],i]
        names(tp)<-object@listOfReferenceRegulons[[i]]
        query[[i]]<-tp
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
  "TNA",
  function(object) {
    status<-tna.get(object, what=c("status"))
    cat("A TNA (transcriptional network analysis) object:\n")
    message("--preprocessing status:")
    print(status$preprocess, quote=FALSE)
    message("--analysis status:")
    print(status$analysis, quote=FALSE)
  }
)

##------------------------------------------------------------------------------
##Master regulator analysis
setMethod(
  "tna.mra",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15,
           tnet="dpi", verbose=TRUE) {
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status$preprocess["hits"]!="[x]")stop("NOTE: input 'hits' is empty and/or need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="tnet",para=tnet)    
    tnai.checks(name="verbose",para=verbose)
    object@para$mra<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                          minRegulonSize=minRegulonSize, tnet=tnet)
    object@summary$para$mra[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, tnet) 
    ##-----get tnet
    if(tnet=="dpi"){
      rgcs<-object@listOfRegulons
    } else {
      rgcs<-object@listOfReferenceRegulons
    }
    tnet.universe<-rownames(object@referenceNetwork)
    rgcs<-lapply(rgcs, intersect, y=tnet.universe)
    gs.size <- unlist(lapply(rgcs, length))
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    
    ##-----run mra analysis
    if(verbose)cat("-Performing master regulatory analysis ...\n")
    if(verbose)cat("--For", object@summary$rgc[,"above.min.size"], "regulons ...\n")
    if(object@summary$rgc[,"above.min.size"] > 0){
      MRA.results <- multiHyperGeoTest4RTN(rgcs, universe=tnet.universe, 
                                       hits=object@hits, minGeneSetSize=object@para$mra$minRegulonSize, 
                                       pAdjustMethod=object@para$mra$pAdjustMethod, verbose=verbose)
    } else {
      MRA.results <- matrix(, nrow=0, ncol=7)
    }
    colnames(MRA.results) <- c("Universe.Size", "Regulon.Size", "Total.Hits", "Expected.Hits", 
                               "Observed.Hits", "Pvalue", "Adjusted.Pvalue")
    MRA.results<-data.frame(Regulon=rownames(MRA.results),MRA.results,stringsAsFactors=FALSE)
    #-----format results 
    MRA.results$Expected.Hits<-round(MRA.results$Expected.Hits,2)
    MRA.results$Pvalue<-signif(MRA.results$Pvalue, digits=2)
    MRA.results$Adjusted.Pvalue<-signif(MRA.results$Adjusted.Pvalue, digits=2)
    ##-----add results
    object@results$MRA.results<-MRA.results
    MRA.results<-tna.get(object,what="mra", reportNames=FALSE)
    object@summary$results["MRA",]<-ifelse(is.data.frame(MRA.results),nrow(MRA.results),0)      
    if(verbose)cat("-Master regulatory analysis complete\n\n")    
    ##-----update status and return results
    object@status$analysis["MRA"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##GSEA1
setMethod(
  "tna.gsea1",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, 
           nPermutations=1000, exponent=1, tnet="dpi", orderAbsValue=TRUE, stepFilter=TRUE, 
           tfs=NULL, verbose=TRUE) {
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")stop("NOTE: input 'phenotype' is empty and/or need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="gsea.tnet",para=tnet)
    tnai.checks(name="stepFilter",para=stepFilter)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="orderAbsValue",para=orderAbsValue)
    tnai.checks(name="verbose",para=verbose)    
    object@para$gsea1<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                            minRegulonSize=minRegulonSize, nPermutations=nPermutations, 
                            exponent=exponent,tnet=tnet,orderAbsValue=orderAbsValue)
    object@summary$para$gsea1[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, 
                                     nPermutations,exponent,tnet,orderAbsValue)
    ##-----get regulons
    if(tnet=="cdt"){
      if(length(object@listOfModulators)>0){
        rgcs<-lapply(object@listOfModulators,function(reg){
          names(reg)
        })
        stepFilter=FALSE
      } else {
        stop("NOTE: slot 'listOfModulators' is emmpty!")
      }
    } else if(tnet=="ref"){
      rgcs<-object@listOfReferenceRegulons
    } else {
      rgcs<-object@listOfRegulons
    }
    rgcs<-lapply(rgcs, intersect, y=names(object@phenotype))
    gs.size <- unlist(lapply(rgcs, length))
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    ##-----stepFilter: use significant regulons inferred from MRA analysis
    if(!is.null(tfs)){
      if(sum(tfs%in%object@transcriptionFactors) > sum(tfs%in%names(object@transcriptionFactors) ) ){
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        tfs<-object@transcriptionFactors[names(object@transcriptionFactors)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid name!")
      rgcs<-rgcs[tfs]
    } else if(stepFilter){
      if(object@status$analysis["MRA"]=="[x]" && object@summary$results["MRA",]>0){
        MRAnames<-rownames(tna.get(object,what="mra", reportNames=FALSE))
        rgcs<-rgcs[MRAnames]
      } else {
        cat("obs: input using 'stepFilter'!! \n")
        if(object@status$analysis["MRA"]=="[x]"){
          cat("-MRA results without significant regulons! \n")
        } else {
          cat("-invalid MRA status for 'stepFilter' option! \n")
        }
        cat("-to run GSEA1 analysis independently from previous results try setting 'stepFilter=FALSE' ! \n\n")
        stop("NOTE: GSEA1 analysis can not be executed!")
      }
    }
    
    ##-----get ordered phenotype
    phenotype<-object@phenotype
    if(orderAbsValue)phenotype<-abs(phenotype)
    phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
    
    ##---reset names to integer values
    rgcs<-lapply(rgcs, match, table=names(phenotype))
    names(phenotype)<-1:length(phenotype)
    
    ##-----run gsea1
    GSEA1.results<-run.gsea1(
      listOfRegulons=rgcs,
      phenotype=phenotype,
      pAdjustMethod=object@para$gsea1$pAdjustMethod,
      pValueCutoff=object@para$gsea1$pValueCutoff,
      nPermutations=object@para$gsea1$nPermutations, 
      minRegulonSize=object@para$gsea1$minRegulonSize,
      exponent=object@para$gsea1$exponent,
      verbose=verbose
    )
    GSEA1.results<-data.frame(Regulon=rownames(GSEA1.results),GSEA1.results,stringsAsFactors=FALSE)
    ##-----add results
    object@results$GSEA1.results<-GSEA1.results     
    GSEA1.results<-tna.get(object,what="gsea1", reportNames=FALSE)
    object@summary$results["GSEA1",]<-ifelse(is.data.frame(GSEA1.results),nrow(GSEA1.results),0)       
    ##-----update status and return results
    object@status$analysis["GSEA1"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##GSEA2
setMethod(
  "tna.gsea2",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, 
           nPermutations=1000, exponent=1, tnet="dpi", stepFilter=TRUE, tfs=NULL, verbose=TRUE) {
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")stop("NOTE: input 'phenotype' is empty and/or need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="stepFilter",para=stepFilter)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="verbose",para=verbose)
    object@para$gsea2<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                           minRegulonSize=minRegulonSize, nPermutations=nPermutations, 
                           exponent=exponent,tnet=tnet)
    object@summary$para$gsea2[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, 
                                    nPermutations,exponent,tnet)
    ##------check phenotype for gsea2
    if(!min(object@phenotype)<0 || !max(object@phenotype)>0){
      stop("NOTE: 'phenotype' data should be provided as differential expression values (e.g. logFC)!")
    }
    ##-----get tnet and regulons
    if(tnet=="ref"){
      tnet<-object@referenceNetwork
      listOfRegulonsAndMode<-tna.get(object,what="refregulons.and.mode")
      
      refreg<-tna.get(object,what="regulons.and.mode")
      lapply(1:length(listOfRegulonsAndMode),function(i){
        tp<-setdiff(names(listOfRegulonsAndMode[[i]]),names(refreg[[i]]))
        listOfRegulonsAndMode[[i]]<<-listOfRegulonsAndMode[[i]][tp]
      })
      
    } else {
      tnet<-object@transcriptionalNetwork
      listOfRegulonsAndMode<-tna.get(object,what="regulons.and.mode")
    }
    ##-----Either tfs or stepFilter: use a sublist or significant regulons inferred from MRA analysis
    if(!is.null(tfs)){
      if(sum(tfs%in%object@transcriptionFactors) > sum(tfs%in%names(object@transcriptionFactors) ) ){
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        tfs<-object@transcriptionFactors[names(object@transcriptionFactors)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid names!")
    } else if(stepFilter){
      if(object@status$analysis["MRA"]=="[x]" && object@summary$results["MRA",]>0){
        tfs<-rownames(tna.get(object,what="mra", reportNames=FALSE))
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        cat("obs: input using 'stepFilter'!! \n")
        if(object@status$analysis["MRA"]=="[x]"){
          cat("-MRA results without significant regulons! \n")
        } else {
          cat("-invalid MRA status for 'stepFilter' option! \n")
        }
        cat("-to run GSEA2 analysis independently from previous results try setting 'stepFilter=FALSE' ! \n\n")
        stop("NOTE: GSEA2 analysis can not be executed!")
      }
    } else {
      tfs<-object@transcriptionFactors
    }
    
    ##-----check regulon size
    gs.size <- sapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      min(sum(tp>0),sum(tp<0))
    })
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    ##-----stop when no subset passes the size requirement
    if(all(gs.size<minRegulonSize)){
      tp<-" overlapped genes with the universe!\n The largest number of overlapped genes is: "
      stop(paste("NOTE: no partial regulon has minimum >= ", minRegulonSize, tp, max(gs.size), sep=""))
    }
    ##-----get filtered list
    idx <- which(gs.size >= minRegulonSize)
    gs.size <- gs.size[idx]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[idx]
    tfs<-tfs[tfs%in%names(listOfRegulonsAndMode)]
    
    ##-----get ordered phenotype
    phenotype<-object@phenotype
    phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
    
    ##---reset names to integer values
    lapply(names(listOfRegulonsAndMode), function(i){
      reg<-listOfRegulonsAndMode[[i]]
      names(listOfRegulonsAndMode[[i]])<<-match(names(reg),names(phenotype))
      NULL
    })
    names(phenotype)<-1:length(phenotype)
    
    #possible check-point for cmap!
    #names(phenotype)<-sample(names(phenotype))
    
    ##-----run gsea2
    GSEA2.results<-run.gsea2(
      listOfRegulonsAndMode=listOfRegulonsAndMode[tfs],
      phenotype=phenotype,
      pAdjustMethod=object@para$gsea2$pAdjustMethod,
      pValueCutoff=object@para$gsea2$pValueCutoff,
      nPermutations=object@para$gsea2$nPermutations, 
      exponent=object@para$gsea2$exponent,
      verbose=verbose
    )
    
    ##----run cmap on downstream tnet
    
    #get observed enrichment scores for the input regulons
    dfpheno<-GSEA2.results$differential[,"Observed.Score"]
    names(dfpheno)<-rownames(GSEA2.results$differential)
    dfpheno<-dfpheno[GSEA2.results$differential[,"Adjusted.Pvalue"]<pValueCutoff]
    
    #get the expected phenotype on the observed tnet structure
    idx<-c(names(dfpheno),setdiff(names(listOfRegulonsAndMode),names(dfpheno)))
    tnet<-tnet[idx,names(dfpheno),drop=FALSE]
    tnet<-tnet[rowSums(abs(tnet))>0,,drop=FALSE]
    expectedEffect<-tnet
    for(i in 1:nrow(expectedEffect)){
      tp<-dfpheno/abs(dfpheno);tp[is.nan(tp)]=0
      expectedEffect[i,]<-expectedEffect[i,]*tp
    }
    #get observed enrichment score for all regulons
    obscore<-run.tna.cmap(listOfRegulonsAndMode[rownames(expectedEffect)], phenotype, exponent)
    obscore<-round(obscore,4)
    #check consistency of downstream effects
    dseffect<-list()
    sapply(colnames(expectedEffect),function(reg){
      expeffect<-expectedEffect[,reg]
      expeffect<-expeffect[expeffect!=0]
      expeffect<-expeffect/abs(expeffect)
      obs<-obscore[names(expeffect)]
      obseffect<-obs/abs(obs)
      obseffect[is.nan(obseffect)]=0
      dseffect[[reg]]<<-rbind(expected.effect=expeffect,
                              observed.effect=obseffect,
                              observed.score=obs)
      NULL
    })
    GSEA2.results$dseffect<-dseffect
    
    ##-----add results
    object@results$GSEA2.results<-GSEA2.results
    GSEA2.results<-tna.get(object,what="gsea2", reportNames=FALSE)
    object@summary$results["GSEA2",]<-ifelse(is.data.frame(GSEA2.results),nrow(GSEA2.results),0)       
    ##-----update status and return results
    object@status$analysis["GSEA2"] <- "[x]"
    if(verbose)cat("-GSEA2 analysis complete \n\n")
    return(object)
  }
)
##------------------------------------------------------------------------------
##overlap
setMethod(
  "tna.overlap",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, 
           tnet="ref", tfs=NULL, verbose=TRUE) {  
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="tfs",para=tfs)
    object@para$overlap<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                              minRegulonSize=minRegulonSize, tnet=tnet) 
    object@summary$para$overlap[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, tnet)
    ##-----get regulons
    if(tnet=="dpi"){
      rgcs<-object@listOfRegulons
      tnet<-object@transcriptionalNetwork
    } else {
      rgcs<-object@listOfReferenceRegulons
      tnet<-object@referenceNetwork
    }
    if(!is.null(tfs)){
      if(sum(tfs%in%object@transcriptionFactors) > sum(tfs%in%names(object@transcriptionFactors) ) ){
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        tfs<-object@transcriptionFactors[names(object@transcriptionFactors)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid names!")
      rgcs<-rgcs[tfs]
    }
    ##-----get reg and universe size
    gs.size <- unlist(lapply(rgcs, length)) 
    tnet<-tnet[rowSums(tnet!=0)>0,]
    tnet.universe<-rownames(tnet)
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    if(object@summary$rgc[,"above.min.size"]<2){
      cat("Overlap analysis can not be executed!\n")
      if(length(object@transcriptionFactors)<2){
        cat("-obs: transcription factor slot with only one element!\n")
      }
      if(ncol(object@referenceNetwork)<2){
        cat("-obs: transcriptional network slot with only one regulon!\n")
      }
      stop("NOTE: require at least two regulons for overlap analysis!")
    }
    ##-----run overlap
    overlap.results<-run.overlap(
      listOfRegulons=rgcs,
      universe=tnet.universe,
      pAdjustMethod=object@para$overlap$pAdjustMethod,
      pValueCutoff=object@para$overlap$pValueCutoff,
      minRegulonSize=object@para$overlap$minRegulonSize,
      verbose=verbose
    )
    ##-----add results
    object@results$overlap.results<-overlap.results
    overlap.results<-tna.get(object,what="overlap", reportNames=FALSE)
    object@summary$results["Overlap",]<-ifelse(is.data.frame(overlap.results),nrow(overlap.results),0)   
    ##-----update status and return results
    object@status$analysis["Overlap"] <- "[x]"
    object@status$analysis[c("Synergy", "Shadow")] <- "[ ]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##synergy analysis
setMethod(
  "tna.synergy",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, minIntersectSize=1,
           nPermutations=1000, exponent=1, tnet="ref", orderAbsValue=TRUE, stepFilter=TRUE, tfs=NULL, verbose=TRUE) {
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")stop("NOTE: input 'phenotype' is empty and/or need preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="minIntersectSize",para=minIntersectSize)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="orderAbsValue",para=orderAbsValue)
    tnai.checks(name="stepFilter",para=stepFilter)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="verbose",para=verbose)
    object@para$synergy<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                              minRegulonSize=minRegulonSize, minIntersectSize=minIntersectSize,
                              nPermutations=nPermutations, exponent=exponent, tnet=tnet,orderAbsValue=orderAbsValue)
    object@summary$para$synergy[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, minIntersectSize, 
                                       nPermutations,exponent,tnet,orderAbsValue)
    
    ##-----get regulons
    if(tnet=="dpi"){
      rgcs<-object@listOfRegulons
    } else {
      rgcs<-object@listOfReferenceRegulons
    }
    rgcs<-lapply(rgcs, intersect, y=names(object@phenotype))
    gs.size <- unlist(lapply(rgcs, length))
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    
    ##-----stop when no gene set passes the size requirement
    if(all(gs.size==0)){
      tp<-" overlapped genes with the universe!\n The largest number of overlapped genes is: "
      stop(paste("NOTE: no regulon has >= ", minRegulonSize, tp, max(gs.size), sep=""))
    }
    
    ##-----get filtered list
    rgcs <- rgcs[which(gs.size >= minRegulonSize)]
    
    if(!is.null(tfs)){
      if(sum(tfs%in%object@transcriptionFactors) > sum(tfs%in%names(object@transcriptionFactors) ) ){
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        tfs<-object@transcriptionFactors[names(object@transcriptionFactors)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid name!")
      rgcs<-rgcs[tfs]
    } else if(stepFilter){
      if(object@status$analysis["GSEA1"]=="[x]" && object@summary$results["GSEA1",]>0){
        GSEAnames<-rownames(tna.get(object,what="gsea1", reportNames=FALSE))
        rgcs <- rgcs[which(names(rgcs)%in%GSEAnames)]
      } else {
        cat("obs: input using 'stepFilter'!! \n")
        if(object@status$analysis["GSEA1"]=="[x]"){
          cat("-GSEA1 results without significant regulons! \n")
        } else {
          cat("-invalid GSEA1 status for 'stepFilter' option! \n")
        }
        cat("-to run synergy analysis independently from previous results try setting 'stepFilter=FALSE' ! \n\n")
        stop("NOTE: synergy analysis can not be executed!")
      }
    }
    
    ##-----get regulon pairs
    collectionsOfPairsR1R2<-list()
    regpairs<-NULL
    if(length(rgcs)>2){
      for(i in 1:length(rgcs)){
        r1<-names(rgcs)[i]
          for(j in 1:length(rgcs)){
            if(j>i){
              r2<-names(rgcs)[j]
              tar.r1<-rgcs[[r1]]
              tar.r2<-rgcs[[r2]]
              collectionsOfPairsR1R2[[paste(r1,r2,sep=".")]]$Union<-union(tar.r1,tar.r2)
              collectionsOfPairsR1R2[[paste(r1,r2,sep=".")]]$Intersect<-intersect(tar.r1,tar.r2)
              regpairs<-rbind(regpairs,c(r1,r2))
            }
          }
      }
    } else {
      stop("NOTE: insufficient number of valid regulons! synergy analysis can not be executed!")
    }
    regpairs<-matrix(regpairs,ncol=2)
    colnames(regpairs)<-c("Regulon1","Regulon2")

    ##-----get ordered phenotype
    phenotype<-object@phenotype
    if(orderAbsValue)phenotype<-abs(phenotype)
    phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
    
    ##-----run tna.synergy
    if(length(collectionsOfPairsR1R2)>0){
      synergy.results<-run.synergy(
        collectionsOfPairsR1R2=collectionsOfPairsR1R2,
        labpair=regpairs,phenotype=phenotype,
        pAdjustMethod=object@para$synergy$pAdjustMethod,
        pValueCutoff=object@para$synergy$pValueCutoff,
        minIntersectSize=object@para$synergy$minIntersectSize,
        nPermutations=object@para$synergy$nPermutations,
        exponent=object@para$synergy$exponent,
        verbose=verbose
      )
    } else {
      synergy.results<-NULL
    }
    ##-----add results
    object@results$synergy.results<-synergy.results
    synergy<-tna.get(object,what="synergy", reportNames=FALSE)
    object@summary$results["Synergy",]<-ifelse(is.data.frame(synergy),nrow(synergy),0)
    ##-----update status and return results
    object@status$analysis["Synergy"] <- "[x]"
    return(object)   
  }
)
##------------------------------------------------------------------------------
##shadow analysis
setMethod(
  "tna.shadow",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, minIntersectSize=1,
           nPermutations=1000, exponent=1, tnet="ref", orderAbsValue=TRUE, stepFilter=TRUE, tfs=NULL, verbose=TRUE) {
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")stop("NOTE: input 'phenotype' is empty and/or need preprocessing!")
    if(object@status$analysis["Overlap"]!="[x]"){
      cat("-invalid 'overlap' status! \n")
      stop("NOTE: shadow requires results from 'tna.overlap' analysis! \n\n")
    }
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="minIntersectSize",para=minIntersectSize)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="orderAbsValue",para=orderAbsValue)
    tnai.checks(name="stepFilter",para=stepFilter)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="verbose",para=verbose)
    object@para$shadow<-list(pValueCutoff=pValueCutoff, pAdjustMethod=pAdjustMethod, 
                             minRegulonSize=minRegulonSize, minIntersectSize=minIntersectSize,
                             nPermutations=nPermutations, exponent=exponent, tnet=tnet,orderAbsValue=orderAbsValue)
    object@summary$para$shadow[1,]<-c(pValueCutoff, pAdjustMethod, minRegulonSize, minIntersectSize, 
                                      nPermutations,exponent,tnet,orderAbsValue)    
    ##-----overlap information is always required for shadow!
    if(object@summary$results["Overlap",]==0){
      cat("-Overlap results without significant regulon pairs! \n")
      stop("NOTE: shadow analysis can not be executed!")
    }
    ##-----get regulons
    if(tnet=="dpi"){
      rgcs<-object@listOfRegulons
    } else {
      rgcs<-object@listOfReferenceRegulons
    }
    rgcs<-lapply(rgcs, intersect, y=names(object@phenotype))
    gs.size <- unlist(lapply(rgcs, length))
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>minRegulonSize)
    ##-----stop when no gene set passes the size requirement
    if(all(gs.size==0)){
      tp<-" overlapped genes with the universe!\n The largest number of overlapped genes is: "
      stop(paste("NOTE: no regulon has >= ", minRegulonSize, tp, max(gs.size), sep=""))
    }
    ##-----get filtered list
    rgcs <- rgcs[which(gs.size >= minRegulonSize)]
    ##-----get sig. pairs from overlap analysis
    overlap.results<-tna.get(object,what="overlap", reportNames=FALSE)
    regpairs<-as.matrix(overlap.results[,c(1,2)],drop=FALSE)
    ##-----match sig. pairs with rgcs    
    idx<-regpairs[,1]%in%names(rgcs)+regpairs[,2]%in%names(rgcs)
    regpairs<-regpairs[idx==2,,drop=FALSE]
    
    if(!is.null(tfs)){
      if(sum(tfs%in%object@transcriptionFactors) > sum(tfs%in%names(object@transcriptionFactors) ) ){
        tfs<-object@transcriptionFactors[object@transcriptionFactors%in%tfs]
      } else {
        tfs<-object@transcriptionFactors[names(object@transcriptionFactors)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid name!")
      idx<-regpairs[,1]%in%tfs+regpairs[,2]%in%tfs
      regpairs<-regpairs[idx==2,,drop=FALSE]
    } else if(stepFilter){
      if(object@status$analysis["GSEA1"]=="[x]" && object@summary$results["GSEA1",]>0){
        GSEAnames<-rownames(tna.get(object,what="gsea1", reportNames=FALSE))
        idx<-regpairs[,1]%in%GSEAnames+regpairs[,2]%in%GSEAnames
        regpairs<-regpairs[idx==2,,drop=FALSE]
      } else {
        cat("obs: input using 'stepFilter'!! \n") 
        if(object@status$analysis["GSEA1"]=="[x]"){
          cat("-GSEA1 results without significant regulons! \n")
        } else {
          cat("-invalid GSEA1 status for 'stepFilter' option! \n")
        }
        cat("-to run shadow analysis independently from previous results try setting 'stepFilter=FALSE' ! \n\n")
        stop("NOTE: shadow analysis can not be executed!")
      }
    }
    ##-----get regulon pairs
    collectionsOfPairsR1<-list()
    collectionsOfPairsR2<-list()
    if(nrow(regpairs)>0){
      for(i in 1:nrow(regpairs)){
        r1<-regpairs[i,1]
        r2<-regpairs[i,2]
        tar.r1<-rgcs[[r1]]
        tar.r2<-rgcs[[r2]]
        collectionsOfPairsR1[[paste(r1,r2,sep=".")]]$R1<-tar.r1
        collectionsOfPairsR1[[paste(r1,r2,sep=".")]]$Unique.R1<-setdiff(tar.r1,tar.r2)          
        collectionsOfPairsR2[[paste(r1,r2,sep=".")]]$R2<-tar.r2
        collectionsOfPairsR2[[paste(r1,r2,sep=".")]]$Unique.R2<-setdiff(tar.r2,tar.r1)
      }
    }
    
    ##-----get ordered phenotype
    phenotype<-object@phenotype
    if(orderAbsValue)phenotype<-abs(phenotype)
    phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
    
    ##-----run shadow
    if(length(collectionsOfPairsR1)>0){
      shadow.results<-run.shadow(
        collectionsOfPairsR1=collectionsOfPairsR1,
        collectionsOfPairsR2=collectionsOfPairsR2,
        labpair=regpairs,
        phenotype=phenotype,
        pAdjustMethod=object@para$shadow$pAdjustMethod,
        pValueCutoff=object@para$shadow$pValueCutoff,
        minIntersectSize=object@para$shadow$minIntersectSize,
        nPermutations=object@para$shadow$nPermutations,
        exponent=object@para$shadow$exponent,
        verbose=verbose
      )
    } else {
      shadow.results<-NULL
    }
    ##-----add results
    object@results$shadow.results<-shadow.results
    shadow<-tna.get(object,what="shadow", reportNames=FALSE)
    object@summary$results["Shadow",]<-ifelse(is.data.frame(shadow),nrow(shadow),0)
    ##-----update status and return results
    object@status$analysis["Shadow"] <- "[x]"
    return(object)   
  }
)
##------------------------------------------------------------------------------
##get graph from TNI 
setMethod(
  "tna.graph",
  "TNA",
  function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL){
    # chech igraph compatibility
    b1<-"package:igraph0" %in% search()
    b2<- "igraph0" %in%  loadedNamespaces()
    if( b1 || b2) {
      stop("\n\n ...conflict with 'igraph0': please use the new 'igraph' package!")
    }
    if(object@status$preprocess["integration"]!="[x]")stop("NOTE: input data need preprocessing!")
    tnai.checks(name="tna.gtype",para=gtype)
    ##-----check input arguments
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="amapFilter",para=amapFilter)
    tnai.checks(name="amapCutoff",para=amapCutoff)
    ##-----get tnet
    if(tnet=="ref"){
      tnet<-object@referenceNetwork
    } else if(tnet=="dpi"){
      tnet<-object@transcriptionalNetwork
    }
    if(is.null(tfs)){
      tfs<-object@transcriptionFactors
      minsz<-colnames(tnet)[colSums(tnet!=0)>=minRegulonSize]
      tfs<-tfs[tfs%in%minsz]
    } else {
      tfs<-as.character(tfs)
      idx<-which(names(object@transcriptionFactors)%in%tfs | object@transcriptionFactors%in%tfs)
      if(length(idx)==0){
        stop("NOTE: input 'tfs' contains no useful data!\n")
      }
      tfs<-object@transcriptionFactors[idx]
    }
    tnet<-tnet[,tfs]
    if(gtype=="rmap"){ #get regulatory maps
      g<-tni.rmap(tnet)
      #add annotation
      if(.hasSlot(object, "annotation")){
        if(nrow(object@annotation)>0)g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      }
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
      g<-att.setv(g=g, from="tfs", to='nodeShape')
      g<-att.setv(g=g, from="tfs", to='nodeSize', xlim=c(20,25,1))
      g<-att.setv(g=g, from="tfs", to='nodeFontSize',xlim=c(10,25,1))
    } else { #get association maps
      adjmt<-tni.amap(tnet)
      #-------------------filter J.C.
      if(amapFilter=="phyper"){
        #filter based phyper distribution (remove non-significant overlaps)
        if(is.null(amapCutoff))amapCutoff=0.01
        pvalue<-amapCutoff
        pmat<-tni.phyper(tnet)
        adjmt[pmat>pvalue]=0
      } else if(amapFilter=="quantile"){
        #filter based on quantile distribution
        if(is.null(amapCutoff))amapCutoff=0.75
        jc<-as.integer(amapCutoff*100)+1
        tp<-as.numeric(adjmt)
        jc<-quantile(tp[tp>0],probs = seq(0, 1, 0.01), na.rm=TRUE)[jc]
        adjmt[adjmt<jc]=0
      } else {
        #custom filter
        if(is.null(amapCutoff))amapCutoff=0
        adjmt[adjmt<amapCutoff]=0
      }
      #-------------------
      g<-igraph::graph.adjacency(adjmt, diag=FALSE, mode="undirected", weighted=TRUE)
      if(.hasSlot(object, "annotation")){
        if(nrow(object@annotation)>0)g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      }
      sz<-apply(tnet!=0, 2, sum)
      idx<-match(V(g)$name,tfs)
      V(g)$nodeAlias<-names(tfs)[idx]
      V(g)$degree<-sz[idx]
      #---set main attribs
      g<-att.sete(g=g, from="weight", to='edgeWidth', nquant=10, xlim=c(1,15,1))
      g<-att.setv(g=g, from="degree", to='nodeSize', xlim=c(20,100,1), nquant=10, roundleg=1)
      V(g)$nodeFontSize<-20
    }
    return(g)
  }
)
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------------------------TNA INTERNAL FUNCTIONS-------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
##internal pre-processing (input via tni2tna.preprocess)
tna.preprocess<-function(object, phenoIDs=NULL, duplicateRemoverMethod="max", verbose=TRUE) {
  ##-----data preprocessing
  if(verbose)cat("-Preprocessing for input data ...\n")
  ##-----data preprocessing: phenotype
  if(!is.null(object@phenotype))object<-pheno.preprocess(object, phenoIDs, duplicateRemoverMethod, verbose)
  ##-----data preprocessing: hits
  if(!is.null(object@hits))object<-hits.preprocess(object, phenoIDs, verbose)
  ##-----data preprocessing: integration
  object<-data.integration(object, verbose)
  ##-----update and return preprocessing
  if(verbose)cat("-Preprocessing complete!\n\n")
  object@status$analysis[c("MRA", "Overlap", "GSEA1", "Synergy","Shadow","GSEA2")] <- "[ ]"
  return(object)
}
##This function returns a preprocessed tna object
pheno.preprocess<-function(object, phenoIDs, duplicateRemoverMethod, verbose){
  ##-----input summary         
  object@summary$gl[,"input"]<-length(object@phenotype)    
  ##-----check phenoIDs if available
  if(!is.null(phenoIDs)){
    if(verbose)cat("--Mapping 'phenotype' to 'phenoIDs' ...\n")
    ids<-phenoIDs[,2]
    names(ids)<-phenoIDs[,1]
    if(sum( !(names(object@phenotype) %in% names(ids)) )>0 ){
      stop("NOTE: all names in 'phenotype' should be available in 'phenoIDs'!")
    }
    names(object@phenotype)<-ids[names(object@phenotype)]
  }
  ##-----remove NAs in phenotype
  pheno<-object@phenotype
  idx<-!is.na(pheno) & names(pheno)!="" & !is.na(names(pheno))
  if(any(!idx)){
    if(verbose) cat("--Removing genes without names or values in 'phenotype' ...\n")
    pheno<-pheno[idx]
  }
  object@summary$gl[,"valid"]<-length(pheno) #genes with valid values
  if(length(pheno)==0)stop("NOTE: input 'phenotype' contains no useful data!\n")
  ##-----duplicate remover in phenotype
  uninames<-unique(names(pheno))
  if(length(names(pheno))>length(uninames)){
    if(verbose) cat("--Removing duplicated genes ...\n")
    pheno<-tna.duplicate.remover(phenotype=pheno,method=duplicateRemoverMethod)
  }
  object@summary$gl[,"duplicate.removed"]<-length(pheno)	#genes after removing duplicates
  if(length(pheno)==0)stop("NOTE: input 'phenotype' contains no useful data!\n")
  ##-----phenotype ordering
  #if(verbose) cat("--Ordering gene list decreasingly ...\n")
  #if(orderAbsValue)pheno<-abs(pheno)
  #object@phenotype<-pheno[order(pheno,decreasing=TRUE)]
  object@phenotype<-pheno

  ##-----check phenotype names in transcriptionalNetwork
  #if(verbose)cat("--Checking 'transcriptionalNetwork' targets in 'phenotype' ...  ")
  #idx<-rownames(object@transcriptionalNetwork) %in% names(object@phenotype)
  #checkmatch<-sum(idx)/length(idx)
  #if(checkmatch==0){
  #  stop("NOTE: 'no agreement between names in 'transcriptionalNetwork' targets and 'phenotype'!")
  #} else if(checkmatch<0.9){
  #  warning(paste("Only",round(checkmatch*100,1),"% of 'transcriptionalNetwork' targets can be mapped to 'phenotype'!"))
  #} else {
  #  if(verbose)cat(paste(round(checkmatch*100,1),"% agreement! \n"))
  #}
  
  ##-----update and return
  object@status$preprocess["phenotype"] <- "[x]"
  return(object)
}
##------------------------------------------------------------------------------
##This function returns a preprocessed tna object
hits.preprocess<-function(object, phenoIDs, verbose){
  ##-----make sure vector 'hits' is set to character
  object@hits<-as.character(object@hits)
  ##-----input summary         
  object@summary$hts[,"input"]<-length(object@hits)
  ##-----check phenoIDs if available
  if(!is.null(phenoIDs)){
    if(verbose)cat("--Mapping 'hits' to 'phenoIDs' ...\n")
    ids<-phenoIDs[,2]
    names(ids)<-phenoIDs[,1]
    if(sum( !(names(object@hits) %in% names(ids)) )>0 ){
      stop("NOTE: all names in 'hits' should be available in 'phenoIDs'!")
    }
    object@hits<-ids[object@hits]
  }
  ##-----remove duplicated hits
  if(verbose) cat("--Removing duplicated hits ...\n")
  object@hits<-unique(object@hits)
  object@hits<-object@hits[!is.na(object@hits)]
  object@summary$hts[,"duplicate.removed"]<-length(object@hits)
  if(length(object@hits)==0)stop("NOTE: input 'hits' contains no useful data!\n")
  
  ##-----remove hits not listed in the universe
  #if(verbose) cat("--Removing 'hits' not listed in 'transcriptionalNetwork' universe...\n")
  #object@hits<-intersect(object@hits,rownames(object@referenceNetwork))
  #object@summary$hts[,"valid"]<-length(object@hits)
  #if(length(object@hits)==0)stop("NOTE: input 'hits' contains no useful data!\n")
  
  ##-----update and return
  object@status$preprocess["hits"] <- "[x]"
  return(object)
}
##------------------------------------------------------------------------------
##This function returns the final preprocessed tna object
data.integration<-function(object, verbose){
  
  ##-----input summary
  object@summary$tar[,"input"]<-nrow(object@transcriptionalNetwork)
  object@summary$rgc[,"input"]<-length(object@transcriptionFactors)
  
  ##-----check annotation if available  
  if(nrow(object@annotation)>0){
    if(verbose)cat("--Mapping 'transcriptionalNetwork' annotation to 'phenotype' ...\n")
    annot<-object@annotation
    #col with possible current refids
    col0<-sapply(1:ncol(annot),function(i){
      sum(rownames(annot)%in%annot[,i],na.rm=TRUE)
    })
    if(max(col0)==nrow(annot)){
      col0 <- which(col0==max(col0))[1]
    } else {
      col0<-0
    }
    #col with possible phenoIDs (phenotype and hits)
    phenoIDs<-unique(c(names(object@phenotype),object@hits))
    col1<-sapply(1:ncol(annot),function(i){
      sum(phenoIDs%in%annot[,i],na.rm=TRUE)
    })
    col1<-which(col1==max(col1))[1]
    #col with possible gene symbols
    #..obs. posso chamar SYMBOL aqui sem problema! a entrada e' unica, via objetos TNI,
    #..ja verificados exaustivamente em pipeline anterior.
    col2<-which(colnames(annot)=="SYMBOL")[1]
    col2<-ifelse(!is.na(col2),col2,0)
    col2<-ifelse(col1!=col2,col2,0)
    #other cols
    othercols<-1:ncol(annot)
    othercols<-othercols[!othercols%in%c(col0,col1,col2)]
    if(col0!=col1){
      #set annotation to correct order
      object@annotation<-annot[,c(col1,col2,col0,othercols),drop=FALSE]
      #update rownames (na anotacao somente mais adiante por nao aceitar nomes duplicados)
      ids<-object@annotation[,1]
      names(ids)<-rownames(object@annotation)
      rownames(object@referenceNetwork)<-ids[rownames(object@referenceNetwork)]
      rownames(object@transcriptionalNetwork)<-ids[rownames(object@transcriptionalNetwork)]
      #update TFs and colnames
      tfs<-object@transcriptionFactors
      coltf<-sapply(1:ncol(object@annotation),function(i){
        sum(tfs%in%object@annotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]
      idx<-match(tfs,object@annotation[,coltf])
      tnames<-object@annotation[idx,1]
      names(tnames)<-names(tfs)
      object@transcriptionFactors<-tnames
      #update transcriptionalNetwork colnames
      tfs<-colnames(object@transcriptionalNetwork)
      coltf<-sapply(1:ncol(object@annotation),function(i){
        sum(tfs%in%object@annotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]
      idx<-match(tfs,object@annotation[,coltf])
      tnames<-object@annotation[idx,1]
      names(tnames)<-names(tfs)  
      colnames(object@transcriptionalNetwork)<-tnames  
      #update referenceNetwork colnames
      tfs<-colnames(object@referenceNetwork)
      coltf<-sapply(1:ncol(object@annotation),function(i){
        sum(tfs%in%object@annotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]      
      idx<-match(tfs,object@annotation[,coltf])
      tnames<-object@annotation[idx,1]
      names(tnames)<-names(tfs)  
      colnames(object@referenceNetwork)<-tnames
      ##-----remove unnamed nodes in the tnets
      ##tnet
      tnames<-rownames(object@transcriptionalNetwork)
      tnames<-!is.na(tnames) & tnames!=""
      object@transcriptionalNetwork<-object@transcriptionalNetwork[tnames,,drop=FALSE]
      ##refnet
      tnames<-rownames(object@referenceNetwork)
      tnames<-!is.na(tnames) & tnames!=""
      object@referenceNetwork<-object@referenceNetwork[tnames,,drop=FALSE]
      ##annotation
      tnames<-object@annotation[,1]
      tnames<-!is.na(tnames) & tnames!=""
      object@annotation<-object@annotation[tnames,,drop=FALSE]
      ##-----remove duplicate nodes in tnet
      uninames<-unique(rownames(object@transcriptionalNetwork))
      if(length(rownames(object@transcriptionalNetwork))>length(uninames)){
        if(verbose) cat("--Removing duplicated targets ...\n")
        abmax<-function(x){
          imax<-max(x);imin<-min(x)
          ifelse(imax>abs(imin),imax,imin)
        }
        tnet<-object@transcriptionalNetwork
        tnet<-aggregate(tnet,by=list(rownames(tnet)),abmax)
        rownames(tnet)<-tnet[,1]
        tnet<-as.matrix(tnet[,-1,drop=FALSE])
        tnet<-tnet[,object@transcriptionFactors,drop=FALSE]
        object@transcriptionalNetwork<-tnet    
      }
      object@summary$tar[,"duplicate.removed"]<-nrow(object@transcriptionalNetwork)
      if(prod(dim(object@transcriptionalNetwork))==0)stop("NOTE: input 'transcriptionalNetwork' contains no useful data!\n")
      ##-----duplicate remover in refnet
      uninames<-unique(rownames(object@referenceNetwork))
      if(length(rownames(object@referenceNetwork))>length(uninames)){
        abmax<-function(x){
          imax<-max(x);imin<-min(x)
          ifelse(imax>abs(imin),imax,imin)
        }
        tnet<-object@referenceNetwork
        tnet<-aggregate(tnet,by=list(rownames(tnet)),abmax)
        rownames(tnet)<-tnet[,1]
        tnet<-as.matrix(tnet[,-1,drop=FALSE])
        tnet<-tnet[,object@transcriptionFactors,drop=FALSE]
        object@referenceNetwork<-tnet    
      }
      ##-----update modulator list if available
      if(length(object@listOfModulators)>0){
        tfs<-names(object@listOfModulators)
        coltf<-sapply(1:ncol(object@annotation),function(i){
          sum(tfs%in%object@annotation[,i],na.rm=TRUE)
        })
        coltf<-which(coltf==max(coltf))[1]    
        idx<-match(tfs,object@annotation[,coltf])
        tnames<-object@annotation[idx,1] 
        names(object@listOfModulators)<-tnames
        lmod<-sapply(object@listOfModulators,function(reg){
          if(length(reg)>0){
            idx<-object@annotation[,coltf]%in%names(reg)
            mnames<-object@annotation[idx,1]
            mnames<-aggregate(reg,by=list(mnames),max)
            reg<-mnames[,2]
            names(reg)<-mnames[,1]
          }
          reg
        })
        object@listOfModulators<-lmod
      }
      ##-----duplicate remover in annotation
      #..get current refids col
      col0<-sapply(1:ncol(object@annotation),function(i){
        sum(rownames(object@annotation)%in%object@annotation[,i],na.rm=TRUE)
      })
      if(max(col0)==nrow(object@annotation)){
        col0 <- which(col0==max(col0))[1]
      } else {
        col0<-0
      }
      uninames<-unique(object@annotation[,1])
      idx<-match(uninames,object@annotation[,1])
      object@annotation<-object@annotation[idx,]
      rownames(object@annotation)<-object@annotation[,1]
      object@annotation<-object@annotation[,-col0,drop=FALSE] #agora da pra remover!
      ##-----check ordering
      tp1<-rownames(object@annotation)
      tp2<-rownames(object@transcriptionalNetwork)
      tp3<-rownames(object@referenceNetwork)
      b1<-all(tp1%in%tp2) & all(tp1%in%tp3)
      b2<-length(tp1)==length(tp2) && length(tp1)==length(tp3)
      if(b1 && b2){
        object@transcriptionalNetwork<-object@transcriptionalNetwork[rownames(object@annotation),]
        object@referenceNetwork<-object@referenceNetwork[rownames(object@annotation),]
      } else {
        warning("NOTE: possible mismatched names between 'transcriptionalNetwork' and 'annotation'!")
      }
    }
  }
  ##update
  object@summary$tar[,"valid"]<-nrow(object@transcriptionalNetwork)
  if(nrow(object@transcriptionalNetwork)==0)stop("NOTE: input 'transcriptionalNetwork' contains no useful data!\n")
  
  #-----check phenotype names in transcriptionalNetwork
  if(!is.null(object@phenotype)){
    if(verbose)cat("--Checking 'transcriptionalNetwork' in 'phenotype' ...  ")
    phenoIDs<-unique(c(names(object@phenotype),object@hits))
    idx<-rownames(object@transcriptionalNetwork) %in% phenoIDs
    checkmatch<-sum(idx)/length(idx)
    if(checkmatch==0){
      stop("NOTE: 'no agreement between 'transcriptionalNetwork' and 'phenotype' names!")
    } else if(checkmatch<0.9){
      warning(paste("Only",round(checkmatch*100,1),"% of 'transcriptionalNetwork' targets can be mapped to 'phenotype'!"))
    } else {
      if(verbose)cat(paste(round(checkmatch*100,1),"% agreement! \n"))
    }
  }
  
  #-----check and remove hits not listed in the universe
  if(!is.null(object@hits)){
    hits.int<-intersect(object@hits,rownames(object@referenceNetwork))
    if(length(hits.int)<length(object@hits)){
      if(verbose) cat("--Removing 'hits' not listed in 'transcriptionalNetwork' universe...\n")
      object@hits<-hits.int
    }
    object@summary$hts[,"valid"]<-length(object@hits)
    if(length(object@hits)==0)stop("NOTE: input 'hits' contains no useful data!\n")
  }
  
  ##-----extracting regulons from 'transcriptionalNetwork' 
  if(verbose) cat("--Extracting regulons ...\n")
  #Regulons from tnet
  listOfRegulons<-list()
  for(i in object@transcriptionFactors){
    idx<-object@transcriptionalNetwork[,i]!=0
    listOfRegulons[[i]]<-rownames(object@transcriptionalNetwork)[idx]
  }
  object@listOfRegulons<-listOfRegulons
  
  #   if(verbose) cat("--Removing empty regulons ...\n")
  #   len<-unlist(lapply(listOfRegulons,length))
  #   object@listOfRegulons<-listOfRegulons[len>0]
  #   #Regulons refnet
  #   listOfReferenceRegulons<-list()
  #   for(i in names(object@listOfRegulons)){
  #     idx<-object@referenceNetwork[,i]!=0
  #     listOfReferenceRegulons[[i]]<-rownames(object@referenceNetwork)[idx]
  #   }
  #   object@listOfReferenceRegulons<-listOfReferenceRegulons
  #   object@summary$rgc[,"empty.removed"]<-length(object@listOfRegulons)
  #   if(length(object@listOfRegulons)==0)stop("NOTE: derived 'listOfRegulons' contains no useful data!\n")
  
  #Regulons refnet
  listOfReferenceRegulons<-list()
  for(i in object@transcriptionFactors){
    idx<-object@referenceNetwork[,i]!=0
    listOfReferenceRegulons[[i]]<-rownames(object@referenceNetwork)[idx]
  }
  object@listOfReferenceRegulons<-listOfReferenceRegulons
  if(length(object@listOfRegulons)==0)stop("NOTE: derived 'listOfRegulons' contains no useful data!\n")
  
  ##-----update and return
  object@status$preprocess["integration"] <- "[x]"
  return(object)
}

##------------------------------------------------------------------------------
##This function returns the overlap among all regulons in the tnet
run.overlap <- function(listOfRegulons, universe, pAdjustMethod="BH", 
                        pValueCutoff=0.05, minRegulonSize=15, 
                        verbose=TRUE) {
  ##-----check arguments
  tnai.checks("listOfRegulons",listOfRegulons)
  tnai.checks("universe",universe)
  tnai.checks("pAdjustMethod",pAdjustMethod)
  tnai.checks("pValueCutoff",pValueCutoff)
  tnai.checks("minRegulonSize",minRegulonSize)
  tnai.checks("verbose",verbose)
  ##-----filter 'listOfRegulons' by 'minRegulonSize'
  gs.size <- unlist(
    lapply(lapply(listOfRegulons, intersect, y = universe), length)
  )
  max.size <- max(gs.size)
  ##-----stop when no gene set passes the size requirement
  if(all(unlist(lapply(listOfRegulons,length))==0)){
    tp<-" overlapped genes with the universe!\n The largest number of overlapped genes is: "
    stop(paste("NOTE: no regulon has >= ", minRegulonSize, tp, max.size, sep=""))    
  }
  ##-----get filtered list
  gs.id <- which(gs.size >= minRegulonSize)
  n.gs.discarded <- length(listOfRegulons) - length(gs.id)
  listOfRegulons <- listOfRegulons[gs.id]
  
  if(verbose)cat("-Performing overlap analysis ...\n")
  if(verbose)cat("--For", length(listOfRegulons), "regulons ...\n")
  if(length(listOfRegulons) > 0){
    HGT.results <- tna.hyper.pairs(
      listOfRegulons, universe=universe, 
      minRegulonSize = minRegulonSize, 
      pAdjustMethod = pAdjustMethod, verbose = verbose)
  } else {
    HGT.results <- matrix(, nrow=0, ncol=7)
    colnames(HGT.results) <- c("Universe.Size", "Size.Regulon1", "Size.Regulon2", "Expected.Hits", 
                              "Observed.Hits", "Pvalue", "Adjusted.Pvalue")
  }
  rownames(HGT.results)<-paste(HGT.results[[1]],HGT.results[[2]],sep=".")
  ##-----adjustment of pvalues
  HGT.results$Adjusted.Pvalue<-p.adjust(HGT.results$Adjusted.Pvalue,method=pAdjustMethod)
  #-----format results
  HGT.results[,"Expected.Overlap"]<-round(HGT.results[,"Expected.Overlap"],2)
  HGT.results[,"Pvalue"]<-signif(HGT.results[,"Pvalue"], digits=4)
  HGT.results[,"Adjusted.Pvalue"]<-signif(HGT.results[,"Adjusted.Pvalue"], digits=4)
  if(verbose)cat("-Overlap analysis complete\n\n")
  return(HGT.results)
}

##------------------------------------------------------------------------------
##This function takes a list of gene sets (regulons), a named phenotype 
##vector, and returns the results of gene set enrichment analysis for all 
##regulons (with multiple hypothesis testing correction).
run.gsea1 <- function(listOfRegulons, phenotype, pAdjustMethod="BH", 
                     pValueCutoff=0.05, nPermutations=1000, 
                     minRegulonSize=15, exponent=1, verbose=TRUE) {
  ##-----check arguments
  tnai.checks("listOfRegulons",listOfRegulons)
  tnai.checks("phenotype",phenotype)
  tnai.checks("pAdjustMethod",pAdjustMethod)
  tnai.checks("pValueCutoff",pValueCutoff)
  tnai.checks("nPermutations",nPermutations)
  tnai.checks("minRegulonSize",minRegulonSize)
  tnai.checks("exponent",exponent)
  tnai.checks("verbose",verbose)
  ##-----filter 'listOfRegulons' by 'minRegulonSize'
  gs.size <- unlist(
    lapply(lapply(listOfRegulons, intersect, y = names(phenotype)), length)
  )
  max.size <- max(gs.size)
  ##-----stop when no gene set passes the size requirement
  if(all(unlist(lapply(listOfRegulons,length))==0)){
    tp<-" overlapped genes with the universe!\n The largest number of overlapped genes is: "
    stop(paste("NOTE: no regulon has >= ", minRegulonSize, tp, max.size, sep=""))
  }
  ##-----get filtered list
  gs.id <- which(gs.size >= minRegulonSize)
  gs.size <- gs.size[gs.id]
  listOfRegulons <- listOfRegulons[gs.id]
  ##-----calculate enrichment scores for all regulons
  test.collection<-list()
  if(length(listOfRegulons) > 0){
    test.collection <- gsea1tna(listOfRegulons, 
                                phenotype=phenotype,exponent=exponent,
                                nPermutations=nPermutations, 
                                minRegulonSize=minRegulonSize,verbose=verbose)
  } else {
    test.collection<-list(Observed.scores=NULL, Permutation.scores=NULL)
  }
  ##-----compute pvals
  if(length(test.collection$Observed.scores) > 0) {
    test.pvalues.collection <- tna.permutation.pvalues(
      permScores = test.collection$Permutation.scores,
      dataScores = test.collection$Observed.scores)  	
    gsea.adjust.pval <- p.adjust(test.pvalues.collection, 
                                 method = pAdjustMethod)
    GSEA1.results <- cbind(gs.size, test.collection$Observed.scores,
                          test.pvalues.collection, gsea.adjust.pval)			
    colnames(GSEA1.results) <- c("Regulon.Size", "Observed.Score", "Pvalue", "Adjusted.Pvalue")
  } else {
    GSEA1.results <- matrix(, nrow=0, ncol=4)
    colnames(GSEA1.results) <- c("Regulon.Size", "Observed.Score", "Pvalue", "Adjusted.Pvalue")
  }
  #-----format results 
  GSEA1.results[,"Observed.Score"]<-round(GSEA1.results[,"Observed.Score"],2)
  GSEA1.results[,"Pvalue"]<-signif(GSEA1.results[,"Pvalue"], digits=5)
  GSEA1.results[,"Adjusted.Pvalue"]<-signif(GSEA1.results[,"Adjusted.Pvalue"], digits=5)
  if(verbose)cat("-Gene set enrichment analysis complete \n\n")
  return( GSEA1.results ) 
}

##------------------------------------------------------------------------------
##This function takes a list of gene sets (regulons), a named phenotype 
##vector, and returns the results of gene set enrichment analysis for all 
##regulons (with multiple hypothesis testing correction).
run.gsea2 <- function(listOfRegulonsAndMode, phenotype, pAdjustMethod="BH", 
                      pValueCutoff=0.05, nPermutations=1000, 
                      exponent=1, verbose=TRUE) {
  ##-----calculate enrichment scores for all regulons
  test.collection.up<-list()
  test.collection.down<-list()
  if(length(listOfRegulonsAndMode) > 0){
    listOfRegulonsUp <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp>0])
    })
    names(listOfRegulonsUp)<-names(listOfRegulonsAndMode)
    listOfRegulonsDown <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp<0])
    })
    names(listOfRegulonsDown)<-names(listOfRegulonsAndMode)
    gs.size.up <- unlist(lapply(lapply(listOfRegulonsUp, intersect, y = names(phenotype)), length))
    gs.size.down <- unlist(lapply(lapply(listOfRegulonsDown, intersect, y = names(phenotype)), length))
    test.collection.up <- gsea2tna(listOfRegulonsUp, phenotype=phenotype,exponent=exponent,
                                   nPermutations=nPermutations, verbose1=verbose, verbose2=TRUE)
    test.collection.down <- gsea2tna(listOfRegulonsDown, phenotype=phenotype,exponent=exponent,
                                     nPermutations=nPermutations, verbose1=verbose, verbose2=FALSE)
  } else {
    test.collection.up<-list(Observed.scores=NULL, Permutation.scores=NULL)
    test.collection.down<-list(Observed.scores=NULL, Permutation.scores=NULL)
  }
  ##-----compute pvals
  b1<-length(test.collection.up$Observed.scores)>0 && length(test.collection.down$Observed.scores)>0
  b2<-length(test.collection.up$Observed.scores)==length(test.collection.down$Observed.scores)
  if(b1 && b2) {
    pvalues.up <- tna.permutation.pvalues(permScores = test.collection.up$Permutation.scores,
                                          dataScores = test.collection.up$Observed.scores)    
    adjust.pval.up <- p.adjust(pvalues.up, method = pAdjustMethod, n=length(pvalues.up)*2)
    GSEA2.results.up <- cbind(gs.size.up, test.collection.up$Observed.scores, pvalues.up, adjust.pval.up)
    
    pvalues.down <- tna.permutation.pvalues(permScores = test.collection.down$Permutation.scores,
                                            dataScores = test.collection.down$Observed.scores)    
    adjust.pval.down <- p.adjust(pvalues.down, method = pAdjustMethod, n=length(pvalues.down)*2)
    GSEA2.results.down <- cbind(gs.size.down,test.collection.down$Observed.scores,pvalues.down, adjust.pval.down)
    
    test.collection.both<-list(Permutation.scores=test.collection.up$Permutation.scores-test.collection.down$Permutation.scores,
                               Observed.scores=test.collection.up$Observed.scores-test.collection.down$Observed.scores)
    pvalues.both <- tna.permutation.pvalues(permScores = test.collection.both$Permutation.scores,
                                            dataScores = test.collection.both$Observed.scores) 
    adjust.pval.both <- p.adjust(pvalues.both, method = pAdjustMethod)
    GSEA2.results.both <- cbind(gs.size.up+gs.size.down,test.collection.both$Observed.scores,pvalues.both, adjust.pval.both)
  } else {
    GSEA2.results.up <- matrix(, nrow=0, ncol=4)
    GSEA2.results.down <- matrix(, nrow=0, ncol=4)
    GSEA2.results.both <- matrix(, nrow=0, ncol=4)
  }
  colnames(GSEA2.results.up) <- c("Regulon.Size", "Observed.Score", "Pvalue", "Adjusted.Pvalue")
  colnames(GSEA2.results.down) <- c("Regulon.Size", "Observed.Score", "Pvalue", "Adjusted.Pvalue")
  colnames(GSEA2.results.both) <- c("Regulon.Size", "Observed.Score", "Pvalue", "Adjusted.Pvalue")
  #-----format results
  GSEA2.results.up[,"Observed.Score"]<-round(GSEA2.results.up[,"Observed.Score"],2)
  GSEA2.results.up[,"Pvalue"]<-signif(GSEA2.results.up[,"Pvalue"], digits=5)
  GSEA2.results.up[,"Adjusted.Pvalue"]<-signif(GSEA2.results.up[,"Adjusted.Pvalue"], digits=5)
  GSEA2.results.down[,"Observed.Score"]<-round(GSEA2.results.down[,"Observed.Score"],2)
  GSEA2.results.down[,"Pvalue"]<-signif(GSEA2.results.down[,"Pvalue"], digits=5)
  GSEA2.results.down[,"Adjusted.Pvalue"]<-signif(GSEA2.results.down[,"Adjusted.Pvalue"], digits=5)
  GSEA2.results.both[,"Observed.Score"]<-round(GSEA2.results.both[,"Observed.Score"],2)
  GSEA2.results.both[,"Pvalue"]<-signif(GSEA2.results.both[,"Pvalue"], digits=5)
  GSEA2.results.both[,"Adjusted.Pvalue"]<-signif(GSEA2.results.both[,"Adjusted.Pvalue"], digits=5)
  #-----pack and return
  GSEA2.results.up<-data.frame(Regulon=rownames(GSEA2.results.up),GSEA2.results.up,stringsAsFactors=FALSE)
  GSEA2.results.down<-data.frame(Regulon=rownames(GSEA2.results.down),GSEA2.results.down,stringsAsFactors=FALSE)
  GSEA2.results.both<-data.frame(Regulon=rownames(GSEA2.results.both),GSEA2.results.both,stringsAsFactors=FALSE)
  GSEA2.results<-list(positive=GSEA2.results.up,negative=GSEA2.results.down,differential=GSEA2.results.both)
  return( GSEA2.results ) 
}

##------------------------------------------------------------------------------
##This function takes a list of regulons, a named phenotype vector,
##and returns results from synergy analysis
run.synergy <- function(collectionsOfPairsR1R2, labpair, phenotype, pAdjustMethod="BH", 
                       pValueCutoff=0.05, minIntersectSize=1, nPermutations=1000, exponent=1, 
                       verbose=TRUE) {
  ##check arguments
  tnai.checks("rgcs",collectionsOfPairsR1R2)
  tnai.checks("labpair",labpair)
  tnai.checks("phenotype",phenotype)
  tnai.checks("pAdjustMethod",pAdjustMethod)
  tnai.checks("pValueCutoff",pValueCutoff)
  tnai.checks(name="minIntersectSize",para=minIntersectSize)
  tnai.checks("nPermutations",nPermutations)
  tnai.checks("exponent",exponent)
  tnai.checks("verbose",verbose)
  ##calculate enrichment scores
  regize<-sapply(1:length(collectionsOfPairsR1R2),function(i){
    unlist(lapply(collectionsOfPairsR1R2[[i]],length))
  })
  test.pairs <- pairwiseSynergy(
    collectionsOfPairs=collectionsOfPairsR1R2,
    phenotype=phenotype,exponent=exponent,
    nPermutations=nPermutations, 
    minIntersectSize=minIntersectSize, 
    verbose=verbose)
  ES.Intersect <- test.pairs$Subset.scores
  ES.Union <- apply(test.pairs$Allset.scores,1,median)
  EffectSize <- log2(abs(ES.Intersect)/abs(ES.Union))
  Pvalue <- tna.permutation.pvalues(test.pairs$Allset.scores,test.pairs$Subset.scores)
  synergy.test<-cbind(t(regize),ES.Union,ES.Intersect,EffectSize,Pvalue)
  #add names
  synergy.test<-data.frame(labpair,synergy.test,stringsAsFactors=FALSE)
  colnames(synergy.test)<-c("Regulon1","Regulon2","Union","Intersect",
                            "ES.Union","ES.Intersect", "EffectSize", "Pvalue")
  #adjust pvalue
  padj <- p.adjust(synergy.test[,"Pvalue"], method=pAdjustMethod)
  synergy.results<-cbind(synergy.test, Adjusted.Pvalue=padj)
  #-----format and return results
  synergy.results[,"ES.Union"]<-round(synergy.results[,"ES.Union"],2)
  synergy.results[,"ES.Intersect"]<-round(synergy.results[,"ES.Intersect"],2)
  synergy.results[,"EffectSize"]<-round(synergy.results[,"EffectSize"],4)
  synergy.results[,"Pvalue"]<-signif(synergy.results[,"Pvalue"], digits=4)
  synergy.results[,"Adjusted.Pvalue"]<-signif(synergy.results[,"Adjusted.Pvalue"], digits=4)
  if(verbose)cat("-Synergy analysis complete \n\n")
  return(synergy.results)
}

##------------------------------------------------------------------------------
##This function takes a list of regulons, a named phenotype vector,
##and returns results from shadow analysis
run.shadow <- function(collectionsOfPairsR1, collectionsOfPairsR2, labpair, phenotype, 
                       pAdjustMethod="BH", pValueCutoff=0.05, nPermutations=1000, 
                       minIntersectSize=1, exponent=1, verbose=TRUE) {
  ##check arguments
  tnai.checks("rgcs",collectionsOfPairsR1)
  tnai.checks("rgcs",collectionsOfPairsR2)
  tnai.checks("phenotype",phenotype)
  tnai.checks("pAdjustMethod",pAdjustMethod)
  tnai.checks("pValueCutoff",pValueCutoff)
  tnai.checks("nPermutations",nPermutations)
  tnai.checks(name="minIntersectSize",para=minIntersectSize)
  tnai.checks("exponent",exponent)
  tnai.checks("verbose",verbose)
  tnai.checks("labpair",labpair)
  if(length(collectionsOfPairsR1)!=length(collectionsOfPairsR2)){
    stop("NOTE: 'collectionsOfPairsR1' and 'collectionsOfPairsR2' should have the same length!")
  }
  
  ##check if package snow has been loaded and 
  ##a cluster object has been created
  b1<-"package:snow" %in% search()
  b2<-tryCatch({
    cl<-getOption("cluster")
    cl.check<-FALSE
    if(is(cl, "cluster")){
      cl.check <- all( sapply(1:length(cl),function(i)isOpen(cl[[i]]$con) ) == TRUE )
    }
    cl.check
  }, error=function(e){ FALSE 
  })
  if( b1 && b2) {
    if(verbose)cat("-Performing shadow analysis (parallel version - ProgressBar not available) ...\n")
  } else {
    if(verbose)cat("-Performing shadow analysis ...\n")
  }
  if(verbose)cat("--For", length(collectionsOfPairsR1), "regulon pairs ...\n")
  
  ##compute enrichment scores for 1st shadow candidates
  regize<-sapply(1:length(collectionsOfPairsR1),function(i){
    unlist(lapply(collectionsOfPairsR1[[i]],length))
  })
  test.pairs <- pairwiseShadow(
    collectionsOfPairs=collectionsOfPairsR1,
    phenotype=phenotype,exponent=exponent,
    nPermutations=nPermutations, 
    minIntersectSize=minIntersectSize,verbose=verbose)
  ES.Unique <- test.pairs$Subset.scores
  ES.All <- apply(test.pairs$Allset.scores,1,median)
  EffectSize <- log2(abs(ES.Unique)/abs(ES.All))
  Pvalue <- tna.permutation.pvalues.shadow(test.pairs$Allset.scores,test.pairs$Subset.scores)
  shadow.test<-cbind(t(regize),ES.All,ES.Unique,EffectSize,Pvalue)
  #add names
  shadowR1.test<-data.frame(labpair,shadow.test,stringsAsFactors=FALSE)[c(2,1,3:8)]
  colnames(shadowR1.test)<-c("Regulon","Shadow","All","Unique",
                             "ES.All","ES.Unique", "EffectSize", "Pvalue")  
  
  ##compute enrichment scores for 2nd shadow candidates
  regize<-sapply(1:length(collectionsOfPairsR2),function(i){
    unlist(lapply(collectionsOfPairsR2[[i]],length))
  })  
  test.pairs <- pairwiseShadow(
    collectionsOfPairs=collectionsOfPairsR2,
    phenotype=phenotype,exponent=exponent,
    nPermutations=nPermutations, 
    minIntersectSize=minIntersectSize,verbose=verbose)
  ES.Unique <- test.pairs$Subset.scores
  ES.All <- apply(test.pairs$Allset.scores,1,median)
  EffectSize <- log2(abs(ES.Unique)/abs(ES.All))
  Pvalue <- tna.permutation.pvalues.shadow(test.pairs$Allset.scores,test.pairs$Subset.scores)
  shadow.test<-cbind(t(regize),ES.All,ES.Unique,EffectSize,Pvalue)
  #add names
  shadowR2.test<-data.frame(labpair,shadow.test,stringsAsFactors=FALSE)
  colnames(shadowR2.test)<-c("Regulon","Shadow","All","Unique",
                             "ES.All","ES.Unique","EffectSize", "Pvalue")
  #adjust pvalue
  p.adj <- p.adjust(shadowR1.test[,"Pvalue"], method = pAdjustMethod)
  shadowR1.test<-cbind(shadowR1.test,Adjusted.Pvalue=p.adj)  
  p.adj <- p.adjust(shadowR2.test[,"Pvalue"], method = pAdjustMethod)
  shadowR2.test<-cbind(shadowR2.test,Adjusted.Pvalue=p.adj)
  
  #identify shadowing with sig. pvalues < pValueCutoff
  shadow.results<-matrix(NA,nrow(shadowR1.test),ncol(shadowR1.test))
  colnames(shadow.results)<-colnames(shadowR1.test)
  shadow.results<-data.frame(shadow.results,check.names=FALSE, stringsAsFactors=FALSE)
  rownames(shadow.results)<-paste(shadowR1.test[[1]],shadowR1.test[[2]],sep=".")
  if(nrow(shadowR1.test)>0) {
    a1<-which(shadowR1.test[,"Adjusted.Pvalue"]<=pValueCutoff & shadowR1.test[,"EffectSize"]<0)
    a2<-which(shadowR2.test[,"Adjusted.Pvalue"]<=pValueCutoff & shadowR2.test[,"EffectSize"]<0)
    #setdiff excludes mutual shadow!
    shadow.a1<-setdiff(a1,a2)
    shadow.a2<-setdiff(a2,a1)  
    if(length(shadow.a2)>0){
      shadow.results[shadow.a2,"Regulon"]<-shadowR2.test[shadow.a2,"Regulon"]
      shadow.results[shadow.a2,"Shadow"]<-shadowR2.test[shadow.a2,"Shadow"]
      shadow.results[shadow.a2,"All"]<-shadowR2.test[shadow.a2,"All"]
      shadow.results[shadow.a2,"Unique"]<-shadowR2.test[shadow.a2,"Unique"]
      shadow.results[shadow.a2,"ES.All"]<-shadowR2.test[shadow.a2,"ES.All"]
      shadow.results[shadow.a2,"ES.Unique"]<-shadowR2.test[shadow.a2,"ES.Unique"]
      shadow.results[shadow.a2,"EffectSize"]<-shadowR2.test[shadow.a2,"EffectSize"]
      shadow.results[shadow.a2,"Pvalue"]<-shadowR2.test[shadow.a2,"Pvalue"]
      shadow.results[shadow.a2,"Adjusted.Pvalue"]<-shadowR2.test[shadow.a2,"Adjusted.Pvalue"]
    }
    if(length(shadow.a1)>0){
      shadow.results[shadow.a1,"Regulon"]<-shadowR1.test[shadow.a1,"Regulon"]
      shadow.results[shadow.a1,"Shadow"]<-shadowR1.test[shadow.a1,"Shadow"]
      shadow.results[shadow.a1,"All"]<-shadowR1.test[shadow.a1,"All"]
      shadow.results[shadow.a1,"Unique"]<-shadowR1.test[shadow.a1,"Unique"]
      shadow.results[shadow.a1,"ES.All"]<-shadowR1.test[shadow.a1,"ES.All"]
      shadow.results[shadow.a1,"ES.Unique"]<-shadowR1.test[shadow.a1,"ES.Unique"]
      shadow.results[shadow.a1,"EffectSize"]<-shadowR1.test[shadow.a1,"EffectSize"]
      shadow.results[shadow.a1,"Pvalue"]<-shadowR1.test[shadow.a1,"Pvalue"]
      shadow.results[shadow.a1,"Adjusted.Pvalue"]<-shadowR1.test[shadow.a1,"Adjusted.Pvalue"]     
    }
    shadow.results<-shadow.results[!is.na(shadow.results$Shadow),,drop=FALSE]
  }
  #-----format and return results
  shadow.results[,"ES.All"]<-round(shadow.results[,"ES.All"],2)
  shadow.results[,"ES.Unique"]<-round(shadow.results[,"ES.Unique"],2)
  shadow.results[,"EffectSize"]<-round(shadow.results[,"EffectSize"],4)
  shadow.results[,"Pvalue"]<-signif(shadow.results[,"Pvalue"], digits=4)
  shadow.results[,"Adjusted.Pvalue"]<-signif(shadow.results[,"Adjusted.Pvalue"], digits=4)
  if(verbose)cat("-Shadow analysis complete \n\n")
  return(list(shadowR1=shadowR1.test, shadowR2=shadowR2.test, results=shadow.results))
}

##------------------------------------------------------------------------------
##CMAP for GSEA2
run.tna.cmap <- function(listOfRegulonsAndMode, phenotype, exponent) {
  if(length(listOfRegulonsAndMode) > 0){
    phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
    listOfRegulonsUp <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp>0])
    })
    names(listOfRegulonsUp)<-names(listOfRegulonsAndMode)
    listOfRegulonsDown <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp<0])
    })
    names(listOfRegulonsDown)<-names(listOfRegulonsAndMode)
    scoreMatUp<-sapply(names(listOfRegulonsUp),function(reg){
      gseaScores4CMAP(phenotype,listOfRegulonsUp[[reg]],exponent)
    })
    scoreMatDown<-sapply(names(listOfRegulonsDown),function(reg){
      gseaScores4CMAP(phenotype,listOfRegulonsDown[[reg]],exponent)
    })
    diffScoreMat<-scoreMatUp-scoreMatDown
  } else {
    diffScoreMat<-numeric()
  }
  return(diffScoreMat)
}

