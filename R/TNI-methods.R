################################################################################
##########################         TNI Class        ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "TNI",
          function(.Object, gexp, transcriptionFactors) {
            ##-----check arguments
            if(missing(gexp))stop("NOTE: 'gexp' is missing!")    
            if(missing(transcriptionFactors))stop("NOTE: 'transcriptionFactors' is missing!")            
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
            colnames(sum.info.para$cdt)<-c("sampling", "pValueCutoff","pAdjustMethod","statFilter","statUniverse",
                                           "minRegulonSize","minIntersectSize","miThreshold")
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
    if(verbose)cat("-Preprocessing for input data ...\n")
    ##-----check gexpIDs if available
    if(!is.null(gexpIDs)){
      if(verbose)cat("--Mapping 'gexp' to 'gexpIDs' ...\n")
      if( any(!rownames(object@gexp)%in%rownames(gexpIDs)) ){
        stop("NOTE: all rownames in 'gexp' should be available in 'gexpIDs'!")
      }
      if(cvfilter){
        if(verbose)cat("--Removing duplicated genes (keep max coefficient of variation!) ...\n")
        #i.e. col1=probe, col2=gene (collapse cv by col2)
        cvres<-cv.filter(object@gexp, gexpIDs)
        object@gexp<-cvres$gexp
        object@annotation<-cvres$ids
      } else {
        #or leave by the user!!!
        object@annotation<-gexpIDs[rownames(object@gexp),]
      }
      #correct 'symbol' name if a valid label
      labs<-c("SYMBOL","symbol","Symbol")
      idx<-which(labs%in%colnames(object@annotation))
      if(length(idx)>0){
        if(all(idx!=1)){
          labs<-labs[idx]
          idx<-which(colnames(object@annotation)==labs)
          colnames(object@annotation)[idx]<-"SYMBOL"
        }
        #..remove any empty space or NA from SYMBOL!!!
        idx<-is.na(object@annotation$SYMBOL)
        object@annotation$SYMBOL[idx]<-rownames(object@annotation)[idx]
        idx<-object@annotation$SYMBOL==""|object@annotation$SYMBOL=="NA"
        object@annotation$SYMBOL[idx]<-rownames(object@annotation)[idx]
      }
    }
    #-----check TFs in gexp
    object@summary$tfs[,"input"]<-length(object@transcriptionFactors)
    if(verbose) cat("--Checking TFs in the dataset ...\n")
    idxtfs<-object@transcriptionFactors%in%rownames(object@gexp)
    object@transcriptionFactors<-object@transcriptionFactors[idxtfs]
    object@summary$tfs[,"valid"]<-length(object@transcriptionFactors)
    if(length(object@transcriptionFactors)==0)stop("NOTE: input 'transcriptionFactors' contains no useful data!\n")
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
          warning("one or more symbols in the named vector 'transcriptionFactors' seem to differ from 'gexpIDs'!")
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
    if(!is.null(object@results$conditional)>0){
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
  function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE) {
    ##-----reset compatibility with old args
    if(what=="tn.dpi")what="tnet"
    if(what=="tn.ref")what="refnet"
    ##-----check input arguments
    tnai.checks(name="tni.what",para=what)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="report",para=report)
    ##-----get query
    query<-NULL  
    if(what=="gexp"){
      query<-object@gexp
    } else if(what=="tfs"){
      query<-object@transcriptionFactors
    } else if(what=="para"){
      query<-object@para
    } else if(what=="refnet"){
      query<-object@results$tn.ref
    } else if(what=="tnet"){
      query<-object@results$tn.dpi
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
      }
    } else if(what=="cdt"){
      query<-object@results$conditional$count
      if(is.null(ntop)){
        for(md in names(query)){
          qry<-query[[md]]
          qry<-qry[qry[,"IConstraint"]==0,-2,drop=FALSE]
          qry<-qry[qry[,"RConstraint"]==0,-2,drop=FALSE]
          qry<-qry[qry[,"AdjustedPvalue"] <= object@para$cdt$pValueCutoff,,drop=FALSE]
          if(order && nrow(qry)>1)qry<-qry[sort.list(qry[,"Pvalue"]),,drop=FALSE]
          query[[md]]<-qry
        }
      } else {
        for(md in names(query)){
          qry<-query[[md]]
          qry<-qry[qry[,"IConstraint"]==0,-2,drop=FALSE]
          qry<-qry[qry[,"RConstraint"]==0,-2,drop=FALSE]
          qryntop<-ntop
          if(nrow(qry)>1){
            if(qryntop>nrow(qry) || qryntop<0)qryntop=nrow(qry)
            idx<-sort.list(qry[,"Pvalue"]) 
            qry<-qry[idx[1:qryntop],,drop=FALSE]
            query[[md]]<-qry
          }
        }
      }
      if(reportNames){
        for(md in names(query)){
          if(nrow(query[[md]])>0){
            idx<-match(query[[md]][,"Modulator"],object@modulators)
            query[[md]][,"Modulator"]<-names(object@modulators)[idx]
            idx<-match(query[[md]][,"TF"],object@transcriptionFactors)
            query[[md]][,"TF"]<-names(object@transcriptionFactors)[idx]
          }
        }
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
##get graph from TNI 
setMethod(
  "tni.graph",
  "TNI",
  function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, amapFilter="quantile", amapCutoff=NULL){
    # chech igraph compatibility
    b1<-"package:igraph0" %in% search()
    b2<- "igraph0" %in%  loadedNamespaces()
    if( b1 || b2) {
      stop("\n\n ...conflict with 'igraph0': please use the new 'igraph' package!")
    }
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!")
    tnai.checks(name="tni.gtype",para=gtype)
    if(gtype=="mmap"){ #get modulatory maps
      ##-----check input arguments
      if(object@status["Conditional"]!="[x]")stop("NOTE: input need conditional analysis!")
      #get tfs and modulators
      cdt<-tni.get(object,what="cdt")
      testedtfs<-names(cdt)
      testedtfs<-object@transcriptionFactors[object@transcriptionFactors%in%testedtfs]
      modulators<-sapply(testedtfs,function(tf){
        rownames(cdt[[tf]])
      })
      modulators<-unlist(modulators)
      modulators<-object@modulators[object@modulators%in%modulators]
      othertfs<-object@transcriptionFactors
      othertfs<-othertfs[!othertfs%in%testedtfs]
      othertfs<-othertfs[othertfs%in%modulators]
      ##-----get tnet
      if(tnet=="ref"){
        tnet<-object@results$tn.ref
      } else if(tnet=="dpi"){
        tnet<-object@results$tn.dpi
      }
      #get adjmt
      tnet<-tnet[unique(c(testedtfs,setdiff(modulators,testedtfs))),testedtfs,drop=FALSE]
      mnet<-tnet;mnet[,]=0
      sapply(colnames(mnet),function(i){
        tp<-cdt[[i]]
        mnet[rownames(tp),i]<<-tp$Mode
        NULL
      })
      #get graph
      g<-tni.mmap(mnet,tnet,othertfs)
      #add annotation
      if(.hasSlot(object, "annotation")){
        if(nrow(object@annotation)>0)g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      }
      #set names if available
      if(!is.null(V(g)$SYMBOL)){
        g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
      } else {
        V(g)$nodeAlias<-V(g)$name
      }
      #map tested tfs and modulators
      V(g)$tfs<-as.numeric(V(g)$name%in%testedtfs)
      V(g)$modulators<-as.numeric(V(g)$name%in%modulators)
      #get regulon size (downstream targets)
      rsize<-tni.get(object,what="regulons")
      rsize<-sapply(rsize,length)
      rsize<-rsize[names(rsize)%in%V(g)$name]
      rsize<-rsize[names(rsize)%in%testedtfs]
      idx<-match(names(rsize),V(g)$name)
      V(g)$regulonSize<-NA
      V(g)$regulonSize[idx]<-rsize
      #test if node size can be represented by quantiles
      nquant=5
      nq<-length(unique(quantile(V(g)$regulonSize,na.rm=TRUE)))
      if(nq<3)nquant=NULL
      g<-att.setv(g=g, from="regulonSize", to='nodeSize',xlim = c(35, 70, 15),nquant=nquant,title="Downstream targets (n)")
      #set node attr
      V(g)$nodeFontSize=15
      V(g)$nodeLineWidth<-1.5
      V(g)$nodeColor<-"grey"
      V(g)$nodeLineColor<-"grey"
      V(g)$nodeFontSize[V(g)$tfs==1]=25
      #set node shape and legend
      g<-att.setv(g=g, from="tfs", to='nodeShape',shapes=c("DIAMOND","ELLIPSE"),title="")
      g$legNodeShape$legend<-c("Modulator","TF")
      if(ecount(g)>0){
        #set edge attr
        E(g)$edgeWidth<-1.5
        #c("#7DA7FF","#63DB93","#FE7276")
        g<-att.sete(g=g, from="modeOfAction", to='edgeColor',cols=c("#96D1FF","grey80","#FF8E91"), 
                    title="ModeOfAction",categvec=-1:1)
        g$legEdgeColor$legend<-c("Down","NA","Up")
        #map upstring edges
        dg<-degree(g)
        el<-get.edgelist(g,names=FALSE)
        ew<-sapply(1:nrow(el),function(i){
          min(dg[el[i,1]],dg[el[i,2]])
        })
        E(g)$terminalEdge<-0
        E(g)$terminalEdge[ew==1]<-1
        #map modeOfAction to node attribute (compute the median mode of action)
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
    } else {
      ##-----check input arguments
      if(object@status["DPI.filter"]!="[x]")stop("NOTE: input data need dpi analysis!")
      tnai.checks(name="tnet",para=tnet)
      tnai.checks(name="tfs",para=tfs)
      tnai.checks(name="minRegulonSize",para=minRegulonSize)
      tnai.checks(name="amapFilter",para=amapFilter)
      tnai.checks(name="amapCutoff",para=amapCutoff)
      ##-----get tnet
      if(tnet=="ref"){
        tnet<-object@results$tn.ref
      } else if(tnet=="dpi"){
        tnet<-object@results$tn.dpi
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
    }
    return(g)
  }
)
##------------------------------------------------------------------------------
##run conditional mutual information analysis
setMethod(
  "tni.conditional",
  "TNI",
  function(object, modulators=NULL, tfs=NULL, sampling=35, pValueCutoff=0.01, 
           pAdjustMethod="bonferroni", statFilter="phyper", statUniverse="all", 
           minRegulonSize=15, minIntersectSize=5, miThreshold=NULL,
           parChunks=10, verbose=TRUE){
    ##-----check input arguments
    if(object@status["DPI.filter"]!="[x]")stop("NOTE: input data need dpi analysis!")
    tnai.checks(name="modulators",para=modulators)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="sampling",para=sampling)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="statFilter",para=statFilter)
    tnai.checks(name="statUniverse",para=statUniverse)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="minIntersectSize",para=minIntersectSize)
    tnai.checks(name="miThreshold",para=miThreshold)
    tnai.checks(name="parChunks",para=parChunks)
    tnai.checks(name="verbose",para=verbose)
    ##-----par info
    object@para$cdt<-list(sampling=sampling, pValueCutoff=pValueCutoff,
                          pAdjustMethod=pAdjustMethod, statFilter=statFilter,
                          statUniverse=statUniverse,
                          minRegulonSize=minRegulonSize,
                          minIntersectSize=minIntersectSize,
                          miThreshold=NA)
    ##-----summary info
    cdt<-unlist(object@para$cdt)
    object@summary$para$cdt<-matrix(cdt,nrow=1,ncol=8)
    rownames(object@summary$para$cdt)<-"Parameter"
    colnames(object@summary$para$cdt)<-names(cdt)

    if(verbose)cat("-Preprocessing for input data ...\n")
    ##-----make sure all 'tfs' are valid
    if(is.null(tfs)){
      tfs<-object@transcriptionFactors
    } else {
      if(verbose)cat("--Checking TFs in the dataset ...\n")
      tfs<-as.character(tfs)
      idx<-which(!tfs%in%object@transcriptionFactors)
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
      if(verbose)cat("--Checking modulators in the dataset ...\n")
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
      if(length(modulators)==1)stop("NOTE: incorrect number of dimensions: range constraint step requires >2 valid modulators!")
      #final check: remove any empty space or NA from md names!!!
      mdnames<-names(modulators)
      idx<-mdnames==""|mdnames=="NA"
      names(modulators)[idx]<-modulators[idx] 
    }
    object@modulators<-modulators
    
    #     ##---
    #     if(verbose && dpiFilter && length(modulators)>30){
    #       cat("\nNote: 'dpiFilter=TRUE' can be computationally costly!\n")
    #       cat("Continue (y/n)? \n")
    #       b1<-readLines(n=1)
    #       if(b1!="y"){
    #         cat("Analysis canceled! \n")
    #         return(object)
    #       }
    #     }
    
    ##-----get TF-targets from tnet
    if(verbose)cat("--Extracting TF-targets ...\n")
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
    ##-----Checking independence of modulators and TFs
    if(verbose)cat("--Applying modulator independence constraint ...\n")
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
    if(verbose)cat("--Applying modulator range constraint ...\n")
    gxtemp<-object@gexp
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)
    gxtemp<-t(apply(gxtemp[modulators,,drop=FALSE],1,sort))[,c(idxLow,idxHigh),drop=FALSE]
    ##--run limma on both tails to filter-out modulators
    t <- factor(c(rep("low",spsz),rep("high",spsz)))
    design <- model.matrix(~0+t)
    fit <- lmFit(gxtemp,design)
    thigh=tlow=NULL
    contrasts <- makeContrasts(thigh-tlow, levels=design)
    ct.fit <- eBayes(contrasts.fit(fit, contrasts))
    res.fit<-unclass(decideTests(ct.fit, adjust.method=object@para$cdt$pAdjustMethod, p.value=object@para$cdt$pValueCutoff))
    RConstraintList<-rownames(res.fit)[res.fit<=0]
    ##--get samples (sorted index) for each pre-selected modulator
    if(verbose)cat("--Selecting subsamples ...\n")
    gxtemp<-object@gexp
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)    
    idx<-t(apply(gxtemp[modulators,],1,sort.list))
    idxLow<-idx[,idxLow]
    idxHigh<-idx[,idxHigh]
    ##-----estimate mutual information threshold
    if(is.null(miThreshold)){
      if(verbose)cat("\n")
      mmark<-tni.conditional.threshold(object, ntfs=length(tfs), nsamples=spsz, parChunks=parChunks, verbose=verbose)
      if(verbose)cat("\n")
    } else {
      #0.2765982/0.03652542/0.04254874
      mmark<-miThreshold
    }
    ##-----update para/summary
    object@para$cdt$miThreshold<-mmark
    object@para$cdt$miThreshold<-mmark
    object@summary$para$cdt[,"miThreshold"]<-mmark
    ##-----start conditional mutual information analysis
    if(verbose)cat("\n")
    if(verbose)cat("-Performing conditional mutual information analysis ...\n")
    if(verbose)cat("--For", length(tfs), "tfs and" , length(modulators), "candidate modulators \n")
    #set data object to save results
    reseffect<-lapply(tfTargets,function(tar){
      data.frame(tagets=tar,stringsAsFactors=FALSE)
    })
    rescount<-lapply(tfTargets,function(tar){
      res<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,stringsAsFactors=FALSE)
      colnames(res)<-c("Modulator","IConstraint","RConstraint","TF","UniverseSize","EffectSize","RegulonSize","Expected",
                       "Observed","Negative","Positive","Mode","Pvalue","AdjustedPvalue")
      res
    })
    #run analysis
    minExpectedSize<-ifelse(statFilter=="phyper",1,5)
    if(verbose) pb <- txtProgressBar(style=3)
    for(i in 1:length(modulators)){
      md<-modulators[i]
      #get sample ordering
      lw<-idxLow[md,]
      hg<-idxHigh[md,]
      #compute conditional on both tails
      milow<-tni.pmim(gxtemp[,lw],tfs,estimator=object@para$perm$estimator)
      mihigh<-tni.pmim(gxtemp[,hg],tfs,estimator=object@para$perm$estimator)
      #get delta mi
      miDelta<-mihigh-milow
      #remove associations below mi threshold
      milow[milow<=mmark]=0
      mihigh[mihigh<=mmark]=0
      #identify significant negative and positive modulations
      neg <- milow > mmark & miDelta < -mmark & mihigh <= mmark
      pos <- mihigh > mmark & miDelta > mmark & milow <= mmark
      #remove modulation below mi threshold
      miDelta[!neg & !pos]=0 
      #...used to compute the sample space and expected values!
      mTargets<-milow+mihigh
      mTargets[mTargets>0]=1
      #get final results
      sapply(tfs,function(tf){
        dtmtx<-miDelta[,tf]
        mtmtx<-mTargets[,tf]
        tftar<-tfTargets[[tf]]
        tfalltar<-tfAllTargets[[tf]]
        #---get SZs
        if(statUniverse=="all"){
          #todos modulaveis
          EffectSZ<-sum(mtmtx!=0)
          #todos testados
          UniSZ<-length(mtmtx)
        } else {
          #todos modulaveis na tnet
          EffectSZ<-sum(mtmtx[tfalltar]!=0)
          #todos testados na tnet
          UniSZ<-sum(rowSums(object@results$tn.dpi!=0)>0)
        }
        RegSZ<-length(tftar)
        ExpOV<-(EffectSZ*RegSZ)/UniSZ
        ObsOV<-sum(dtmtx[tftar]!=0)
        if((ObsOV/RegSZ*100)>=minIntersectSize && ExpOV>=minExpectedSize){
          #apply exclusion list (independence/range constraint)
          iconst<-ifelse(md%in%IConstraintList[[tf]],1,0)
          rconst<-ifelse(md%in%RConstraintList,1,0)
          if(iconst || rconst){
            EffectSZ<-NA;ExpOV<-NA;ObsOV<-NA;ObsPos<-NA
            ObsNeg<-NA;Mode<-NA;pvalue<-NA
          } else {
            #---get mode of actions
            ObsPos<-sum(dtmtx[tftar]>0)
            ObsNeg<-sum(dtmtx[tftar]<0)
            Mode<-if(ObsNeg>ObsPos) -1 else if(ObsNeg<ObsPos) 1 else 0
            #---resolve obs
            Obs<-if(Mode==-1) abs(ObsNeg) else if(Mode==1) ObsPos else ObsOV
            if(statFilter=="phyper"){
              #---run FET with phyper (obs-1)
              pvalue <- phyper(Obs-1, RegSZ, UniSZ-RegSZ, EffectSZ, lower.tail=FALSE)
            } else {
              #---chisq test
              if(Obs>ExpOV){
                yt <- min(0.5,Obs-ExpOV)
                chist<-((Obs-ExpOV)-yt)^2/ExpOV
              } else {
                chist<-0
              }
              pvalue<-pchisq(chist, df=1, lower.tail=FALSE)
            }
          }
          #---add results to a list
          reseffect[[tf]][[md]]<<-dtmtx[tftar]
          rescount[[tf]][md,]<<-c(NA,NA,NA,NA,UniSZ,EffectSZ,RegSZ,ExpOV,ObsOV,ObsNeg,ObsPos,Mode,pvalue,NA)
          rescount[[tf]][md,c(1,2,3,4)]<<-c(md,iconst,rconst,tf)
        }
        NULL
      })
      if(verbose) setTxtProgressBar(pb, i/length(modulators))
    }
    if(verbose)close(pb)
    #set data format
    sapply(names(rescount),function(tf){
      results<-rescount[[tf]][-1,,drop=FALSE]
      if(nrow(results)>0){
        results[,"AdjustedPvalue"] <- p.adjust(results[,"Pvalue"], method=pAdjustMethod)
        results[,"Expected"]<-round(results[,"Expected"],2)
        results[,"Pvalue"]<-signif(results[,"Pvalue"], digits=4)
        results[,"AdjustedPvalue"]<-signif(results[,"AdjustedPvalue"], digits=4)
        results<-results[sort.list(results[,"Observed"], decreasing=TRUE),,drop=FALSE]
      }
      rescount[[tf]]<<-results
      NULL
    })
    ##update summary and return results
    object@results$conditional$count<-rescount
    object@results$conditional$effect<-reseffect
    object@status["Conditional"] <- "[x]"
    if(verbose)cat("-Conditional analysis complete! \n\n")
    return(object)
  }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "TNI",
  function(object) {
    cat("A TNI (transcriptional network inference) object:\n")
    message("--status:")
    print(tni.get(object, what=c("status")), quote=FALSE)
  }
)

