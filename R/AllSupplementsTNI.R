

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##       OPTIMIZED MI FUNCTIONS FOR LARGE-SCALE TNETS
##      -- IMPLEMENTS AN INTERFACE FOR ARACNE/MINET --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##This function takes a gene expression matrix (x), a list of TFs
##and computes a partial mutal information matrix, i.e., between each TF 
##and all potential targets. Sig. mi values are inferred by permutation 
##analysis.
tni.pmim<-function(x,tfs,estimator="pearson",setdg=0,getadj=FALSE){
  x=t(x)
  dg=setdg
  pmim=cor(x[,tfs],x, method=estimator,use="complete.obs")^2
  if(length(tfs)>1){
    diag(pmim[,tfs])=dg
  } else {
    pmim[,tfs]=dg
  }
  maxi=0.999999
  pmim[which(pmim > maxi)]=maxi
  pmim = -0.5 * log(1 - pmim)
  pmim[pmim<0]=0
  if(getadj){
    mim=matrix(dg,ncol=ncol(x),nrow=ncol(x))
    colnames(mim)=colnames(x)
    rownames(mim)=colnames(x)  
    mim[tfs,]=pmim
    mim[,tfs]=t(pmim)
    mim
  } else {
    pmim<-t(pmim)
    colnames(pmim)<-tfs
    pmim
  }
}

##This function takes a gene expression matrix (x), a tnet and
##computes a simple correlation matrix, i.e., between each TF 
##and all targets, and returns the mode of action (-/+)
##of all pre-computed sig. associations.
tni.cor<-function(x,tnet,estimator="pearson",dg=0, asInteger=TRUE, mapAssignedAssociation=TRUE){
  tfs<-colnames(tnet)
  tar<-rownames(tnet)
  ids<-unique(c(tfs,setdiff(tar,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,tar], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  if(length(tfs)>1)diag(pcorm[,tfs])=dg
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  if(mapAssignedAssociation)pcorm[tnet==0]=0
  pcorm
}
tni.tfmdcor<-function(x,tfs, mds, estimator="pearson",dg=0, asInteger=FALSE){
  ids<-unique(c(tfs,setdiff(mds,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,mds], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  pcorm
}

##------------------------------------------------------------------------
##This function takes a partial mutal information matrix and apply the dpi 
##filter (aracne algorithm). TFs should be in cols and potential targets 
##in rows.
tni.dpi<-function(pmim,eps=0){
  tfs<-match(colnames(pmim),rownames(pmim))
  tar<-setdiff(1:nrow(pmim),tfs)
  x=pmim
  for(i in tar){
    idx<-pmim[i,]>0
    if(sum(idx)>1){
      mi<-pmim[c(i,tfs[idx]),idx]
      mi<-cbind(c(0,mi[1,]),mi)
      mi<-aracne(mi,eps=eps)
      tarvec<-mi[-1,1]
      x[i,idx]<-tarvec
    }
  }
  if(ncol(pmim)>1)x[tfs,]<-aracne(pmim[tfs,],eps=eps)
  return(x)
}

##------------------------------------------------------------------------
##permutation analysis
tni.perm<-function(object,verbose=TRUE){
  ##local pmim function for permutation
  perm.pmim<-function(x,tf,estimator="pearson",nPerm=1000){
    x=t(x)
    tf<-which(tf==colnames(x))
    x=cor(replicate(nPerm,sample(x[,tf])),x[,-tf], method=estimator,use="complete.obs")^2
    maxi=0.999999
    x[which(x > maxi)]=maxi
    x = -0.5 * log(1 - x)
    x[x<0]=0
    t(x)
  }
  ##compute partial mi matrix
  pmim<-tni.pmim(object@gexp,object@transcriptionFactors,object@para$perm$estimator)
  #compute null distributions
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
  if( b1 && b2 && length(object@transcriptionFactors)>1){
    if(verbose)cat("-Performing permutation analysis (parallel version)...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    mipval<-parSapply(getOption("cluster"), object@transcriptionFactors, function(tf) {
      midist<-perm.pmim(object@gexp, tf, object@para$perm$estimator, object@para$perm$nPermutations)
      midist<-sort(midist)
      np<-length(midist)
      midist<-np-findInterval(pmim[,tf],midist)
      ##pseudocounts are added to avoid P-values of zero
      midist <- (1 + midist)/(1 + np)
      ##pvalue adjustment (local distributions)
      midist <- p.adjust(midist,method=object@para$perm$pAdjustMethod)
      midist    
    })
    return(mipval)
  } else {
    if(verbose)cat("-Performing permutation analysis...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    if(verbose) pb <- txtProgressBar(style=3)
    mipval<-sapply(1:length(object@transcriptionFactors), function(i){
      tf=object@transcriptionFactors[i]
      midist<-perm.pmim(object@gexp, tf, object@para$perm$estimator, object@para$perm$nPermutations)
      midist<-sort(midist)
      np<-length(midist)
      midist<-np-findInterval(pmim[,tf],midist)
      #midist <- np-getprobs(pmim[,tf],midist) 
      ##pseudocounts are added to avoid P-values of zero
      midist <- (1 + midist)/(1 + np)
      ##pvalue adjustment (local distributions)
      midist <- p.adjust(midist,method=object@para$perm$pAdjustMethod)
      if(verbose) setTxtProgressBar(pb, i/length(object@transcriptionFactors))      
      midist   
    })
    if(verbose)close(pb)
  }
  ##if globalAdjustment...
  if(object@para$perm$globalAdjustment){
    mipval<-p.adjust(mipval,method=object@para$perm$pAdjustMethod)
    mipval<-matrix(mipval,nrow=nrow(pmim),ncol=ncol(pmim))    
  }
  ##decide on the significance and return results
  pmim[mipval>object@para$perm$pValueCutoff]<-0
  return(list(tn.ref=pmim,adjpv=mipval))
}

##------------------------------------------------------------------------
##permutation analysis with pooled null distribution
tni.perm.pooled=function(object, parChunks=10, verbose=TRUE){
  ##local pmim function for permutation
  perm.pmim<-function(x,n, estimator="pearson"){
    x<-matrix(sample(x),ncol=nrow(x),nrow=ncol(x))
    x<-cor(x[,1:n],x, method=estimator,use="complete.obs")^2
    if(n>1){
      diag(x[,1:n])=NA
    } else {
      x[,1:n]=NA
    }
    maxi=0.999999
    x[which(x>maxi)]=maxi
    x = -0.5*log(1-x)
    x[x<0]=0
    t(x)
  }
  ##compute partial mi matrix and get unique values
  uniqueVec<-tni.pmim(object@gexp,object@transcriptionFactors,object@para$perm$estimator)
  uniqueVec<-sort(unique(as.numeric(uniqueVec)))
  ##initialize permutation count
  ctsum<-numeric(length(uniqueVec))
  ##build the null distribution via permutation  
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
  if( b1 && b2 && length(object@transcriptionFactors)>1) {
    if(verbose)cat("-Performing permutation analysis (parallel version)...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    if(object@para$perm$nPermutations<=parChunks){
      stop("NOTE: 'nPermutations' should be multiple and greater than 'parChunks'!")
    }
    nper<-object@para$perm$nPermutations/parChunks
    if(nper>as.integer(nper)){
      stop("NOTE: 'nPermutations' should be multiple of 'parChunks'!")
    }
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:parChunks){
      permdist<-parLapply(getOption("cluster"),1:nper,function(j){
        permt<-perm.pmim(object@gexp,length(object@transcriptionFactors),object@para$perm$estimator)
        permt<-sort(permt)
        length(permt)-findInterval(uniqueVec,permt)
        #length(permt)-getprobs(uniqueVec,permt)
      })
      sapply(1:length(permdist),function(j){
        ctsum<<-ctsum+permdist[[j]]
        NULL
      })
      rm(permdist);gc()
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose)cat("-Performing permutation analysis...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    if(verbose)pb<-txtProgressBar(style=3)
    sapply(1:object@para$perm$nPermutations, function(i){
      permt<-perm.pmim(object@gexp,length(object@transcriptionFactors),object@para$perm$estimator)
      permt<-sort(permt)
      permt<-length(permt)-findInterval(uniqueVec,permt)
      ctsum<<-ctsum+permt
      if(verbose)setTxtProgressBar(pb, i/object@para$perm$nPermutations)
      NULL
    })
  }
  if(verbose)close(pb)
  ##compute pvals
  pmim<-tni.pmim(object@gexp,object@transcriptionFactors,object@para$perm$estimator)
  np <- object@para$perm$nPermutations * ( prod(dim(pmim)) - length(object@transcriptionFactors) )
  mipval <- (1 + ctsum)/(1 + np)
  mipval<-mipval[match(as.numeric(pmim),uniqueVec)]
  ##adjust pvals
  mipval<-p.adjust(mipval,method=object@para$perm$pAdjustMethod)
  mipval<-matrix(mipval,nrow=nrow(pmim),ncol=ncol(pmim))
  ##decide on the significance
  pmim[mipval>object@para$perm$pValueCutoff]=0.0
  return(list(tn.ref=pmim,adjpv=mipval))
}

##------------------------------------------------------------------------
##compute delta mi from pooled null distributions, returns global mi markers
miThresholdMd=function(gexp, nsamples, prob=0.95, nPermutations=1000, 
                       estimator="pearson", verbose=TRUE){
  perm.pmim<-function(x, nsamples, nsc, estimator="pearson"){
    x<-t(x[,sample(ncol(x),nsamples)]) ##i.e. sample modulators
    x<-cor(x[,sample(ncol(x),nsc)],x, method=estimator,use="complete.obs")^2
    diag(x[,1:nsc])=0
    maxi=0.999999
    x[which(x>maxi)]=maxi
    x = -0.5*log(1-x)
    x[x<0]=0
    t(x)
  }
  if(nsamples>ncol(gexp))nsamples=ncol(gexp)
  ##build null distribution
  nsc<-10 #set a scale-up factor for permutation
  if(verbose)pb<-txtProgressBar(style=3)
  nullmark<-sapply(1:nPermutations, function(i){
    if(verbose)setTxtProgressBar(pb, i/nPermutations)
    rmil<-perm.pmim(gexp,nsamples=nsamples,nsc=nsc,estimator=estimator)
    rmih<-perm.pmim(gexp,nsamples=nsamples,nsc=nsc,estimator=estimator)
    permt<-abs(rmih-rmil)
    quantile(permt,probs=prob,na.rm=TRUE,names=FALSE)
  })
  if(verbose)close(pb)
  mimark<-max(nullmark,na.rm=TRUE)
  return(mimark)
}

##------------------------------------------------------------------------
##compute delta mi, returns mi markers for each tf individually
miThresholdMdTf<-function(gexp, tfs, nsamples, prob=0.95, nPermutations=1000, 
                          estimator="pearson", verbose=TRUE){
  #---get low/high sample idxs
  nc<-ncol(gexp)
  idxl<-1:nsamples
  idxh<-(nc-nsamples-1):nc
  #---run permutation
  if(verbose)pb<-txtProgressBar(style=3)
  rmid<-sapply(1:nPermutations,function(i){
    if(verbose)setTxtProgressBar(pb, i/nPermutations)
    ridx<-sample(nc)
    rmil<-tni.pmim(gexp[,ridx[idxl]],tfs=tfs,estimator=estimator)
    rmih<-tni.pmim(gexp[,ridx[idxh]],tfs=tfs,estimator=estimator)
    rmid<-abs(rmih-rmil)
    apply(rmid,2,quantile,probs=prob,na.rm=TRUE,names=FALSE)
  })
  if(verbose)close(pb)
  if(length(tfs)>1){
    dmimark<-apply(rmid,1,max,na.rm=TRUE)
  } else {
    dmimark<-max(rmid,na.rm=TRUE)
  }
  names(dmimark)<-tfs
  dmimark
}

##compute delta mi, returns mi markers for each tf-target individually
miThresholdMdTfTar<-function(gexp, tfs, nsamples, prob=0.95, nPermutations=1000, 
                             estimator="pearson", verbose=TRUE){
  #---get low/high sample idxs
  nc<-ncol(gexp)
  idxl<-1:nsamples
  idxh<-(nc-nsamples-1):nc
  #---run permutation
  dmimark<-sapply(1:length(tfs),function(i){
    if(verbose)cat(paste("-",i,"/",length(tfs),"\n",sep=""))
    if(verbose)pb<-txtProgressBar(style=3)
    rmid<-sapply(1:nPermutations,function(j){
      if(verbose)setTxtProgressBar(pb, j/nPermutations)
      ridx<-sample(nc)
      rmil<-tni.pmim(gexp[,ridx[idxl]],tfs=tfs[i],estimator=estimator)
      rmih<-tni.pmim(gexp[,ridx[idxh]],tfs=tfs[i],estimator=estimator)
      abs(rmih-rmil)
    })
    if(verbose)close(pb)
    apply(rmid,1,quantile,probs=prob,na.rm=TRUE,names=FALSE)
  })
  rownames(dmimark)<-rownames(gexp)
  colnames(dmimark)<-tfs
  dmimark
}

##------------------------------------------------------------------------
#bootstrap analysis
tni.boot<-function(object, parChunks=10, verbose=TRUE){
  #local pmim function
  boot.pmim<-function(x,tfs,estimator="pearson"){
    x<-t(x[,sample(1:ncol(x),replace=TRUE)])
    x<-cor(x[,tfs],x, method=estimator,use="complete.obs")^2
    if(length(tfs)>1){
      diag(x[,tfs])=0
    } else {
      x[,tfs]=0
    }
    maxi=0.999999
    x[which(x>maxi)]=maxi
    x <- -0.5*log(1-x)
    x[x<0]=0
    t(x)
  }  
  ##identify empirical boundaries from tn.ref;
  ##..this strategy turns the bootstrap more tractable for large 
  ##..scale datasets.
  mimark<-sapply(1:ncol(object@results$tn.ref), function(i){
    tp<-abs(object@results$tn.ref[object@results$tn.ref[,i]!=0,i])
    ifelse(length(tp)>0 , min(tp), 0)
  })
  if(object@para$perm$globalAdjustment){
    mimark<-rep(max(mimark),length(mimark))
  }
  ##initialize bootstrap matrix
  bcount<-matrix(0,ncol=ncol(object@results$tn.ref),nrow=nrow(object@results$tn.ref))
  ##run bootstrap
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
  if( b1 && b2 && length(object@transcriptionFactors)>1) {
    if(verbose)cat("-Performing bootstrap analysis (parallel version)...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    if(object@para$boot$nBootstraps<=parChunks){
      stop("NOTE: 'nBootstraps' should be multiple and greater than 'parChunks'!")
    }
    nboot<-object@para$boot$nBootstraps/parChunks
    if(nboot>as.integer(nboot)){
      stop("NOTE: 'nBootstraps' should be multiple of'parChunks'!")
    }
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:parChunks){
      bootdist<-parLapply(getOption("cluster"),1:nboot,function(j){
        boott<-boot.pmim(object@gexp,object@transcriptionFactors)
        sapply(1:ncol(boott),function(k){as.numeric(boott[,k]>mimark[k])})
      })
      sapply(1:length(bootdist),function(j){
        bcount<<-bcount+bootdist[[j]]
        NULL
      })
      rm(bootdist);gc()
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose)cat("-Performing bootstrap analysis...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons...\n")
    if(verbose) pb <- txtProgressBar(style=3)    
    sapply(1:object@para$boot$nBootstraps,function(i){
      if(verbose) setTxtProgressBar(pb, i/object@para$boot$nBootstraps)
      bootdist<-boot.pmim(object@gexp,object@transcriptionFactors)
      bootdist<-sapply(1:ncol(bootdist),function(j){as.numeric(bootdist[,j]>mimark[j])})
      bcount <<- bcount + bootdist
      NULL
    })
  }
  if(verbose)close(pb)
  ##decide on the stability of the inferred associations
  cons<-as.integer( object@para$boot$nBootstraps * (object@para$boot$consensus/100) )
  object@results$tn.ref[bcount<=cons]<-0
  return(object@results$tn.ref)
}

##------------------------------------------------------------------------
##This function takes a row gene expression matrix, a named vector with
##matched ids and remove duplicated genes using the coefficient
##of variation (CV) to decide for the most informative probes.
cv.filter<-function(gexp, ids){ 
  # compute CV
  meangx<-apply(gexp,1,mean,na.rm=TRUE)
  sdgx<-apply(gexp,1,sd,na.rm=TRUE)
  cvgx<-sdgx/meangx
  # aline ids
  ids <- data.frame(ids[rownames(gexp),], CV=cvgx, stringsAsFactors=FALSE)
  # remove probes without annotation
  ids <- ids[!is.na(ids[,2]),]
  ids <- ids[ids[,2]!="",]
  # remove duplicated genes by CV
  ids <- ids[sort.list(ids$CV,decreasing=TRUE),]
  unids<-unique(ids[,2])
  idx<-match(unids,ids[,2])
  ids<-ids[idx,]
  # update and return gexp matrix
  gexp<-gexp[rownames(ids),]
  ids<-ids[,-ncol(ids)]
  return(list(gexp=gexp,ids=ids))
}

##------------------------------------------------------------------------
##---Simplified version of cv.filter---
##This function takes a gene expression matrix and remove duplicated genes 
##using the coefficient variation (CV) to decide for the most informative probes.
##Input format -> rownames: probeID; col[1]: geneID, Col[2...n]:numeric data
cv.filter.simp<-function(gexp){
  x<-as.matrix(gexp[,-1])
  ids<-data.frame(ProbeID=rownames(gexp),geneID=gexp[,1],stringsAsFactors=FALSE)
  rownames(ids)<-ids$ID1
  # compute CV
  meangx<-apply(x,1,mean,na.rm=TRUE)
  sdgx<-apply(x,1,sd,na.rm=TRUE)
  cvgx<-sdgx/meangx
  # aline ids
  ids <- data.frame(ids, CV=cvgx, stringsAsFactors=FALSE)
  # remove proves without annotation
  ids <- ids[!is.na(ids[,2]),]
  ids <- ids[ids[,2]!="",]
  # remove duplicated genes by CV
  ids <- ids[sort.list(ids$CV,decreasing=TRUE),]
  unids<-unique(ids[,2])
  idx<-match(unids,ids[,2])
  ids<-ids[idx,]
  # update and return gexp matrix
  x<-x[rownames(ids),]
  ids<-ids[,-ncol(ids)]
  rownames(x)<-ids$geneID
  #return(list(gexp=x,ids=ids))
  x
}


##------------------------------------------------------------------------
##------------------------------------------------------------------------
##              FUNCTIONS FOR REGULONS AND TNETS
##          -- IMPLEMENTS METHODS FOR REDER/IGRAPH --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##return an adjacence matrix with regulons (tfs and targets)
tni.rmap<-function(tnet){
  tnet[tnet>0]=1;tnet[tnet<0]=-1
  tfs<-colnames(tnet)
  if(sum(!tfs%in%rownames(tnet))>0){
    idx<-!tfs%in%rownames(tnet)
    addtfs<-matrix(0,nrow=sum(idx),ncol=ncol(tnet))
    rownames(addtfs)<-tfs[idx]
    tnet<-rbind(addtfs,tnet)
  }
  idx<-apply(tnet!=0,1,sum)>0
  onlytar<-rownames(tnet)[idx]
  onlytar<-onlytar[!onlytar%in%tfs]
  tnet<-tnet[c(tfs,onlytar),,drop=FALSE]
  tmat<-t(tnet[onlytar,,drop=FALSE])
  zmat<-matrix(0,length(onlytar),length(onlytar))
  bmat<-rbind(tmat,zmat)
  rmap<-cbind(tnet,bmat)
  rmap[onlytar,tfs]<-0
  g<-igraph::graph.adjacency(rmap, diag=FALSE, mode="directed", weighted=TRUE)
  edl<-igraph::get.edgelist(g)
  if(nrow(edl)>0){
    etype<-edl[,1]%in%tfs+edl[,2]%in%tfs
    etype[etype==2]<-0
    E(g)$modeOfAction<-E(g)$weight
    E(g)$modeOfAction[etype==0]=0
    E(g)$arrowType<-E(g)$modeOfAction
  }
  return(g)
}

##-----------------------------------------------------------------------------
##return an adjacence matrix for regulons
##using jaccard coefficient
tni.amap<-function(tnet){
  tnet[tnet!=0]=1
  jc<-function(x,xmat){
    c<-x+xmat
    a<-colSums(c==2)
    b<-colSums(c>0)
    b[b==0]=1
    a/b
  }
  amap<-apply(tnet,2,jc,xmat=tnet)
  colnames(amap)<-colnames(tnet)
  diag(amap)=0
  return(amap)
}

##-----------------------------------------------------------------------------
##return an adjacence matrix
tni.mmap<-function(object,mnet,tnet,pvnet,othertfs,testedtfs,modulators){
  tnet[tnet<0]=-1;tnet[tnet>0]=1
  mnet<-mnet[,colnames(tnet),drop=FALSE]
  mnet<-mnet[rownames(tnet),,drop=FALSE]
  pvnet<-pvnet[,colnames(tnet),drop=FALSE]
  pvnet<-pvnet[rownames(tnet),,drop=FALSE]
  #---merge/simplify tnet and mnet
  mcnet<-mnet*10
  mcnet[mcnet==0]<-tnet[mcnet==0]
  #---get ids
  tfs<-colnames(mcnet)
  mod<-rownames(mcnet)
  onlytar<-setdiff(mod,tfs)
  #---check tfs
  if(sum(!tfs%in%rownames(mcnet))>0){
    idx<-!tfs%in%rownames(mcnet)
    addtfs<-matrix(0,nrow=sum(idx),ncol=ncol(mcnet))
    rownames(addtfs)<-tfs[idx]
    mcnet<-rbind(addtfs,mcnet)
    pvnet<-rbind(addtfs,pvnet)
  }
  #---set ordering
  mcnet<-mcnet[c(tfs,onlytar),,drop=FALSE]
  pvnet<-pvnet[c(tfs,onlytar),,drop=FALSE]
  #---
  elist<-NULL
  emode<-NULL
  etype<-NULL
  pvalue<-NULL
  tpnet<-mcnet[tfs,tfs,drop=FALSE]
  tpnet[tpnet!=0]=1
  sapply(1:nrow(tpnet),function(i){
    sapply(1:ncol(tpnet),function(j){
      tp<-tpnet[i,j] + tpnet[j,i]
      if(j>i){
        if(tp!=0){
          emode<<-c(emode,0)
          etype<<-c(etype,"TF.TF")
          elist<<-rbind(elist,c(tfs[i],tfs[j]))
          pvalue<<-c(pvalue,NA)
          emode<<-c(emode,0)
          etype<<-c(etype,"TF.TF")
          elist<<-rbind(elist,c(tfs[j],tfs[i]))
          pvalue<<-c(pvalue,NA)
        }
      }
      NULL
    })
    NULL
  })
  #others
  if(length(onlytar)>0){
    tpnet<-mcnet[onlytar,tfs,drop=FALSE]
    tpvnet<-pvnet[onlytar,tfs,drop=FALSE]
    sapply(1:nrow(tpnet),function(i){
      sapply(1:ncol(tpnet),function(j){
        tp<-tpnet[i,j]
        tpv<-tpvnet[i,j]
        if(abs(tp)==10){
          emode<<-c(emode,tp)
          tpp<-ifelse(tp<0,"-Md.TF","+Md.TF")
          etype<<-c(etype,tpp)
          elist<<-rbind(elist,c(onlytar[i],tfs[j]))
          pvalue<<-c(pvalue,tpv)
        } else if(abs(tp)==1){
          emode<<-c(emode,tp)
          tpp<-ifelse(tp<0,"-TF.Target","+TF.Target")
          etype<<-c(etype,tpp)
          elist<<-rbind(elist,c(tfs[j],onlytar[i]))
          pvalue<<-c(pvalue,tpv)
        }
        NULL
      })
      NULL
    })
  }
  #Correct TF-TF assigments derived from the MI analysis
  if(length(othertfs)>0){
    alltfs<-unique(c(tfs,othertfs))
    idx<-elist[,1]%in%alltfs+elist[,2]%in%alltfs == 2
    idx<-idx & abs(emode)==1
    emode[idx]<-0
    etype[idx]<-"TF.TF"
  }
  if(!is.null(emode)){
    idx<-abs(emode)==10
    emode[idx]<-emode[idx]/10
  }
  #---get graph and set mode
  if(is.null(elist)){
    g<-graph.empty(n=0, directed=TRUE)
  } else {
    g<-igraph::graph.edgelist(elist, directed=TRUE)
    E(g)$arrowType<-emode
    E(g)$pvalueModulation<-pvalue
    E(g)$interactionType<-etype
    E(g)$modeOfAction<-NA
    E(g)$modeOfAction[E(g)$interactionType=="-Md.TF"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+Md.TF"]<- 1
    E(g)$modeOfAction[E(g)$interactionType=="TF.TF"]<- 0
    E(g)$modeOfAction[E(g)$interactionType=="-TF.Target"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+TF.Target"]<- 1
  }
  #---add disconnected tfs 
  idx<-!tfs%in%V(g)$name
  if(any(idx))g<-g+vertices(tfs[idx])
  #---
  g<-setregs1(object,g,testedtfs,modulators)
  return(g)
}
setregs1<-function(object,g,testedtfs,modulators){
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
  #---interno, so p/ customizar -- linearSize!!
  customSz=TRUE
  if(customSz){
    nquant=NULL
    tp1<-min(V(g)$regulonSize,na.rm=TRUE)
    tp3<-max(V(g)$regulonSize,na.rm=TRUE)
    tp1<-pretty(tp1,n=1,shrink=0.5)[1]
    tp3<-pretty(tp3,n=1,shrink=0.5)[2]
    tp2<-round((tp1+tp3)/2,0)
    breaks<-c(tp1,tp2,tp3)
  } else {
    breaks=NULL
    nq<-unique(V(g)$regulonSize)
    nq<-sum(!is.na(nq))
    nquant<-ifelse(nq>5,5,nq)
    nq<-length(unique(quantile(V(g)$regulonSize,na.rm=TRUE)))
    if(nq<3)nquant=NULL
  }
  xlim = c(25, 70, 10)
  #---
  g<-att.setv(g=g, from="regulonSize", to='nodeSize',xlim=xlim,roundleg=0,
              nquant=nquant,title="Downstream targets (n)",breaks=breaks)
  #set node attr
  V(g)$nodeFontSize=15
  V(g)$nodeLineWidth<-1.5
  V(g)$nodeColor<-"grey"
  V(g)$nodeLineColor<-"grey"
  V(g)$nodeFontSize[V(g)$tfs==1]=25
  #set node shape and legend
  g<-att.setv(g=g, from="tfs", to='nodeShape',shapes=c("DIAMOND","ELLIPSE"),title="", categvec=c(0,1))
  g$legNodeShape$legend<-c("Modulator","TF")
  if(ecount(g)>0){
    #set edge attr
    E(g)$edgeWidth<-1.5
    if(any(E(g)$modeOfAction==0)){
      g<-att.sete(g=g, from="modeOfAction", to='edgeColor',cols=c("#96D1FF","grey80","#FF8E91"), 
                  title="Mode",categvec=c(-1,0,1))
      g$legEdgeColor$legend<-c("Down","NA","Up")
    } else {
      g<-att.sete(g=g, from="modeOfAction", to='edgeColor',cols=c("#96D1FF","grey80","#FF8E91"), 
                  title="Mode",categvec=c(-1,1))
      g$legEdgeColor$legend<-c("Down","Up")
    }
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
    V(g)$nodeColor[V(g)$tfs==1]<-"white"
    V(g)$nodeLineColor[V(g)$tfs==1]<-"black"
  }
  g
}

##-----------------------------------------------------------------------------
##experimental code!!! (opem mmaps)
tni.mmap.detailed<-function(object,mnet,testedtfs,ntop){
  mnet.targets<-list()
  junk<-sapply(testedtfs,function(tf){
    tp1<-object@results$conditional$effect[[tf]]
    tp2<-mnet[,tf];tp2<-tp2[tp2!=0]
    tp3<-tp1[,names(tp2),drop=FALSE]
    if(ncol(tp3)>0){
      res<-sapply(names(tp2),function(md){
        tp<-tp3[,md]
        if(tp2[md]>0) tp[tp<0]<-0 else tp[tp>0]<-0
        tp
      })
      rownames(res)<-tp1$tagets
      mnet.targets[[tf]]<<-res
    }
    NULL
  })
  #---
  gLists<-lapply(names(mnet.targets),function(tf){
    tar<-object@results$tn.dpi[,tf]
    mtar<-mnet.targets[[tf]]
    mds<-colnames(mtar)
    glists<-lapply(mds,function(md){
      mtartp<-mtar[,md]
      mode<-mnet[md,tf]
      tartp<-tar[names(mtartp)]
      elist<-cbind(A=c(md,tf,rep("TF|Md",length(tartp))),B=c("TF|Md","TF|Md",names(tartp)))
      rownames(elist)<-NULL
      #---
      modeOfAction<-c(NA,NA,tartp)
      modulationEffect<-c(mode,1,mtartp)
      modulationType<-modulationEffect
      modulationType[modulationType!=0]<-1
      modulationType[1:2]<-2
      #---
      interactionType<-rep("TF.Target",length(modeOfAction))
      interactionType[1:2]<-"TF|Md"
      elist<-data.frame(elist,interactionType,modeOfAction,modulationEffect,modulationType,stringsAsFactors=FALSE)
      idx<-sort.list(abs(elist$modeOfAction),na.last=FALSE,decreasing=FALSE)
      elist<-elist[idx,]
      if(!is.null(ntop)){
        idx<-c(1,2,rev(3:nrow(elist)))[1:(ntop+2)]
        elist<-elist[sort(idx),]
      }
      #---
      g<-igraph::graph.edgelist(as.matrix(elist[,1:2]), directed=TRUE)
      #set node type
      V(g)$nodeType<-rep("Target",length(V(g)$name))
      idx<-match(c(md,"TF|Md",tf),V(g)$name)
      V(g)$nodeType[idx]<-c("Md","TF|Md","TF")
      #---
      E(g)$interactionType<-elist$interactionType
      E(g)$modeOfAction<-elist$modeOfAction
      E(g)$modulationEffect<-elist$modulationEffect
      E(g)$modulationType<-elist$modulationType
      xy<-drcircle(nv=vcount(g)-3)
      hd<-data.frame(x=c(-0.14,0,0.2),y=c(-0.9,-0.7,-0.8))
      xy<-rbind(hd,xy)
      V(g)$coordX<-xy$x
      V(g)$coordY<-xy$y
      g
    })
    names(glists)<-mds
    glists
  })
  names(gLists)<-names(mnet.targets)
  #---
  gLists<-lapply(gLists,function(gg){
    gg<-lapply(gg,function(g){
      if(.hasSlot(object, "annotation") && nrow(object@annotation)>0){
        g<-att.mapv(g=g,dat=object@annotation,refcol=1)
      } else {
        V(g)$SYMBOL<-V(g)$name
      }
      if(ecount(g)>0)g<-setregs2(g)
    })
    gg
  })
  gLists
}
setregs2<-function(g){
  g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
  V(g)$nodeAlias[V(g)$name=="TF|Md"]<-"TF | Md"
  V(g)$nodeColor<-"white"
  g<-att.setv(g=g, from="nodeType", to='nodeLineColor',categvec=c("Md","TF","TF|Md","Target"),
              cols=c("#FF6666","#66CCFF","#FF9900","#00CC99"),title="")
  g$legNodeLineColor$legend<-c("Md","TF","TF | Md", "Mod. target")
  #--set legNodeColor
  g$legNodeColor<-g$legNodeLineColor
  g<-remove.graph.attribute(g,"legNodeLineColor")
  #--
  V(g)$nodeSize<-30
  V(g)$nodeSize[V(g)$nodeType=="Target"]<-20
  V(g)$nodeFontSize<-12
  V(g)$nodeFontSize[V(g)$nodeType!="Target"]<-16
  V(g)$nodeLineWidth<-2.8
  #set edge attr
  E(g)$edgeColor<-"grey85"
  E(g)$edgeWidth<-2
  E(g)$edgeWidth[E(g)$modulationType==0]<-1
  #--
  g<-att.sete(g=g, from="modulationType", to='edgeColor', title="",
              categvec=c(2,1,0),cols=c("black","grey65","grey85"))
  g$legEdgeColor$legend<-c("conditioned","modulated","non-modulated")
  #--
  el<-data.frame(get.edgelist(g),E(g)$modulationType,stringsAsFactors=FALSE)
  ndm<-el[el[,3]==0,2]
  V(g)$nodeLineColor[V(g)$name%in%ndm]<-"grey85"
  V(g)$nodeFontColor<-"black"
  V(g)$nodeFontColor[V(g)$name%in%ndm]<-"grey85"
  
  
  
  
  g
}
drcircle<-function (nv, ybend=1, xbend=1, ang=1, rotate=0, plot=FALSE) {
  angle.inc <- ang * pi/(nv-1)
  angles <- seq(0, ang*pi, by=angle.inc)+rotate*pi
  radius<-0.5
  xv <- cos(angles) * radius * xbend
  yv <- sin(angles) * radius * ybend
  if(plot){
    plot(-1:1,-1:1,type="n",xlab="",ylab="",main="")
    polygon(xv, yv)
  }
  data.frame(x=xv,y=yv)
}
#drcircle(nv=10,ang=1,plot=TRUE)

#---------------------------------------------------------------
#test overlap of among regulons in a tnet via phyper function
tni.phyper<-function(tnet){
  tnet[tnet!=0]<-1
  tnet<-tnet[rowSums(tnet)>0,]
  phtest<-function(x,xmat){
    c<-x+xmat
    pop<-length(x)
    sample1<-sum(x>0)
    sample2<-colSums(xmat>0)
    overlap<-colSums(c==2)
    resph<-phyper(overlap, sample1, pop - sample1, sample2, lower.tail=FALSE)
    resph
  }
  pmat<-apply(tnet,2,phtest,xmat=tnet)
  colnames(pmat)<-colnames(tnet)
  diag(pmat)=1
  pmat<-matrix(p.adjust(pmat),ncol=ncol(pmat),nrow=nrow(pmat),
               dimnames=list(rownames(pmat),colnames(pmat)))
  return(pmat)
}

#---------------------------------------------------------------
#reverse results from conditional analyses, from TF-MD to MD-TF
cdtReverse<-function(cdt,pAdjustMethod="bonferroni"){
  cdtrev<-list()
  sapply(1:length(cdt),function(i){
    tf<-names(cdt)[i]
    tp<-cdt[[i]]
    if(nrow(tp)>0){
      mds<-rownames(tp)
      sapply(mds,function(md){
        tpp<-tp[md,c(2,1,3:ncol(tp))]
        rownames(tpp)<-tf
        cdtrev[[md]]<<-rbind(cdtrev[[md]],tpp)
      })
    }
  })
  if(length(cdtrev)>0){
    cdtrev<-p.adjust.cdt(cdt=cdtrev,pAdjustMethod=pAdjustMethod, p.name="PvFET",adjp.name="AdjPvFET",sort.name="PvKS")
    cdtrev<-p.adjust.cdt(cdt=cdtrev,pAdjustMethod=pAdjustMethod, p.name="PvFET",adjp.name="AdjPvFET")
  }
  cdtrev
}
#---------------------------------------------------------------
#compute global p.adjustment for conditional analysis
p.adjust.cdt<-function(cdt,pAdjustMethod="bonferroni",p.name="Pvalue",
                       adjp.name="AdjustedPvalue",sort.name=NULL, decreasing=FALSE){
  adjp<-array(,dim=c(1,2))
  sapply(1:length(cdt),function(i){
    tp<-cdt[[i]]
    if(nrow(tp)>0){
      if(is.character(sort.name)){
        tp<-tp[sort.list(tp[[sort.name]],decreasing=decreasing),]
        cdt[[i]]<<-tp
      }
      adjp<<-rbind(adjp,cbind(i,tp[[p.name]]))
    }
    NULL
  })
  adjp<-adjp[-1,,drop=FALSE]
  colnames(adjp)<-c("idx","p")
  adjp[,"p"]<-p.adjust(adjp[,"p"],method=pAdjustMethod)
  adjp[,"p"]<-signif(adjp[,"p"], digits=4)
  #update adjp
  idx<-unique(adjp[,1])
  sapply(idx,function(i){
    tp<-cdt[[i]]
    tp[[adjp.name]]<-adjp[adjp[,1]==i,2]
    cdt[[i]]<<-tp
  })
  cdt
}
sortblock.cdt<-function(cdt,coln="PvFET"){
  #sort blocks
  idx1<-unlist(lapply(cdt,nrow))
  idx1<-idx1/max(idx1)+1
  idx2<-lapply(1:length(cdt),function(i){
    tp<-cdt[[i]][[coln]]
    if(length(tp)>0 && sum(tp,na.rm=TRUE)>0){
      tp<--log10(mean(tp,na.rm=TRUE))
    } else {
      tp<-0
    }
    tp
  })
  idx2<-unlist(idx2)
  idx2<-idx2/sum(idx2)
  idx<-sort.list(idx1+idx2,decreasing=TRUE)
  cdt<-cdt[idx]
  cdt
}
#---------------------------------------------------------------
#test Cohen's Kappa agreement between Modulon and regulons
# pkappa<-function (Modulon,Regulon){
#   ttab<-data.frame(A=Modulon,B=Regulon)
#   ttab<-table(ttab)
#   #---get tabs
#   nc <- ncol(ttab)
#   ns <- sum(ttab)
#   weighttab <- matrix(0, nrow = nc, ncol = nc)
#   diag(weighttab)<-c(1,1)
#   #---kvalue
#   agreeP <- sum(ttab * weighttab)/ns
#   tm1 <- apply(ttab, 1, sum)
#   tm2 <- apply(ttab, 2, sum)
#   eij <- outer(tm1, tm2)/ns
#   chanceP <- sum(eij * weighttab)/ns
#   kvalue <- (agreeP - chanceP)/(1 - chanceP)
#   return(kvalue)
# }
#test overlap of among regulons in a tnet via phyper function
# tni.phyper<-function(tnet){
#   tnet[tnet!=0]=1
#   tnet<-tnet[rowSums(tnet)>0,]
#   labs<-colnames(tnet)
#   pmat<-matrix(NA,ncol=ncol(tnet),nrow=ncol(tnet), dimnames=list(labs, labs))
#   phtest<-function(c1,c2){
#     c<-c1+c2
#     pop<-length(c)
#     sample1<-sum(c1>0)
#     sample2<-sum(c2>0)
#     overlap<-sum(c==2)  	
#     resph<-phyper(overlap, sample1, pop - sample1, sample2, lower.tail=FALSE)
#     resph
#   }
#   for(i in 1:ncol(tnet)){
#     v1<-tnet[,i]
#     for(j in 1:ncol(tnet)){
#       if(j>i){
#         v2<-tnet[,j]
#         ph<-phtest(as.integer(v1),as.integer(v2))
#         pmat[i,j]<-ph
#       }
#     }
#   }
#   lt<-lower.tri(pmat)
#   pmat[lt]=t(pmat)[lt]
#   diag(pmat)=1
#   pmat<-matrix(p.adjust(pmat),ncol=ncol(pmat),nrow=nrow(pmat),
#                dimnames=list(rownames(pmat),colnames(pmat)))           
#   pmat
# }

