

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
tni.cor<-function(x,tnet,estimator="pearson",dg=0, asInteger=TRUE){
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
  pcorm[tnet==0]=0
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
    if(verbose)cat("-Performing permutation analysis (parallel version - ProgressBar not available) ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
    mipval<-parSapply(getOption("cluster"), object@transcriptionFactors, function(tf) {
      midist<-perm.pmim(object@gexp, tf, object@para$perm$estimator, object@para$perm$nPermutations)
      midist<-sort(midist)
      np<-length(midist)
      midist<-np-findInterval(pmim[,tf],midist)
      #midist<-np-getprobs(pmim[,tf],midist)
      ##pseudocounts are added to avoid P-values of zero
      midist <- (1 + midist)/(1 + np)
      ##pvalue adjustment (local distributions)
      midist <- p.adjust(midist,method=object@para$perm$pAdjustMethod)
      midist    
    })
    return(mipval)
  } else {
    if(verbose)cat("-Performing permutation analysis ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
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
    if(verbose)cat("-Performing permutation analysis (parallel version) ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
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
    if(verbose)cat("-Performing permutation analysis ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
    if(verbose)pb<-txtProgressBar(style=3)
    sapply(1:object@para$perm$nPermutations, function(i){
      permt<-perm.pmim(object@gexp,length(object@transcriptionFactors),object@para$perm$estimator)
      permt<-sort(permt)
      permt<-length(permt)-findInterval(uniqueVec,permt)
      #permt<-length(permt)-getprobs(uniqueVec,permt)
      ctsum<<-ctsum+permt
      if(verbose)setTxtProgressBar(pb, i/object@para$perm$nPermutations)
      NULL
    })
  }
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
##permutation analysis with pooled null distribution used to infer
##statistical boundaries for the conditional mutual information analysis
tni.conditional.threshold=function(object, ntfs, nsamples, parChunks=10, verbose=TRUE){
  ##local pmim function for permutation
  perm.pmim<-function(x,ntfs, nsamples, estimator="pearson"){
    x<-x[,sample(1:ncol(x))[1:nsamples]]
    x<-matrix(sample(x),ncol=nrow(x),nrow=ncol(x))
    x<-cor(x[,1:ntfs],x, method=estimator,use="complete.obs")^2
    if(ntfs>1){
      diag(x[,1:ntfs])=0
    } else {
      x[,1:ntfs]=0
    }
    maxi=0.999999
    x[which(x>maxi)]=maxi
    x = -0.5*log(1-x)
    x[x<0]=0
    t(x)
  }
  #set ntfs to a minimun!
  #to reduce fluctuations in small jobs!
  if(ntfs<15)ntfs=15
  #set cutoffmarks
  nullmark=NULL
  nmarks= (nrow(object@gexp) * ntfs )
  tp<-object@para$cdt$pValueCutoff*1/2 #threshold position!
  idxmark<-as.integer(nmarks * tp)
  if(idxmark<2)idxmark=2
  idxmark=nmarks-idxmark
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
  if( b1 && b2) {
    if(verbose)cat("-Estimating mutual information threshold (parallel version) ...\n")
    if(object@para$perm$nPermutations<=parChunks){
      stop("NOTE: 'nPermutations' should be multiple and greater than 'parChunks'!")
    }
    nper<-object@para$perm$nPermutations/parChunks
    if(nper>as.integer(nper)){
      stop("NOTE: 'nPermutations' should be multiple of 'parChunks'!")
    }
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:parChunks){
      permmark<-parSapply(getOption("cluster"),1:nper,function(j){
        permt<-perm.pmim(object@gexp,ntfs=ntfs,nsamples=nsamples,estimator=object@para$perm$estimator)
        sort(permt)[idxmark]
      })
      nullmark<-c(nullmark,permmark)
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose)cat("-Estimating mutual information threshold ...\n")
    if(verbose)pb<-txtProgressBar(style=3)
    nullmark<-sapply(1:object@para$perm$nPermutations, function(i){
      if(verbose)setTxtProgressBar(pb, i/object@para$perm$nPermutations)
      permt<-perm.pmim(object@gexp,ntfs=ntfs,nsamples=nsamples,estimator=object@para$perm$estimator)
      sort(permt)[idxmark]
    })
  }
  ##get mi threshold
  mmark<-max(nullmark,na.rm=TRUE)
  return(mmark)
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
  mmark<-sapply(1:ncol(object@results$tn.ref), function(i){
    tp<-abs(object@results$tn.ref[object@results$tn.ref[,i]!=0,i])
    ifelse(length(tp)>0 , min(tp), 0)
  })
  if(object@para$perm$globalAdjustment){
    mmark<-rep(max(mmark),length(mmark))
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
    if(verbose)cat("-Performing bootstrap analysis (parallel version) ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
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
        sapply(1:ncol(boott),function(k){as.numeric(boott[,k]>mmark[k])})
      })
      sapply(1:length(bootdist),function(j){
        bcount<<-bcount+bootdist[[j]]
        NULL
      })
      rm(bootdist);gc()
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose)cat("-Performing bootstrap analysis ...\n")
    if(verbose)cat("--For", length(object@transcriptionFactors), "regulons \n")
    if(verbose) pb <- txtProgressBar(style=3)    
    sapply(1:object@para$boot$nBootstraps,function(i){
      if(verbose) setTxtProgressBar(pb, i/object@para$boot$nBootstraps)
      bootdist<-boot.pmim(object@gexp,object@transcriptionFactors)
      bootdist<-sapply(1:ncol(bootdist),function(j){as.numeric(bootdist[,j]>mmark[j])})
      bcount <<- bcount + bootdist
      NULL
    })
    if(verbose)close(pb)
  }
  ##decide on the stability of the inferred associations
  cons<-as.integer( object@para$boot$nBootstraps * (object@para$boot$consensus/100) )
  object@results$tn.ref[bcount<=cons]<-0
  return(object@results$tn.ref)
}

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##      OPTIMIZED GSEA FUNCTIONS FOR REGULONS AND TNETS
##       -- IMPLEMENTS AN INTERFACE FOR HTSanalyzeR --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##This function computes observed and permutation-based scores associated 
##with a gene set enrichment analysis for a collection of regulons.
gsea1tna <- function(listOfRegulons, phenotype, exponent=1, nPermutations=1000, 
                     minRegulonSize=15, verbose=TRUE) {	 
  ##convert names to integers if parallel  
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
  isPar <- b1 && b2
#   if(isPar){
#     intnames<-c(1:length(phenotype))
#     names(intnames)<-names(phenotype)
#     names(phenotype)<-intnames[names(phenotype)]
#     for(i in names(listOfRegulons)){
#       listOfRegulons[[i]]<-as.integer(intnames[listOfRegulons[[i]]])
#     }    
#   }
  #OBS:names should be provided as integer values
  pheno.names <- as.integer(names(phenotype))
  ##tag the gene sets that can be used in the analysis, i.e. those 
  ##that are smaller than the size of the gene list and that have more 
  ##than 'minRegulonSize' elements that can be found in the phenotype	
  nRegulons <- length(listOfRegulons)
  tagRegulons <- rep(FALSE, nRegulons)
  tagRegulons[which(unlist(lapply(listOfRegulons, length)) < 
    length(phenotype))] <- TRUE
  tagRegulons[which(unlist(lapply(lapply(listOfRegulons, 
    intersect, y=pheno.names), length)) < minRegulonSize)] <- FALSE
  ##check that there are actually some gene sets that pass the max 
  ##and min cutoffs
  n.tagRegulons <- sum(tagRegulons)
  if(n.tagRegulons == 0) 
    warning(paste("There are no gene sets in your collection",
                  " that pass the cutoffs on size", sep=""))
  if(n.tagRegulons > 0) {
    ##Generate a matrix to store the permutation-based scores, with 
    ##one row for each gene set (that has been tagged) and one column 
    ##for each permutation	
    scoresperm <- matrix(rep(0, (nPermutations * n.tagRegulons)), nrow=n.tagRegulons)
    rownames(scoresperm) <- names(listOfRegulons)[which(tagRegulons)]
    ##Generate a vector to store the experimental scores
    ##one entry for each gene set (that has been tagged)
    scoresObserved <- rep(0, n.tagRegulons)
    names(scoresObserved) <- names(listOfRegulons)[which(tagRegulons)]
    ##Compute the scores	
    ##create permutation gene list
    perm.gL <- sapply(1:nPermutations, function(n) pheno.names[
      sample(1:length(phenotype), length(phenotype),replace=FALSE)])
    perm.gL<-cbind(pheno.names,perm.gL)
    ##check if package snow has been loaded and a cluster object 
    ##has been created for HTSanalyzeR
    if(isPar && n.tagRegulons>1) {
      if(verbose)cat("-Performing GSEA (parallel version  - ProgressBar not available) ...\n")
      if(verbose)cat("--For", length(listOfRegulons[which(tagRegulons)]), "regulons... \n")      
      scores <- gseaScoresBatchParallel4RTN(geneList=phenotype, geneNames.perm = perm.gL,
                                        collectionOfGeneSets=listOfRegulons[which(tagRegulons)],
                                        exponent=exponent,nPermutations=nPermutations)
      sapply(1:n.tagRegulons, function(i){
        scoresperm[i,]<<-unlist(scores["scoresperm",i])
        scoresObserved[i]<<-unlist(scores["scoresObserved",i])
        NULL
      })
    } else {
      if(verbose) cat("-Performing gene set enrichment analysis ...\n")
      if(verbose) cat("--For", length(listOfRegulons[which(tagRegulons)]), "regulons \n")
      if(verbose) pb <- txtProgressBar(style=3)
      for(i in 1:n.tagRegulons) {
        scores <- gseaScoresBatch4RTN(geneList=phenotype, geneNames.perm=perm.gL, 
                                  geneSet=listOfRegulons[[which(tagRegulons)[i]]],
                                  exponent=exponent, nPermutations=nPermutations)
        scoresObserved[i] <- scores$scoresObserved
        scoresperm[i,] <- scores$scoresperm
        if(verbose) setTxtProgressBar(pb, i/n.tagRegulons)
      }	
      if(verbose) close(pb)
    }
  } else {
    scoresObserved <- NULL
    scoresperm <- NULL
  }
  return(list("Observed.scores" = scoresObserved , "Permutation.scores" = scoresperm))	
}

##This function computes observed and permutation-based scores 
gsea2tna <- function(listOfRegulons, phenotype, exponent=1, nPermutations=1000, verbose1=TRUE, verbose2=TRUE) {   
  ##convert names to integers if parallel  
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
  isPar <- b1 && b2

  #OBS:names should be provided as integer values
  pheno.names <- as.integer(names(phenotype))
  nRegulons <- length(listOfRegulons)
  if(nRegulons == 0) 
    warning(paste("There are no gene sets in your collection",
                  " that pass the cutoffs on size", sep=""))
  if(nRegulons > 0){
    ##Generate a matrix to store the permutation-based scores, with 
    ##one row for each gene set (that has been tagged) and one column 
    ##for each permutation	
    scoresperm <- matrix(rep(0, (nPermutations * nRegulons)), nrow=nRegulons)
    rownames(scoresperm) <- names(listOfRegulons)
    ##Generate a vector to store the experimental scores
    ##one entry for each gene set (that has been tagged)
    scoresObserved <- rep(0, nRegulons)
    names(scoresObserved) <- names(listOfRegulons)
    ##Compute the scores	
    ##create permutation gene list
    perm.gL <- sapply(1:nPermutations, function(n) pheno.names[
      sample(1:length(phenotype), length(phenotype),replace=FALSE)])
    perm.gL<-cbind(pheno.names,perm.gL)
    ##check if package snow has been loaded and a cluster object 
    ##has been created for HTSanalyzeR
    if(isPar && nRegulons>1) {
      if(verbose1 && verbose2)cat("-Performing two-tailed GSEA (parallel version  - ProgressBar not available) ...\n")
      if(verbose1 && verbose2)cat("--For", length(listOfRegulons), "regulons... \n")      
      scores <- gseaScoresBatchParallel4RTN(geneList=phenotype, geneNames.perm = perm.gL,
                                            collectionOfGeneSets=listOfRegulons,
                                            exponent=exponent,nPermutations=nPermutations)
      sapply(1:nRegulons, function(i){
        scoresperm[i,]<<-unlist(scores["scoresperm",i])
        scoresObserved[i]<<-unlist(scores["scoresObserved",i])
        NULL
      })
    } else {
      if(verbose1 && verbose2) cat("-Performing two-tailed GSEA analysis ...\n")
      if(verbose1 && verbose2) cat("--For", length(listOfRegulons), "regulons \n")
      if(verbose1) pb <- txtProgressBar(style=3)
      for(i in 1:nRegulons) {
        scores <- gseaScoresBatch4RTN(geneList=phenotype, geneNames.perm=perm.gL, 
                                      geneSet=listOfRegulons[[i]],
                                      exponent=exponent, nPermutations=nPermutations)
        scoresObserved[i] <- scores$scoresObserved
        scoresperm[i,] <- scores$scoresperm
        if(verbose1) setTxtProgressBar(pb, i/nRegulons)
      }	
      if(verbose1) close(pb)
    }
  } else {
    scoresObserved <- NULL
    scoresperm <- NULL
  }
  return(list("Observed.scores" = scoresObserved , "Permutation.scores" = scoresperm))	
}

##This function computes observed and permutation-based scores 
# gsea2cmap <- function(listOfRegulons, phenotype, exponent=1, nPermutations=1000, verbose1=TRUE, verbose2=TRUE) {   
#   ##convert names to integers if parallel  
#   b1<-"package:snow" %in% search()
#   b2<-tryCatch({
#     cl<-getOption("cluster")
#     cl.check<-FALSE
#     if(is(cl, "cluster")){
#       cl.check <- all( sapply(1:length(cl),function(i)isOpen(cl[[i]]$con) ) == TRUE )
#     }
#     cl.check
#   }, error=function(e){ FALSE 
#   })
#   isPar <- b1 && b2
#   
#   #OBS:names should be provided as integer values
#   pheno.names <- as.integer(names(phenotype))
#   nRegulons <- length(listOfRegulons)
#   if(nRegulons == 0) 
#     warning(paste("There are no gene sets in your collection",
#                   " that pass the cutoffs on size", sep=""))
#   if(nRegulons > 0){
#     ##Generate a matrix to store the permutation-based scores, with 
#     ##one row for each gene set (that has been tagged) and one column 
#     ##for each permutation  
#     scoresperm <- matrix(rep(0, (nPermutations * nRegulons)), nrow=nRegulons)
#     rownames(scoresperm) <- names(listOfRegulons)
#     ##Generate a vector to store the experimental scores
#     ##one entry for each gene set (that has been tagged)
#     scoresObserved <- rep(0, nRegulons)
#     names(scoresObserved) <- names(listOfRegulons)
#     ##Compute the scores	
#     ##create permutation gene list
#     perm.gL <- sapply(1:nPermutations, function(n) pheno.names[
#       sample(1:length(phenotype), length(phenotype),replace=FALSE)])
#     perm.gL<-cbind(pheno.names,perm.gL)
#     ##check if package snow has been loaded and a cluster object 
#     ##has been created for HTSanalyzeR
#     if(isPar && nRegulons>1) {
#       if(verbose1 && verbose2)cat("-Performing two-tailed GSEA (parallel version  - ProgressBar not available) ...\n")
#       if(verbose1 && verbose2)cat("--For", length(listOfRegulons), "regulons... \n")      
#       scores <- gseaScoresBatchParallel4RTN(geneList=phenotype, geneNames.perm = perm.gL,
#                                             collectionOfGeneSets=listOfRegulons,
#                                             exponent=exponent,nPermutations=nPermutations)
#       sapply(1:nRegulons, function(i){
#         scoresperm[i,]<<-unlist(scores["scoresperm",i])
#         scoresObserved[i]<<-unlist(scores["scoresObserved",i])
#         NULL
#       })
#     } else {
#       if(verbose1 && verbose2) cat("-Performing two-tailed GSEA analysis ...\n")
#       if(verbose1 && verbose2) cat("--For", length(listOfRegulons), "regulons \n")
#       if(verbose1) pb <- txtProgressBar(style=3)
#       for(i in 1:nRegulons) {
#         scores <- gseaScoresBatch4RTN(geneList=phenotype, geneNames.perm=perm.gL, 
#                                       geneSet=listOfRegulons[[i]],
#                                       exponent=exponent, nPermutations=nPermutations)
#         scoresObserved[i] <- scores$scoresObserved
#         scoresperm[i,] <- scores$scoresperm
#         if(verbose1) setTxtProgressBar(pb, i/nRegulons)
#       }	
#       if(verbose1) close(pb)
#     }
#   } else {
#     scoresObserved <- NULL
#     scoresperm <- NULL
#   }
#   return(list("Observed.scores" = scoresObserved , "Permutation.scores" = scoresperm))	
# }

##------------------------------------------------------------------------
##This function computes pairwise ES scores for synergy analysis
##..here synergy is inferred from the 'effectsize' by
##..comparing the enrichment score of the intesect related to the union!
##..i.e. the null distribution is computed by resampling the union!
pairwiseSynergy <- function(collectionsOfPairs, phenotype, exponent=1, 
                      nPermutations=1000, minIntersectSize=1, verbose=TRUE) {
  
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
  ##min intersect
  minter<-as.integer(minIntersectSize)*0.01
  if(b1 && b2 && length(collectionsOfPairs)>1) {
    if(verbose)cat("-Performing synergy analysis (parallel version - ProgressBar not available) ...\n")
    if(verbose)cat("--For", length(collectionsOfPairs), "regulon pairs ...\n")
    gseaScores<-gseaScores
    #get permutation scores
    permScores<-parSapply(getOption("cluster"),1:length(collectionsOfPairs), function(i){
      gsPair<-collectionsOfPairs[[i]]
      bl1 <- length(gsPair[[2]]) >= 1
      bl2 <- length(gsPair[[2]]) >= minter*length(gsPair[[1]])
      if( bl1 && bl2 ){
        sapply(1:nPermutations, function(j){
          gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                     geneSet=sample(gsPair[[1]],length(gsPair[[2]])) )
        })
      } else {
        rep(NA,nPermutations)
      }
    })
    permScores<-t(permScores)
  } else {
    if(verbose)cat("-Performing synergy analysis ...\n")
    if(verbose)cat("--For", length(collectionsOfPairs), "regulon pairs ...\n")
    if(verbose) pb <- txtProgressBar(style=3)
    permScores<-matrix(rep(0,length(collectionsOfPairs) * nPermutations), ncol=nPermutations)
    #get permutation scores
    sapply(1:length(collectionsOfPairs), function(i){
      gsPair<-collectionsOfPairs[[i]]
      bl1 <- length(gsPair[[2]]) >= 1
      bl2 <- length(gsPair[[2]]) >= minter*length(gsPair[[1]])
      if(bl1 && bl2){
        permScores[i,]<<-sapply(1:nPermutations, function(j){
          gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                     geneSet=sample(gsPair[[1]],length(gsPair[[2]])) )
        })
      } else {
        permScores[i,]<<-rep(NA,nPermutations)
      }
      if(verbose) setTxtProgressBar(pb, i/length(collectionsOfPairs))
      NULL
    })
    if(verbose) close(pb)
  }
  rownames(permScores)<-names(collectionsOfPairs)
  #get observed scores
  observedScores<-sapply(1:length(collectionsOfPairs), function(i){
    gsPair<-collectionsOfPairs[[i]]
    bl1 <- length(gsPair[[2]])>= 1
    bl2 <- length(gsPair[[2]])>= minter*length(length(gsPair[[1]]))
    if(bl1 && bl2){
      gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                 geneSet=gsPair[[2]])
    } else {
      return(NA)
    }
  })
  names(observedScores)<-names(collectionsOfPairs)
  #return scores
  return(list("Subset.scores"=observedScores, "Allset.scores"=permScores))
}

##------------------------------------------------------------------------
##This function computes pairwise ES scores for shadow analysis
##..it compares the enrichment score of the unique genes of a regulon with
##..all genes of the same regulon! i.e. the null distribution is computed 
##..by taking subsamples from the regulon!
##..(obs.: subsamples of the same size of the unique set!)
pairwiseShadow <- function(collectionsOfPairs, phenotype, exponent=1, 
                           nPermutations=1000, minIntersectSize=1, verbose=TRUE) {
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
  ##min intersect
  minter<-as.integer(minIntersectSize)*0.01
  if( b1 && b2 && length(collectionsOfPairs)>1) {
    gseaScores<-gseaScores
    #get permutation scores
    permScores<-parSapply(getOption("cluster"),1:length(collectionsOfPairs), function(i){
      gsPair<-collectionsOfPairs[[i]]
      bl1 <- length(gsPair[[1]])-length(gsPair[[2]]) >= 1
      bl2 <- length(gsPair[[1]])-length(gsPair[[2]]) >= minter*length(gsPair[[1]])
      if(bl1 && bl2){
        sapply(1:nPermutations, function(j){
          gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                     geneSet=sample(gsPair[[1]],length(gsPair[[2]])) )
        })
      } else {
        rep(NA,nPermutations)
      }
    })
    permScores<-t(permScores)
  } else {
    #get permutation scores
    if(verbose) pb <- txtProgressBar(style=3)
    permScores<-matrix(rep(NA,length(collectionsOfPairs) * nPermutations), ncol=nPermutations)
    sapply(1:length(collectionsOfPairs), function(i){
      gsPair<-collectionsOfPairs[[i]]
      bl1 <- length(gsPair[[1]])-length(gsPair[[2]]) >= 1
      bl2 <- length(gsPair[[1]])-length(gsPair[[2]]) >= minter*length(gsPair[[1]])
      if(bl1 && bl2){
        permScores[i,]<<-sapply(1:nPermutations, function(j){
          gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                     geneSet=sample(gsPair[[1]],length(gsPair[[2]])) )
        })
      }
      if(verbose) setTxtProgressBar(pb, i/length(collectionsOfPairs))
      NULL
    })
    if(verbose) close(pb)
  }
  rownames(permScores)<-names(collectionsOfPairs)
  #get observed scores
  observedScores<-sapply(1:length(collectionsOfPairs), function(i){
    gsPair<-collectionsOfPairs[[i]]
    bl1 <- length(gsPair[[1]])-length(gsPair[[2]]) >= 1
    bl2 <- length(gsPair[[1]])-length(gsPair[[2]]) >= minter*length(gsPair[[1]])   
    if(bl1 && bl2){
      gseaScores4RTN(geneList=phenotype,exponent=exponent, mode="score", 
                 geneSet=gsPair[[2]])
    } else {
      return(NA)
    }
  })
  names(observedScores)<-names(collectionsOfPairs)
  #return scores
  return(list("Subset.scores"=observedScores, "Allset.scores"=permScores))
}
##------------------------------------------------------------------------
##This function gets rid of the duplicates in a gene list.
tna.duplicate.remover <- function(phenotype, method = "max") {
  ##Get the unique names and create a vector that will store the 
  ##processed values corresponding to those names
  pheno.names <- names(phenotype)
  datanames <- unique(pheno.names)
  data.processed <- rep(0, length(datanames))
  names(data.processed) <- datanames
  data.processed2 <- data.processed
  l.data.processed <- length(data.processed)	
  ##If the absolute value of the min is bigger than the absolute 
  ##value of the max, then it is the min that is kept
  if(method=="max") {
    abmax<-function(x){
      imax<-max(x);imin<-min(x)
      ifelse(imax>abs(imin),imax,imin)
    }
    tp<-aggregate(phenotype,by=list(names(phenotype)),abmax)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1] 
  } else if(method=="min") {
    abmin<-function(x){
      imax<-max(x);imin<-min(x)
      ifelse(imax>abs(imin),imin,imax)
    }
    tp<-aggregate(phenotype,by=list(names(phenotype)),abmin)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1] 
  } else if(method=="average") {
    tp<-aggregate(phenotype,by=list(names(phenotype)),mean)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1]
  }	
  data.processed
}

##------------------------------------------------------------------------
##This function compute the nominal p-value associated with a GSEA for a
##collection of regulons, from the outputs of collectionGsea
tna.permutation.pvalues <- function(permScores, dataScores){
	##check 'permScores'
	if(!is.matrix(permScores)) 
		stop("NOTE: the argument permScores should be a matrix!")
	#check 'dataScores'
	if(!is.vector(dataScores) || is.null(names(dataScores))) 
		stop("NOTE: the argument dataScores should be a named vector!")
	if(!is.integer(dataScores) && !is.numeric(dataScores)) 
		stop("NOTE: the argument dataScores should be a numerical vector!")
	##check that the dimensions of permScores and dataScores match to 
	##what is expected
	l.dataScores <- length(dataScores)
	nPerm <- ncol(permScores)
	if(nrow(permScores) != l.dataScores) 
		warning("The number of rows of the permScores matrix is not ",
			"equal to the length of dataScores")
	##initialize a pvalues vector	
	pval <- rep(1, l.dataScores)
	##go through each element of dataScores and see how many permScores 
	##in the corresponding row are higher (in magnitude)
	##this is done separately for negative and positive scores
	##the scores of zero are simply not treated because they are left 
	##at the initial pvalue of 1.
	##obs. pseudocounts are added to avoid P-values of zero.
	valid.id <- which(!is.na(dataScores) & dataScores!=0)
	sapply(valid.id, function(i) {
	  pval[i]<<-ifelse(
			dataScores[i] > 0,
			(sum(permScores[i, ] > dataScores[i])+1)/(nPerm+1),
			(sum(permScores[i, ] < dataScores[i])+1)/(nPerm+1)
		)
    NULL
	})
	names(pval)<-names(dataScores)		
	return(pval)
}
tni.permutation.pvalues <- function(permScores, dataScores){
  if(!is.matrix(permScores)) 
    stop("NOTE: the argument permScores should be a matrix!")
  if(!is.matrix(dataScores)) 
    stop("NOTE: the argument dataScores should be a matrix!")
  if(!is.integer(dataScores) && !is.numeric(dataScores)) 
    stop("NOTE: the argument dataScores should be a numerical vector!")
  if(nrow(permScores) != ncol(dataScores)) 
    warning("The number of rows of the permScores matrix is not ",
            "equal to the ncols of dataScores")
  pvalMat<-sapply(1:nrow(dataScores), function(i){
    dtScores<-dataScores[i,]
    l.dataScores <- length(dtScores)
    nPerm <- ncol(permScores)
    ##initialize a pvalues vector  
    pval <- rep(1, l.dataScores)
    valid.id <- which(!is.na(dtScores) & dtScores!=0)
    sapply(valid.id, function(j) {
      pval[j]<<-ifelse(
        dtScores[j] > 0,
        (sum(permScores[j, ] > dtScores[j])+1)/(nPerm+1),
        (sum(permScores[j, ] < dtScores[j])+1)/(nPerm+1)
      )
      NULL
    })
    pval
  })
  colnames(pvalMat)<-colnames(dataScores)
  rownames(pvalMat)<-rownames(dataScores)
  return(pvalMat)
}

##------------------------------------------------------------------------
##This function compute the nominal p-value associated with a GSEA for a
##collection of regulons, from the outputs of collectionGsea
##minor changes for shadow analysis!
tna.permutation.pvalues.shadow <- function(permScores, dataScores){
  ##check 'permScores'
  if(!is.matrix(permScores)) 
    stop("NOTE: the argument permScores should be a matrix!")
  #check 'dataScores'
  if(!is.vector(dataScores) || is.null(names(dataScores))) 
    stop("NOTE: the argument dataScores should be a named vector!")
  if(!is.integer(dataScores) && !is.numeric(dataScores)) 
    stop("NOTE: the argument dataScores should be a numerical vector!")
  ##check that the dimensions of permScores and dataScores match to 
  ##what is expected
  l.dataScores <- length(dataScores)
  nPerm <- ncol(permScores)
  if(nrow(permScores) != l.dataScores) 
    warning("The number of rows of the permScores matrix is not ",
            "equal to the length of dataScores")
  ##initialize a pvalues vector	
  pval <- rep(1, l.dataScores)
  ##go through each element of dataScores and see how many permScores 
  ##in the corresponding row are higher (in magnitude)
  ##this is done separately for negative and positive scores
  ##the scores of zero are simply not treated because they are left 
  ##at the initial pvalue of 1.
  ##obs. pseudocounts are added to avoid P-values of zero.
  valid.id <- which(!is.na(dataScores))
  sapply(valid.id, function(i) {
    pval[i]<<-ifelse(
      dataScores[i] > 0,
      (sum(permScores[i, ] <= dataScores[i])+1)/(nPerm+1),
      (sum(permScores[i, ] >= dataScores[i])+1)/(nPerm+1)
    )
    NULL
  })
  names(pval)<-names(dataScores)		
  return(pval)
}

##------------------------------------------------------------------------
##This function applies the hyperGeoTest function to an entire list of 
##regulons and returns a data frame with the overlap among regulon pairs
tna.hyper.pairs <- function(listOfRegulons, universe, 
  minRegulonSize = 15, pAdjustMethod = "BH", verbose = TRUE) {
  l.Regulons <- length(listOfRegulons)
  regulon.size <- unlist(
    lapply(
      lapply(listOfRegulons, intersect, y = universe), length
    )
  )
  if(all(regulon.size < minRegulonSize))
    stop(paste("NOTE: the largest number of overlapped genes ",
               "with universe is: ", max(regulon.size), ", which is < ", 
               minRegulonSize, "!\n", sep = ""))
  regulon.filtered <- which(regulon.size >= minRegulonSize)
  ##if verbose, create a progress bar to monitor computation progress
  if(verbose) pb <- txtProgressBar(style=3)
  results <- matrix(, nrow=0, ncol=8)
  colnames(results) <- c("Regulon1","Regulon2","Universe.Size", "R1.Size", "R2.Size", 
                         "Expected.Overlap", "Observed.Overlap", "Pvalue")
  for(i in 1:length(regulon.filtered)){
      if(verbose) setTxtProgressBar(pb, i/length(regulon.filtered))
      r1<-regulon.filtered[i]
      r.filter<-which(regulon.filtered > r1)
      res <- t(
            sapply(r.filter, 
                   function(r2) {
                     tna.hyper(listOfRegulons[r1], universe, listOfRegulons[r2])
                   }
            )
      )
      if(ncol(res)==6){
        res<-data.frame(Regulon1=names(r1),Regulon2=rownames(res),res,stringsAsFactors=FALSE)
        results<-rbind(results,res)       
      }
  }
  if(verbose) 
    close(pb)
  if(nrow(results) > 0) {
    ##Adjust pvalues
    adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
    results <- cbind(results, adjPvals)
    colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
    results <- results[order(results[, "Adjusted.Pvalue"]), , drop=FALSE]	
    rownames(results)<-1:nrow(results)
  } else {
    results <- matrix(, nrow=0, ncol=7)
    colnames(results) <- c("Universe.Size", 
                           "R1.Size", "R2.Size", "Expected.Overlap", 
                           "Observed.Overlap", "Pvalue", "Adjusted.Pvalue")
    
  }
  return(results)
}

##------------------------------------------------------------------------
##This function takes two gene sets (i.e. regulons), a vector 
##containing the size of the gene universe, and compute the number of 
##genes expected to occur in both regulons, the actual observed overlap, 
##and the pvalue from a hypergeometric test.
tna.hyper <- function(regulon1, universe, regulon2) {
  ##number of genes in universe
  N <- length(universe)			
  ##remove genes from gene set that are not in universe			
  regulon1 <- intersect(regulon1[[1]], universe) 
  regulon2 <- intersect(regulon2[[1]], universe)
  ##size of gene set	
  m <- length(regulon1) 							
  Nm <- N-m	
  ##regulon2 in gene set
  overlap <- intersect(regulon1, regulon2) 	
  ##number of hits between regulons		
  k <- length(overlap) 							
  n <- length(regulon2)	
  HGTresults <- phyper(k, m, Nm, n, lower.tail = F)
  ex <- (n/N)*m
  if(m == 0) HGTresults <- NA
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe.Size", "R1.Size", "R2.Size", 
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(hyp.vec)
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
  # remove proves without annotation
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
#--------------------------------------------------------------------
#The next 2 functions were extracted from the old version of HTSanalyzeR
#due to a compatibility issue:
#...the authors introduced a correction in "phyper" function that is not
#...required for RTN, and as result the corrected p-values are no longer
#...reproducible! here the functions have minor changes only to keep the 
#...compatibility!
#--------------------------------------------------------------------
##This function performs hypergeometric tests for over-representation 
##of hits, on a list of gene sets. This function applies the 
##hyperGeoTest function to an entire list of gene sets and returns a 
##data frame.
multiHyperGeoTest4RTN <- function(collectionOfGeneSets, universe, hits, 
                              minGeneSetSize = 15, pAdjustMethod = "BH", verbose = TRUE) {
  geneset.size <- unlist(
    lapply(
      lapply(collectionOfGeneSets, intersect, y = universe), length
    )
  )
  if(all(geneset.size < minGeneSetSize)){
    stop(paste("NOTE: the largest number of overlapped genes of gene ",
               "sets with universe is: ", max(geneset.size), ", which is < ", 
               minGeneSetSize, "!\n", sep = ""))
  }
  geneset.filtered <- which(geneset.size >= minGeneSetSize)
  l.GeneSet <- length(geneset.filtered)
  ##if verbose, create a progress bar to monitor computation progress
  if(verbose) pb <- txtProgressBar(style=3)
  results <- t(
    sapply(geneset.filtered, 
           function(i) {
             if(verbose) 
               setTxtProgressBar(pb, i/l.GeneSet)		
             hyperGeoTest4RTN(collectionOfGeneSets[i], universe, hits)
           }
    )
  )
  if(verbose) 
    close(pb)
  if(length(results) > 0) {
    ##results <- t(as.data.frame(results))
    ##rownames(results) <- names(collectionOfGeneSets)
    ##remove gene sets with genes < minGeneSetSize in the geneList
    ##Adjust pvalues
    adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
    results <- cbind(results, adjPvals)
    colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
    results <- results[order(results[, "Adjusted.Pvalue"]), , drop=FALSE]		
  } else {
    reuslts <- matrix(, nrow=0, ncol=7)
    colnames(results) <- c("Universe Size", "Gene Set Size", 
                           "Total Hits", "Expected Hits", "Observed Hits", "Pvalue", 
                           "Adjusted.Pvalue")
  }
  return(results)
}
##This function takes in a single gene set (GeneSet), a vector 
##(GeneList) of gene symbols for all tested genes, a vector of "hits" 
##(hits), and a p-value adjustment method. It outputs a vector 
##containing the size of the gene universe, the size of the gene set 
##within this universe (i.e. how many genes from the universe map to 
##this gene set), the total number of hits, the number of hits expected 
##to occur in the gene set, the actual hits observed in the gene set, 
##and the pvalue from a hypergeometric test.
hyperGeoTest4RTN <- function(geneSet, universe, hits) {
  ##number of genes in universe
  N <- length(universe) 			
  ##remove genes from gene set that are not in universe			
  geneSet <- intersect(geneSet[[1]], universe) 
  ##size of gene set	
  m <- length(geneSet) 							
  Nm <- N-m	
  ##hits in gene set
  overlap <- intersect(geneSet, hits) 	
  ##number of hits in gene set		
  k <- length(overlap) 							
  n <- length(hits)	
  HGTresults <- phyper(k, m, Nm, n, lower.tail = F)
  ex <- (n/N)*m
  if(m == 0) HGTresults <- NA
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe Size", "Gene Set Size", "Total Hits", 
                      "Expected Hits", "Observed Hits", "Pvalue")
  return(hyp.vec)
}
#--------------------------------------------------------------------
#The next 3 functions have being extracted from HTSanalyzeR
#due to a compatibility issue between igraph/igraph0 versions!
#...in HTSanalyzeR, the dependency of igraph was changed to igraph0 
#...to adapt to the change in BioNet. However, RTN requires the new
#...igraph version and, as result, no other HTSanalyzeR function is
#...required.
#--------------------------------------------------------------------
##This function computes enrichment scores for GSEA, running score and 
##position of hits for a gene set.
gseaScores4RTN <- function(geneList, geneSet, exponent=1, mode="score", geneSetMode=NULL) {
  ##The geneSet should be a subset of the gene universe, i.e. we keep 
  ##only those element of the gene set that appear in the geneList		
  geneSet<-intersect(names(geneList), geneSet)
  ##Compute the size of the gene set and of the genelist	
  nh <- length(geneSet)
  N <- length(geneList)
  ##Initialize the ES, runningES and the Phit and Pmiss by position 
  ##(the actual values of Phit and Pmiss are actually cumulative sums 
  ##of these 'by position' values)	
  ES <- 0
  Phit <- rep(0, N)
  Pmiss <- rep(0, N)
  runningES <- rep(0, N)
  hits <- rep(FALSE, N)
  hitsmode <- rep(0, N)
  ##Stop if the geneSet is larger than the gene universe	
  if(nh > N) {
    stop("NOTE: gene set is larger than gene list!")
  } else {
    ##Compute the positions of the hits in the geneList (0 if there 
    ##is no match, 1 if there is a match)	
    hits[which(!is.na(match(names(geneList), geneSet)))] <- TRUE
    if(!is.null(geneSetMode)){
      hitsmode[which(!is.na(match(names(geneList), names(geneSetMode[geneSetMode>0]))))] <-  1
      hitsmode[which(!is.na(match(names(geneList), names(geneSetMode[geneSetMode<0]))))] <- -1
    }
    ##If sum(hits)=0 then there is no match between geneList and 
    ##geneSet, and all scores stay at 0.		
    if(sum(hits)!=0) {
      ##Fill the Phit by position		
      Phit[which(hits)]<-abs(geneList[which(hits)])^exponent
      NR=sum(Phit)
      ##Fill the Pmiss by positions			
      Pmiss[which(!hits)]<-1/(N-nh)
      ##Do the cumulative sums	and compute the runningES		
      Phit=cumsum(Phit/NR)
      Pmiss=cumsum(Pmiss)
      runningES<-Phit-Pmiss
      ##Compute the maximal (positive) and minimal (or maximal 
      ##negative) values of the ES, and choose which one is kept			
      ESmax<-max(runningES)
      ESmin<-min(runningES)
      ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
    }
  }
  ##Return the relevant information according to mode  	
  if(mode=="score"){
    return(ES)
  } else if(mode=="graph"){
    if(!is.null(geneSetMode))hits=hitsmode
    return(list("enrichmentScore"=ES, "runningScore"=runningES, "positions"=as.integer(hits)))
  }
}
##------------------------------------------------------------------------------
gseaScores4CMAP <- function(geneList, geneSet, exponent) {
  nh <- length(geneSet)
  N <- length(geneList)  
  ES <- 0
  Phit <- rep(0, N)
  Pmiss <- rep(0, N)
  runningES <- rep(0, N)
  hits <- rep(FALSE, N)
  hits[which(names(geneList)%in%geneSet)] <- TRUE  
  if(sum(hits)!=0) {
    Phit[which(hits)]<-abs(geneList[which(hits)])^exponent
    NR=sum(Phit)	
    Pmiss[which(!hits)]<-1/(N-nh)	
    Phit=cumsum(Phit/NR)
    Pmiss=cumsum(Pmiss)
    runningES<-Phit-Pmiss		
    ESmax<-max(runningES)
    ESmin<-min(runningES)
    ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
  }
  return(ES)
}
#--------------------------------------------------------------------
##This function computes enrichment score for both input 'geneList' 
##and its permutations for one gene set.
gseaScoresBatch4RTN <- function(geneList, geneNames.perm, geneSet, 
                            exponent=1, nPermutations=1000) {
  if(!is.matrix(geneNames.perm))
    stop("NOTE: 'geneNames.perm' should be a matrix!\n")
  if(ncol(geneNames.perm) != (nPermutations+1))
    stop("NOTE: the No of columns of 'geneNames.perm' should be equal to 'nPermutations'!\n")
  geneList.names <- as.integer(names(geneList)) #ja foi transformado em inteiros!
  ##Compute the size of the gene set and of the genelist	
  nh<-length(geneSet)
  N<-length(geneList)
  ##The geneSet should be a subset of the gene universe, i.e. we keep 
  ##only those element of the gene set that appear in the geneList		
  geneSet<-intersect(geneList.names,geneSet)
  ES<-rep(0,nPermutations+1)
  Phit<-matrix(0,nrow=N,ncol=nPermutations+1)
  Pmiss<-Phit
  runningES<-NULL
  
  if(nh>N) {
    stop("NOTE: gene set is larger than gene list!")
  } else {
    hits <- matrix(FALSE, nrow = N, ncol = nPermutations+1)
    hits[which(!is.na(match(geneNames.perm, geneSet)))] <- TRUE	
    hits <- matrix(hits,ncol = nPermutations+1 , byrow = FALSE)		
    if(sum(hits[,1]) > 0) {
      junk <- sapply(1:(nPermutations+1), function(i) 
        Phit[which(hits[, i]), i] <<- 
          abs(geneList[which(hits[, i])])^exponent)	
      NR <- colSums(Phit)		
      Pmiss[which(!hits)] <- 1/(N-nh)		
      Pmiss <- sapply(1:(nPermutations+1), function(i) 
        cumsum(Pmiss[, i]))
      Phit <- sapply(1:(nPermutations+1), function(i) 
        cumsum(Phit[, i])/NR[i])		
      runningES <- Phit - Pmiss		
      ESrange <- sapply(1:(nPermutations+1), function(i) 
        range(runningES[, i]))
      ES <- sapply(1:(nPermutations+1), function(i) 
        ESrange[which.max(abs(ESrange[, i])), i])	
      if(is.list(ES)) ES<-unlist(ES)
    }
  }
  #Return the relevant information according to mode		
  ES<-list(scoresObserved=ES[1], scoresperm=ES[2:(nPermutations+1)])
  return(ES)	
}
#--------------------------------------------------------------------
##This function computes enrichment scores for both input 'geneList' 
##and their permutations for multiple gene sets in parallel
gseaScoresBatchParallel4RTN <- function(geneList, geneNames.perm, 
                                    collectionOfGeneSets, 
                                    exponent = 1, nPermutations = 1000) {
  if(!is.matrix(geneNames.perm))
    stop("NOTE: 'geneNames.perm' should be a matrix!\n")
  if(ncol(geneNames.perm)!=(nPermutations+1))
    stop("NOTE: the No of columns of 'geneNames.perm' should be equal to 'nPermutations'!\n")	
  ##local function for computation of gsea scores with a single core
  gseaScoresBatchLocal <- function(geneList, geneNames.perm, geneSet, 
                                   exponent, nPermutations) {	
    geneList.names <- as.integer(names(geneList))
    ##Compute the size of the gene set and of the genelist	
    nh <- length(geneSet)
    N <- length(geneList)
    ##The geneSet should be a subset of the gene universe, i.e. we 
    ##keep only those element of the gene set that appear in the 
    ##geneList		
    geneSet <- intersect(geneList.names, geneSet)
    ES <- rep(0, nPermutations+1)
    Phit <- matrix(0, nrow = N, ncol = nPermutations+1)
    Pmiss <- Phit
    runningES <- NULL
    
    if(nh > N)
      stop("NOTE: gene set is larger than gene list!")
    
    hits <- matrix(FALSE, nrow = N, ncol = nPermutations+1) 	
    hits[which(!is.na(match(geneNames.perm, geneSet)))] <- TRUE	
    hits <- matrix(hits, ncol = nPermutations+1, byrow = FALSE)		
    if(sum(hits[,1]) > 0) {
      junk <- sapply(1:(nPermutations+1), function(i) 
        Phit[which(hits[, i]), i] <<- 
          abs(geneList[which(hits[, i])])^exponent)	
      NR <- colSums(Phit)		
      Pmiss[which(!hits)] <- 1/(N-nh)		
      Pmiss <- sapply(1:(nPermutations+1), function(i) 
        cumsum(Pmiss[, i]))
      Phit <- sapply(1:(nPermutations+1), function(i) 
        cumsum(Phit[, i])/NR[i])		
      runningES <- Phit-Pmiss		
      ESrange <- sapply(1:(nPermutations+1), function(i) 
        range(runningES[, i]))
      ES <- sapply(1:(nPermutations+1), function(i) 
        ESrange[which.max(abs(ESrange[,i])),i])	
      if(is.list(ES)) ES <- unlist(ES)
    }	
    ##Return the relevant information according to mode		
    ES <- list(scoresObserved = ES[1], scoresperm = ES[2:(nPermutations+1)])
    return(ES)	
  }
  #parallel computing
  scores <- parSapply(getOption("cluster"), 1:length(collectionOfGeneSets), 
                      function(i) {
                        gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm, 
                                             geneSet = collectionOfGeneSets[[i]], 
                                             exponent = exponent, 
                                             nPermutations = nPermutations)
                      }
  )
  return(scores)
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
  tnet<-tnet[c(tfs,onlytar),]
  tmat<-t(tnet[onlytar,])
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
tni.mmap<-function(mnet,cnet,othertfs){
  cnet[cnet<0]=-1;cnet[cnet>0]=1
  mnet<-mnet[,colnames(cnet)]
  mnet<-mnet[rownames(cnet),]
  #---merge/simplify cnet and mnet
  mcnet<-mnet*10
  mcnet[mcnet==0]<-cnet[mcnet==0]
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
  }
  #---set ordering
  mcnet<-mcnet[c(tfs,onlytar),]
  #---
  elist<-NULL
  emode<-NULL
  etype<-NULL  
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
          emode<<-c(emode,0)
          etype<<-c(etype,"TF.TF")
          elist<<-rbind(elist,c(tfs[j],tfs[i]))
        }
      }
      NULL
    })
    NULL
  })
  #others
  if(length(onlytar)>0){
    tpnet<-mcnet[onlytar,tfs,drop=FALSE]
    sapply(1:nrow(tpnet),function(i){
      sapply(1:ncol(tpnet),function(j){
        tp<-tpnet[i,j]
        if(abs(tp)==10){
          emode<<-c(emode,tp)
          tpp<-ifelse(tp<0,"-MD.TF","+MD.TF")
          etype<<-c(etype,tpp)
          elist<<-rbind(elist,c(onlytar[i],tfs[j]))
        } else if(abs(tp)==1){
          emode<<-c(emode,tp)
          tpp<-ifelse(tp<0,"-TF.Target","+TF.Target")
          etype<<-c(etype,tpp)
          elist<<-rbind(elist,c(tfs[j],onlytar[i]))
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
  idx<-abs(emode)==10
  emode[idx]<-emode[idx]/10
  #---get graph and set mode
  if(is.null(elist)){
    g<-graph.empty(n=0, directed=TRUE)
  } else {
    g<-igraph::graph.edgelist(elist, directed=TRUE)
    E(g)$arrowType<-emode
    E(g)$interactionType<-etype
    E(g)$modeOfAction<-NA
    E(g)$modeOfAction[E(g)$interactionType=="-MD.TF"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+MD.TF"]<- 1
    E(g)$modeOfAction[E(g)$interactionType=="TF.TF"]<- 0
    E(g)$modeOfAction[E(g)$interactionType=="-TF.Target"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+TF.Target"]<- 1
  }
  #---add disconnected tfs 
  idx<-!tfs%in%V(g)$name
  if(any(idx))g<-g+vertices(tfs[idx])
  return(g)
}

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
