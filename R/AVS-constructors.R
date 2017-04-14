
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##                      -- AVS constructor --
##                      --  LDHapMapRel27  --
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
##pre-processing for an avs
avs.preprocess.LDHapMapRel27 <- function(markers, nrand=1000, mergeColinked=TRUE, 
                                         snpop="HapMapRel27", verbose=TRUE){
  
  #---check arguments
  if(missing(markers))stop("NOTE: 'markers' is missing!") 
  markers=avs.checks(name="markers",markers)
  avs.checks(name="nrand",para=nrand)
  avs.checks(name="mergeColinked",para=mergeColinked)
  snpop=avs.checks(name="snpop",para=snpop)
  avs.checks(name="verbose",para=verbose)
  
  #---check available source for random AVS
  if(is.character(snpop)){
    opts<-c("dbSNP","HapMapRel27")
    if(!(snpop %in% opts))
      stop(paste("available 'snpop' options:", paste(opts,collapse = ", ") ),call.=FALSE)
  }
  
  #---initialize
  object <- .avs.initialize(markers,nrand,snpop,"RTNdata.LDHapMapRel27")
  object <- .checkUniverse.LDHapMapRel27(object, verbose=verbose)
  
  #---build AVS
  object <- .buildAVS.LDHapMapRel27(object, mergeColinked, verbose=verbose)
  
  #---build random AVS
  object <- .buildRandomAVS.LDHapMapRel27(object, snpop=snpop, verbose=verbose)
  
  #---get IRanges, update status and return results
  object <- .avs.mapranges(object,verbose)
  return(object)
  
}

##------------------------------------------------------------------------------
## remove markers not listed in the "marker's universe" of the LD call
.checkUniverse.LDHapMapRel27 <- function(object, verbose=TRUE){
  if(verbose){
    cat("Checking markers at HapMap data (CEU population, release #27, NCBI-B36/hg18)...\n")
  }
  obj <- data("SNPsHapMapRel27", package="RTNdata.LDHapMapRel27", envir=environment())
  SNPsHapMapRel27 <- get(obj)
  idx <- object@validatedMarkers$rsid %in% SNPsHapMapRel27
  object@validatedMarkers <- object@validatedMarkers[idx,]
  object@summary$markers[,"universe.removed"]<-nrow(object@validatedMarkers)
  return(object)
}

##------------------------------------------------------------------------------
## build AVS
.buildAVS.LDHapMapRel27<-function(object, mergeColinked, verbose=TRUE){
  reldata <- "RTNdata.LDHapMapRel27"
  mysnps <- object@validatedMarkers
  #set valid chroms
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(verbose){
    cat("Building AVS from HapMap LD data (CEU population, release #27, NCBI-B36/hg18)...\n")
  }
  variantSet<-list()
  for(chr in chrs){
    snps<-mysnps[mysnps$chrom==chr,,drop=FALSE]
    if(verbose) cat("...",nrow(snps)," marker(s) from ",chr,"!\n",sep="")
    if(nrow(snps)>0){
      ld_data<-getldata(chr,package=reldata) 
      marker1<-data.table(marker=ld_data$marker1,ord=1:nrow(ld_data))
      marker2<-data.table(marker=ld_data$marker2,ord=1:nrow(ld_data))
      setkey(marker1,'marker')
      setkey(marker2,'marker')
      variantSet[[chr]]<-lapply(1:nrow(snps),function(i){
        marker<-snps$rsid[i]
        markpos<-snps$position[i]
        ld<-ld_data[c(marker1[marker,nomatch=0]$ord,marker2[marker,nomatch=0]$ord),]
        if(nrow(ld)>0){
          tp1<-ld[ld$marker2==marker,c("marker1","position1"),drop=FALSE]
          tp2<-ld[ld$marker1==marker,c("marker2","position2"),drop=FALSE]
          markpos<-c(markpos,tp1$position1,tp2$position2)
          names(markpos)<-c(marker,tp1$marker1,tp2$marker2)
          markpos<-sort(markpos)
        } else {
          names(markpos)<-marker  
        }
        markpos
      })
      names(variantSet[[chr]])<-snps$rsid
    }
  }
  #---check co-linked AVS
  if(mergeColinked){
    if(verbose)cat("Checking co-linked markers in the AVS...\n")
    variantSet<-check.colinked(variantSet, verbose=verbose)
    mkrs<-getMarkers.vset(variantSet,getlinked=FALSE)
    object@validatedMarkers<-object@validatedMarkers[mkrs,,drop=FALSE]
    object@summary$markers[,"colinked.removed"]<-nrow(object@validatedMarkers)
    if(verbose)cat("\n")
  }
  if(verbose)cat("\n")
  object@variantSet <- variantSet
  return(object)
}

##------------------------------------------------------------------------------
## build random AVS
.buildRandomAVS.LDHapMapRel27<-function(object, snpop, verbose=TRUE){
  reldata <- "RTNdata.LDHapMapRel27"
  variantSet <- object@variantSet
  nrand <- object@para$avs$nrand
  #set valid chroms
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(verbose){
    cat("Building matched random AVS from HapMap LD data (CEU population, release #27, NCBI-B36/hg18)...\n")
  }
  #count observed clusters in the variantSet
  clCount<-unlist(sapply(1:length(variantSet),function(i){
    unlist(lapply(variantSet[[i]],length))
  }))
  ncl1<-sum(clCount==1);ncl2<-sum(clCount>1)
  #get all markers in the variantSet
  allvset<-getMarkers.vset(variantSet,TRUE)
  #get popsnp
  if(is.data.frame(snpop)){
    #custom option: 'snp population' is defined by the user!
    popsnp<-snpop[snpop$chrom%in%chrs,]
    if(nrow(popsnp)<=(length(allvset)*2))
      stop("NOTE: not enough valid markers in the input 'snpop'!",call.=FALSE)
  } else {
    #default option: pre-computed 'snp population' based on a reference dataset (eg dbSNP, HapMap)
    if(snpop=='dbSNP'){
      obj=data('popsnp1', package=reldata, envir=environment())
      popsnp <- get(obj)
    } else if(snpop=='HapMapRel27'){
      obj=data('popsnp2', package=reldata, envir=environment())
      popsnp <- get(obj)
    }
    if(nrow(popsnp)<=(length(allvset)*2))
      stop("NOTE: not enough markers to build the random sets for the input data/options!",call.=FALSE)
  }

  #remove variantSet from popsnp
  popsnp<-popsnp[!popsnp$rsid%in%allvset,]

  #get risk markers only
  nvset<-length(getMarkers.vset(variantSet,FALSE))

  #---build randomSet -- by chroms, to speed the process
  if(verbose) cat("Computing random sets...\n")
  randomSet<-list()
  for(chr in chrs){
    if(verbose) cat("",chr,"\n",sep="")
    ld_data<-getldata(chr,package=reldata)
    marker1<-data.table(marker=ld_data$marker1,ord=1:nrow(ld_data))
    marker2<-data.table(marker=ld_data$marker2,ord=1:nrow(ld_data))
    setkey(marker1,'marker')
    setkey(marker2,'marker')
    #------------------------------------------------
    if(verbose) pb <- txtProgressBar(style=3)
    rSet<-lapply(1:nrand,function(i){
      if(verbose) setTxtProgressBar(pb, i/nrand)
      snps<-popsnp[sample.int(n=nrow(popsnp),size=nvset),]
      snps<-snps[snps$chrom==chr,,drop=FALSE]
      if(nrow(snps)>0){
        rset<-list()
        rset[[chr]]<-lapply(1:nrow(snps),function(j){
          marker<-snps$rsid[j]
          markpos<-snps$position[j]
          ld<-ld_data[c(marker1[marker,nomatch=0]$ord,marker2[marker,nomatch=0]$ord),]
          if(nrow(ld)>0){
            tp1<-ld[ld$marker2==marker,c("marker1","position1"),drop=FALSE]
            tp2<-ld[ld$marker1==marker,c("marker2","position2"),drop=FALSE]
            linkedMarkers<-c(markpos,tp1$position1,tp2$position2)
            names(linkedMarkers)<-c(marker,tp1$marker1,tp2$marker2)
            linkedMarkers<-sort(linkedMarkers)
            return(linkedMarkers)
          } else {
            names(markpos)<-marker
            return(markpos)
          }
        })
        names(rset[[chr]])<-snps$rsid
      } else {
        rset<-NA
      }
      return(rset)
    })
    if(verbose) close(pb)
    #update randomSet
    if(length(randomSet)==0){
      randomSet<-rSet
    } else {
      junk<-sapply(1:length(randomSet),function(i){
        randomSet[[i]]<<-c(randomSet[[i]],rSet[[i]])
        NULL
      })
    }
  }
  #remove any chrom not selected!
  junk<-sapply(1:length(randomSet),function(i){
    tp<-randomSet[[i]]
    tp<-tp[!is.na(tp)]
    randomSet[[i]]<<-tp
    NULL
  })
  object@randomSet <- randomSet
  return(object)
}


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##                      --   AVS constructor    --
##                      --  LD1000gRel20130502  --
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##pre-processing for an avs
avs.preprocess.LD1000gRel20130502 <- function(markers, nrand=1000, mergeColinked=TRUE, 
                                           snpop="1000g", ld_window_kb=200, 
                                           ld_threshold_pval=1e-7, pAdjustMethod="bonferroni", 
                                           ld_threshold_r2=NULL, verbose=TRUE){
  
  #---check arguments
  if(missing(markers))stop("NOTE: 'markers' is missing!") 
  markers=avs.checks(name="markers",markers)
  avs.checks(name="nrand",para=nrand)
  avs.checks(name="mergeColinked",para=mergeColinked)
  snpop=avs.checks(name="snpop",para=snpop)
  avs.checks(name="LD1000gRel20130502_KBW",para=ld_window_kb)
  avs.checks(name="LD1000gRel20130502_pvalueR2",para=ld_threshold_pval)
  avs.checks(name="pAdjustMethod",para=pAdjustMethod)
  avs.checks(name="LD1000gRel20130502_R2",para=ld_threshold_r2)
  avs.checks(name="verbose",para=verbose)
  
  #---check available source for random AVS
  if(is.character(snpop)){
    opts<-c("dbSNP","1000g")
    if(!(snpop %in% opts))
      stop(paste("available 'snpop' options:", paste(opts,collapse = ", ") ),call.=FALSE)
  }
  
  #---initialize
  object <- .avs.initialize(markers,nrand,snpop,"RTNdata.LD1000gRel20130502")
  object <- .checkUniverse.LD1000gRel20130502(object, verbose=verbose)
  
  #---build AVS
  object <- .buildAVS.LD1000gRel20130502(object, ld_window_kb, ld_threshold_pval, pAdjustMethod,
                                         ld_threshold_r2, mergeColinked, verbose=verbose)
  
  #---build random AVS
  object <- .buildRandomAVS.LD1000gRel20130502(object, snpop, ld_window_kb, ld_threshold_pval, 
                                               pAdjustMethod, ld_threshold_r2, verbose=verbose)
  
  #---get IRanges, update status and return results
  object <- .avs.mapranges(object,verbose)
  return(object)
  
}

##------------------------------------------------------------------------------
## remove markers not listed in the "marker's universe" of the LD call
.checkUniverse.LD1000gRel20130502 <- function(object, verbose=TRUE){
  if(verbose){
    cat("Checking markers at 1000 Genomes (CEU population, release #20130502, GRCh37/hg19)...\n")
  }
  obj <- data("SNPs1000gRel20130502_A", package="RTNdata.LD1000gRel20130502", envir=environment())
  SNPs1000gRel20130502_A <- get(obj)
  idx1 <- object@validatedMarkers$rsid %in% SNPs1000gRel20130502_A
  rm(SNPs1000gRel20130502_A)
  obj <- data("SNPs1000gRel20130502_B", package="RTNdata.LD1000gRel20130502", envir=environment())
  SNPs1000gRel20130502_B <- get(obj)
  idx2 <- object@validatedMarkers$rsid %in% SNPs1000gRel20130502_B
  rm(SNPs1000gRel20130502_B)
  object@validatedMarkers <- object@validatedMarkers[idx1 | idx2, ]
  object@summary$markers[,"universe.removed"]<-nrow(object@validatedMarkers)
  return(object)
}


##------------------------------------------------------------------------------
## build AVS
.buildAVS.LD1000gRel20130502<-function(object, ld_window_kb, ld_threshold_pval, pAdjustMethod, 
                                       ld_threshold_r2, mergeColinked, verbose=TRUE){
  #---------
  reldata <- "RTNdata.LD1000gRel20130502"
  nCases <- 99 # n. cases in 1000 genome, Rel20130502, CEU population
  #---------
  mysnps <- object@validatedMarkers
  #set valid chroms
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(verbose){
    cat("Building AVS from 1000 Genomes LD data (CEU population, release #20130502, GRCh37/hg19)...\n")
  }
  variantSet<-list()
  for(chr in chrs){
    snps<-mysnps[mysnps$chrom==chr,,drop=FALSE]
    if(verbose) cat("...",nrow(snps)," marker(s) from ",chr,"!\n",sep="")
    if(nrow(snps)>0){
      ld_data<-getldata(chr,package=reldata) 
      #---apply thresholds
      ld_data <- ld_data[abs(ld_data$position1-ld_data$position2)<ld_window_kb*1000,]
      if(is.null(ld_threshold_r2)){
        pvalue <- 1 - pchisq(ld_data$R2*nCases,1)
        if(pAdjustMethod=="bonferroni"){
          ldw <- ld_window_kb*1000*2 #window around each marker!
          pvalue <- pvalue*ldw
        } else {
          pvalue <- p.adjust(pvalue, method=pAdjustMethod)
        }
        ld_data <- ld_data[pvalue<ld_threshold_pval,]
      } else {
        ld_data <- ld_data[ld_data$R2>ld_threshold_r2,]
      }
      #-------------------
      marker1<-data.table(marker=ld_data$marker1,ord=1:nrow(ld_data))
      marker2<-data.table(marker=ld_data$marker2,ord=1:nrow(ld_data))
      setkey(marker1,'marker')
      setkey(marker2,'marker')
      variantSet[[chr]]<-lapply(1:nrow(snps),function(i){
        marker<-snps$rsid[i]
        markpos<-snps$position[i]
        ld<-ld_data[c(marker1[marker,nomatch=0]$ord,marker2[marker,nomatch=0]$ord),]
        if(nrow(ld)>0){
          tp1<-ld[ld$marker2==marker,c("marker1","position1"),drop=FALSE]
          tp2<-ld[ld$marker1==marker,c("marker2","position2"),drop=FALSE]
          markpos<-c(markpos,tp1$position1,tp2$position2)
          names(markpos)<-c(marker,tp1$marker1,tp2$marker2)
          markpos<-sort(markpos)
        } else {
          names(markpos)<-marker  
        }
        markpos
      })
      names(variantSet[[chr]])<-snps$rsid
    }
  }
  #---check co-linked AVS
  if(mergeColinked){
    if(verbose)cat("Checking co-linked markers in the AVS...\n")
    variantSet<-check.colinked(variantSet, verbose=verbose)
    mkrs<-getMarkers.vset(variantSet,getlinked=FALSE)
    object@validatedMarkers<-object@validatedMarkers[mkrs,,drop=FALSE]
    object@summary$markers[,"colinked.removed"]<-nrow(object@validatedMarkers)
    if(verbose)cat("\n")
  }
  if(verbose)cat("\n")
  object@variantSet <- variantSet
  return(object)
}

##------------------------------------------------------------------------------
## build random AVS
.buildRandomAVS.LD1000gRel20130502<-function(object, snpop, ld_window_kb, ld_threshold_pval, 
                                             pAdjustMethod, ld_threshold_r2, verbose=TRUE){
  #---------
  reldata <- "RTNdata.LD1000gRel20130502"
  nCases <- 99 # n. cases in 1000 genome, Rel20130502, CEU population
  #---------
  variantSet <- object@variantSet
  nrand <- object@para$avs$nrand
  #set valid chroms
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(verbose){
    cat("Building matched random AVS from 1000 Genomes LD data (CEU population, release #20130502, GRCh37/hg19)...\n")
  }
  
  #count observed clusters in the variantSet
  clCount<-unlist(sapply(1:length(variantSet),function(i){
    unlist(lapply(variantSet[[i]],length))
  }))
  ncl1<-sum(clCount==1);ncl2<-sum(clCount>1)
  
  #get all markers in the variantSet
  allvset<-getMarkers.vset(variantSet,TRUE)
  
  #get popsnp
  if(is.data.frame(snpop)){
    #custom option: 'snp population' is defined by the user!
    popsnp<-snpop[snpop$chrom%in%chrs,]
    if(nrow(popsnp)<=(length(allvset)*2))
      stop("NOTE: not enough valid markers in the input 'snpop'!",call.=FALSE)
  } else {
    #default option: pre-computed 'snp population' based on a reference dataset (eg dbSNP, 1000g)
    if(snpop=='dbSNP'){
      obj=data('popsnp1', package=reldata, envir=environment())
      popsnp <- get(obj)
    } else if(snpop=='1000g'){
      obj=data('popsnp2', package=reldata, envir=environment())
      popsnp <- get(obj)
    }
    if(nrow(popsnp)<=(length(allvset)*2))
      stop("NOTE: not enough markers to build the random sets for the input data/options!",call.=FALSE)
  }
  
  #remove variantSet from popsnp
  popsnp<-popsnp[!popsnp$rsid%in%allvset,]
  
  #get risk markers only
  nvset<-length(getMarkers.vset(variantSet,FALSE))
  
  #---build randomSet -- by chroms, to speed the process
  if(verbose) cat("Computing random sets...\n")
  randomSet<-list()
  for(chr in chrs){
    if(verbose) cat("",chr,"\n",sep="")
    ld_data<-getldata(chr,package=reldata)
    #---apply thresholds
    ld_data <- ld_data[abs(ld_data$position1-ld_data$position2)<ld_window_kb*1000,]
    if(is.null(ld_threshold_r2)){
      pvalue <- 1 - pchisq(ld_data$R2*nCases,1)
      if(pAdjustMethod=="bonferroni"){
        ldw <- ld_window_kb*1000*2 #window around each marker!
        pvalue <- pvalue*ldw
      } else {
        pvalue <- p.adjust(pvalue, method=pAdjustMethod)
      }
      ld_data <- ld_data[pvalue<ld_threshold_pval,]
    } else {
      ld_data <- ld_data[ld_data$R2>ld_threshold_r2,]
    }
    #-------------------
    marker1<-data.table(marker=ld_data$marker1,ord=1:nrow(ld_data))
    marker2<-data.table(marker=ld_data$marker2,ord=1:nrow(ld_data))
    setkey(marker1,'marker')
    setkey(marker2,'marker')
    #------------------------------------------------
    if(verbose) pb <- txtProgressBar(style=3)
    rSet<-lapply(1:nrand,function(i){
      if(verbose) setTxtProgressBar(pb, i/nrand)
      snps<-popsnp[sample.int(n=nrow(popsnp),size=nvset),]
      snps<-snps[snps$chrom==chr,,drop=FALSE]
      if(nrow(snps)>0){
        rset<-list()
        rset[[chr]]<-lapply(1:nrow(snps),function(j){
          marker<-snps$rsid[j]
          markpos<-snps$position[j]
          ld<-ld_data[c(marker1[marker,nomatch=0]$ord,marker2[marker,nomatch=0]$ord),]
          if(nrow(ld)>0){
            tp1<-ld[ld$marker2==marker,c("marker1","position1"),drop=FALSE]
            tp2<-ld[ld$marker1==marker,c("marker2","position2"),drop=FALSE]
            linkedMarkers<-c(markpos,tp1$position1,tp2$position2)
            names(linkedMarkers)<-c(marker,tp1$marker1,tp2$marker2)
            linkedMarkers<-sort(linkedMarkers)
            return(linkedMarkers)
          } else {
            names(markpos)<-marker
            return(markpos)
          }
        })
        names(rset[[chr]])<-snps$rsid
      } else {
        rset<-NA
      }
      return(rset)
    })
    if(verbose) close(pb)
    #update randomSet
    if(length(randomSet)==0){
      randomSet<-rSet
    } else {
      junk<-sapply(1:length(randomSet),function(i){
        randomSet[[i]]<<-c(randomSet[[i]],rSet[[i]])
        NULL
      })
    }
  }
  #remove any chrom not selected!
  junk<-sapply(1:length(randomSet),function(i){
    tp<-randomSet[[i]]
    tp<-tp[!is.na(tp)]
    randomSet[[i]]<<-tp
    NULL
  })
  object@randomSet <- randomSet
  return(object)
}



##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##                          -- Supplements --
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## initialize an avs object
.avs.initialize <- function(markers,nrand,snpop,reldata){
  spflag<-ifelse(is.data.frame(snpop),"custom",snpop)
  object <- new("AVS", markers=markers)
  object@para$avs<-list(nrand=nrand,reldata=reldata,snpop=spflag)
  object@summary$para$avs[1,]<-c(nrand,reldata,spflag)
  object@summary$markers[,"input"]<-length(object@markers)
  object@validatedMarkers<-validateMarkers(object@validatedMarkers)
  object@summary$markers[,"valid"]<-nrow(object@validatedMarkers)
  if(nrow(object@validatedMarkers)<3)stop("not enough valid rs# markers!")
  object@validatedMarkers<-sortPosition(object@validatedMarkers)
  return(object)
}

##------------------------------------------------------------------------------
.avs.mapranges <- function(object,verbose){
  if(verbose)cat("-Mapping AVS to ranges of integer values...\n")
  object@variantSet<-getAvsRanges(object@variantSet)
  object@randomSet<-getRandomAvsRanges(object@randomSet, verbose=verbose)
  object@status["Preprocess"] <- "[x]"
  return(object)
}

##------------------------------------------------------------------------------
## get ld data from the RTNdata package
getldata<-function(chrom,package){
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(chrom==chrs[1]){data('ldatachr1', package=package, envir=environment());ldata<-get("ldatachr1")
  } else if(chrom==chrs[2]){data('ldatachr2', package=package, envir=environment());ldata<-get("ldatachr2")
  } else if(chrom==chrs[3]){data('ldatachr3', package=package, envir=environment());ldata<-get("ldatachr3")
  } else if(chrom==chrs[4]){data('ldatachr4', package=package, envir=environment());ldata<-get("ldatachr4")
  } else if(chrom==chrs[5]){data('ldatachr5', package=package, envir=environment());ldata<-get("ldatachr5")
  } else if(chrom==chrs[6]){data('ldatachr6', package=package, envir=environment());ldata<-get("ldatachr6")
  } else if(chrom==chrs[7]){data('ldatachr7', package=package, envir=environment());ldata<-get("ldatachr7")
  } else if(chrom==chrs[8]){data('ldatachr8', package=package, envir=environment());ldata<-get("ldatachr8")
  } else if(chrom==chrs[9]){data('ldatachr9', package=package, envir=environment());ldata<-get("ldatachr9")
  } else if(chrom==chrs[10]){data('ldatachr10', package=package, envir=environment());ldata<-get("ldatachr10")
  } else if(chrom==chrs[11]){data('ldatachr11', package=package, envir=environment());ldata<-get("ldatachr11")
  } else if(chrom==chrs[12]){data('ldatachr12', package=package, envir=environment());ldata<-get("ldatachr12")
  } else if(chrom==chrs[13]){data('ldatachr13', package=package, envir=environment());ldata<-get("ldatachr13")
  } else if(chrom==chrs[14]){data('ldatachr14', package=package, envir=environment());ldata<-get("ldatachr14")
  } else if(chrom==chrs[15]){data('ldatachr15', package=package, envir=environment());ldata<-get("ldatachr15")
  } else if(chrom==chrs[16]){data('ldatachr16', package=package, envir=environment());ldata<-get("ldatachr16")
  } else if(chrom==chrs[17]){data('ldatachr17', package=package, envir=environment());ldata<-get("ldatachr17")
  } else if(chrom==chrs[18]){data('ldatachr18', package=package, envir=environment());ldata<-get("ldatachr18")
  } else if(chrom==chrs[19]){data('ldatachr19', package=package, envir=environment());ldata<-get("ldatachr19")
  } else if(chrom==chrs[20]){data('ldatachr20', package=package, envir=environment());ldata<-get("ldatachr20")
  } else if(chrom==chrs[21]){data('ldatachr21', package=package, envir=environment());ldata<-get("ldatachr21")
  } else if(chrom==chrs[22]){data('ldatachr22', package=package, envir=environment());ldata<-get("ldatachr22")
  } else if(chrom==chrs[23]){data('ldatachrX', package=package, envir=environment());ldata<-get("ldatachrX")
  }
  return(ldata)
}

##------------------------------------------------------------------------------
## Check co-linked SNP blocks
check.colinked<-function(variantSet, verbose=TRUE){
  variantSet<-lapply(variantSet,function(vset){
    res<-list()
    remove<-NULL
    junk<-lapply(1:length(vset),function(i){
      nm<-names(vset[i])
      if(!nm%in%remove){
        tp<-NULL
        junk<-sapply(1:length(vset),function(j){
          if( nm %in% names(vset[[j]]) ){
            tp<<-c(tp,vset[[j]])
          }
        })
        tp<-tp[!duplicated(names(tp))]
        res[[nm]]<<-sort(tp)
        tp<-names(tp)
        names(tp)<-rep(nm,length(tp))
        remove<<-c(remove,tp)
      } else {
        if(verbose)cat("...marker ",nm," found in cluster ",names(remove)[which(remove==nm)],"!\n",sep="")
      }
    })
    res
  })
  return(variantSet)
}

##------------------------------------------------------------------------------
##validate rs# markers (use for both input 'markers' and eventual 'snpop')
validateMarkers<-function(markers){
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  markers<-markers[markers$chrom%in%chrs,]
  markers<-markers[complete.cases(markers),]
  if(is.null(markers$position))markers$position<-markers$start
  markers<-markers[!base::duplicated(markers$rsid),]
  b1<-length(grep( "^rs[0-9]+$", markers$rsid))!=length(markers$rsid)
  if ( b1 )stop("all makers should be prefixed with 'rs' (e.g. rs10490113)!")
  rownames(markers)<-markers$rsid
  return(markers)
}
validateSnpop<-function(snpop){
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  snpop<-snpop[snpop$chrom%in%chrs,]
  snpop<-snpop[complete.cases(snpop),]
  if(is.null(snpop$position))snpop$position<-snpop$start
  snpop<-snpop[!base::duplicated(snpop$rsid),]
  tp <- length(grep( "^rs[0-9]+$", snpop$rsid))
  tp <- (tp/length(snpop$rsid))*100
  if ( tp < 95)warning("makers in 'snpop' should be prefixed with 'rs' (e.g. rs10490113)!")
  return(snpop)
}

##------------------------------------------------------------------------------
## sort markers
sortPosition<-function(markers){
  #sort markers
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  sortmrks<-data.frame(stringsAsFactors=FALSE)
  for(chr in chrs){
    mrk<-markers[markers$chrom==chr,,drop=FALSE]
    if(nrow(mrk)>0){
      idx<-sort.list(mrk$position)
      sortmrks<-rbind(sortmrks,mrk[idx,,drop=FALSE])
    }
  }
  rownames(sortmrks)<-sortmrks$rsid
  sortmrks
}

##------------------------------------------------------------------------------
##This function is used for argument checking
avs.checks <- function(name, para) {
  if(name=="nrand") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nrand' should be an integer >=1 !",call.=FALSE)
  }
  else if(name=="LD1000gRel20130502_R2") {
    if(!is.null(para)){
      if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0.2 || para>1)
        stop("'LD1000gRel20130502_R2' should be numeric value >=0.2 and <=1.0!\n",call.=FALSE)
    }
  }
  else if(name=="LD1000gRel20130502_KBW") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || para>1000)
      stop("'LD1000gRel20130502_KBW' should be an integer value >=1 and <=1000 kb!\n",call.=FALSE)
  }
  else if(name=="ld_threshold_pval") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'ld_threshold_pval' should be an integer or numeric value >=0 and <=1  !",call.=FALSE)
  }
  else if(name=="mergeColinked") {
    if(!is.logical(para) || length(para)!=1)
      stop("'mergeColinked' should be a logical value!",call.=FALSE)
  }
  else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
       !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm','hochberg','hommel','bonferroni','BH','BY','fdr' and 'none'!",call.=FALSE)
  }
  else if(name=="verbose") {
    if(!is.logical(para) || length(para)!=1)
      stop("'verbose' should be a logical value!",call.=FALSE)
  } else if(name=="markers") {
    if( !is.data.frame(para) || ncol(para)<4 ){
      stop("'markers' should be a dataframe, a 'BED file' format with ncol >= 4 !",call.=FALSE)
    }
    para<-para[,1:4]
    b1 <- !is.numeric(para[,2]) && !is.integer(para[,2])
    b2 <- !is.numeric(para[,3]) && !is.integer(para[,3])
    if(b1 || b2){
      stop("'markers' should have a 'BED file' format, with chromosomal positions as integer values!",call.=FALSE)
    }
    para$start<-as.integer(para$start)
    para$end<-as.integer(para$end)
    colnames(para)<-c("chrom","start","end","rsid")
    if(is.numeric(para$chrom) || is.integer(para$chrom)){
      para$chrom <- paste("chr",para$chrom,sep="")
    }
    para$chrom<-as.character(para$chrom)
    para$rsid<-as.character(para$rsid)
    return(para)
  } else if(name=="snpop") {
    if(is.character(para)){
      if(length(para)!=1)
        stop("'snpop' should be a single string!",call.=FALSE)
    } else if(is.data.frame(para)){
      if(all(colnames(para)%in%c("rsid","chrom","position"))){
        if(!is.numeric(para$position) || !is.integer(para$position)){
          stop("chromosomal positions in 'snpop' should be numerical or integer values!",call.=FALSE)
        }
      } else if(ncol(para)<4 ){
        stop("'snpop' should be a dataframe, a 'BED file' format with ncol >= 4 !",call.=FALSE)
      } else {
        para<-para[,1:4]
        colnames(para)<-c("chrom","start","end","rsid")
        if(!is.numeric(para$start) || !is.integer(para$start)){
          stop("chromosomal positions in 'snpop' should be numerical or integer values!",call.=FALSE)
        }
      }
      if(is.numeric(para$chrom) || is.integer(para$chrom)){
        para$chrom <- paste("chr",para$chrom,sep="")
      }
      para$chrom<-as.character(para$chrom)
      para$rsid<-as.character(para$rsid)
      para<-validateSnpop(para)
    } else {
      stop("'snpop' should be a dataframe, a 'BED file' format with ncol = 4 !",call.=FALSE)
    }
    return(para)
  }
  
}
