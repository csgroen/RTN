

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##                      -- AVS supplements --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##-------------------------------------------------------------------------
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

##-------------------------------------------------------------------------
## sort annotation
sortAnnotation<-function(annotation){
  #sort annotation
  chrs<-c(paste("chr",1:22,sep=""),"chrX","chrY")
  sortannot<-data.frame(stringsAsFactors=FALSE)
  for(chr in chrs){
    annot<-annotation[annotation$chrom==chr,,drop=FALSE]
    if(nrow(annot)>0){
      idx<-sort.list(annot$start)
      sortannot<-rbind(sortannot,annot[idx,,drop=FALSE])
    }
  }
  rownames(sortannot)<-NULL
  sortannot
}

##-------------------------------------------------------------------------
## build AVS
buildAVS<-function(mysnps, reldata="rel27CEU-NCBIB36", ldfilter="DprimeLOD", verbose=TRUE){
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(reldata=="rel27CEU-NCBIB36"){
    if(verbose) cat("Building AVS from HapMap LD data (CEU, rel #27, NCBI B36)...\n")
  }
  variantSet<-list()
  for(chr in chrs){
    snps<-mysnps[mysnps$chrom==chr,,drop=FALSE]
    if(verbose) cat("...",nrow(snps)," marker(s) from ",chr,"!\n",sep="")
    if(nrow(snps)>0){
      ld_data<-getldata(chr,ldfilter)
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
  return(variantSet)
}

##-------------------------------------------------------------------------
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

##-------------------------------------------------------------------------
## build random AVS
#for(i in 1:length(variantSet))print(length(variantSet[[i]]))
buildRandomAVS<-function(variantSet, nrand=100, reldata="rel27CEU-NCBIB36", ldfilter="DprimeLOD", snpop=NULL, verbose=TRUE){
  #set valid chroms
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  #get HapMap
  if(reldata=="rel27CEU-NCBIB36"){
    if(verbose) cat("Building random AVS from HapMap LD data (CEU, rel #27, NCBI B36)...\n")
    data('popsnps.hg18', envir=environment())
    popsnps<-get("popsnps")
    popsnps<-popsnps[popsnps$chrom%in%chrs,]
  }
  #get all markers from the variantSet
  vset<-getMarkers.vset(variantSet,TRUE)
  #filter popsnps
  if(!is.null(snpop)){
    popsnps<-popsnps[popsnps$rsid%in%snpop,]
    if(nrow(popsnps)<=(length(vset)*2)){
      stop("NOTE: not enough valid markers in 'snpop'!",call.=FALSE)
    } else if(nrow(popsnps)<length(snpop)){
      n<-length(snpop)-nrow(popsnps)
      warning("NOTE: ",n," marker(s) from 'snpop' not represented in the '",reldata,"' LD dataset!")
    }
  }
  #remove variantSet from popsnp
  popsnps<-popsnps[!popsnps$rsid%in%vset,]
  #get risk markers
  vset<-getMarkers.vset(variantSet,FALSE)
  #get random markers
  randmarkers<-sapply(1:nrand,function(i){
    sample.int(n=nrow(popsnps),size=length(vset))
  })
  #---build randomSet by chrom
  randomSet<-list()
  for(chr in chrs){
    if(verbose) cat("",chr,"\n",sep="")
    ld_data<-getldata(chr,ldfilter)
    marker1<-data.table(marker=ld_data$marker1,ord=1:nrow(ld_data))
    marker2<-data.table(marker=ld_data$marker2,ord=1:nrow(ld_data))
    setkey(marker1,'marker')
    setkey(marker2,'marker')
    #------------------------------------------------
    if(verbose) pb <- txtProgressBar(style=3)
    rSet<-lapply(1:nrand,function(i){
      if(verbose) setTxtProgressBar(pb, i/nrand)
      snps<-popsnps[randmarkers[,i],]
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
  return(randomSet)
}

##------------------------------------------------------------------------
##get markers from computed AVS
getMarkers.vset<-function(variantSet,getlinked=TRUE){
  markers<-NULL
  lkmarkers<-unique(unlist(
    lapply(variantSet,function(vset){
      if(class(vset)=="IRanges"){
        markers<<-c(markers,names(vset@metadata$blocks))
        tp<-lapply(vset@metadata$blocks,names)
      } else {
        markers<<-c(markers,names(vset))
        tp<-lapply(vset,names)
      }
      tp
    })
  ))
  if(getlinked){
    return(lkmarkers)
  } else {
    return(markers)
  }
}

##------------------------------------------------------------------------
##get markers from computed random AVS
getMarkers.rset<-function(randomSet,getlinked=TRUE){
  lkmarkers<-lapply(1:length(randomSet),function(i){
    res<-getMarkers.vset(randomSet[[i]],getlinked)
  })
  unique(unlist(lkmarkers))
}

##-------------------------------------------------------------------------
##map AVS to snpdate (speed-up the permutation step)
mapvset<-function(vSet,snpnames){
  snpnames <- data.table(x=snpnames, y=1:length(snpnames),key="x")
  y=NULL
  junk<-lapply(1:length(vSet),function(i){
    tp<-metadata(vSet[[i]])
    mappedMarkers<-list()
    junk<-lapply(1:length(tp$blocks),function(j){
      tpp<-tp$blocks[[j]]
      mapm<-snpnames[data.table(names(tpp))][,y]
      mapm<-mapm[!is.na(mapm)]
      mappedMarkers[[names(tp$blocks[j])]]<<-mapm
      NULL
    })
    vSet[[i]]@metadata$mappedMarkers<<-mappedMarkers
    NULL
  })
  vSet
}

##-------------------------------------------------------------------------
##map ramdom AVS to snpdate (speed-up the permutation step)
maprset<-function(rSet,snpnames,verbose=TRUE){
  snpnames <- data.table(x=snpnames, y=1:length(snpnames),key="x")
  nr<-length(rSet)
  if(verbose) pb <- txtProgressBar(style=3)
  y=NULL
  resrset<-lapply(1:nr,function(i){
    vSet<-rSet[[i]]
    if(verbose) setTxtProgressBar(pb, i/nr)
    junk<-lapply(1:length(vSet),function(i){
      tp<-metadata(vSet[[i]])
      mappedMarkers<-list()
      junk<-lapply(1:length(tp$blocks),function(j){
        tpp<-tp$blocks[[j]]
        mapm<-snpnames[data.table(names(tpp))][,y]
        mapm<-mapm[!is.na(mapm)]
        mappedMarkers[[names(tp$blocks[j])]]<<-mapm
        NULL
      })
      vSet[[i]]@metadata$mappedMarkers<<-mappedMarkers
      NULL
    })
    vSet
  })
  resrset
}

##-------------------------------------------------------------------------
##get rs# markers from internal dataset
getmarkers<-function(markers,reldata="rel27CEU-NCBIB36"){
  if(reldata=="rel27CEU-NCBIB36"){
    data('popsnps.hg18', envir=environment())
    popsnps<-get("popsnps")
    markers<-popsnps[popsnps$rsid%in%markers,,drop=FALSE]
  } else {
    #eventually, other LD dataset!
  }
  return(markers)
}

##-------------------------------------------------------------------------
##get IRanges for the AVS
getAvsRanges<-function(vSet, continued=FALSE){
  clustersRanges<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    blocks<-vSet[[i]]
    clRanges<-IRanges()
    index<-NULL
    junk<-lapply(1:length(blocks),function(j){
      pos<-as.integer(blocks[[j]])
      if(continued){
        query<-IRanges(start=min(pos), end=max(pos),names=names(blocks)[j])
        index<<-c(index,j)
      } else {
        query<-IRanges(start=pos, end=pos,names=names(blocks[[j]]))
        index<<-c(index,rep(j,length(query)))
      }
      clRanges<<-c(clRanges,query)
    })
    #if(getTree)clRanges<-IntervalTree(clRanges)
    clRanges@metadata$chr<-chr
    clRanges@metadata$markers<-names(blocks)
    clRanges@metadata$blocks<-blocks
    clRanges@metadata$index<-index
    clRanges
  })
  names(clustersRanges)<-names(vSet)
  return(clustersRanges)
}
##get IRanges for the random AVS
getRandomAvsRanges<-function(rSet, continued=FALSE, verbose=TRUE){
  if(verbose) pb <- txtProgressBar(style=3)
  clustersRanges<-lapply(1:length(rSet),function(i){
    if(verbose) setTxtProgressBar(pb, i/length(rSet)) 
    getAvsRanges(rSet[[i]], continued)
  })
  if(verbose)close(pb)
  clustersRanges
}
##coerce avs: IRanges to IntervalTree
coerceAvsRanges<-function(vSet){
  clustersRanges<-sapply(1:length(vSet),function(i){
    clRanges<-IntervalTree(vSet[[i]])
    clRanges@metadata$chr<-vSet[[i]]@metadata$chr
    clRanges@metadata$markers<-vSet[[i]]@metadata$markers
    clRanges@metadata$blocks<-vSet[[i]]@metadata$blocks
    clRanges@metadata$index<-vSet[[i]]@metadata$index
    clRanges
  })
  names(clustersRanges)<-names(vSet)
  return(clustersRanges)
}
##coerce random avs: IRanges to IntervalTree
coerceRandomAvsRanges<-function(rSet,verbose=TRUE){
  if(verbose) pb <- txtProgressBar(style=3)
  clustersRanges<-lapply(1:length(rSet),function(i){
    if(verbose) setTxtProgressBar(pb, i/length(rSet)) 
    coerceAvsRanges(rSet[[i]])
  })
  if(verbose)close(pb)
  clustersRanges
}
##get IRanges/IntervalTree for annotation
getAnnotRanges<-function(annotation,maxgap=0, getTree=TRUE, 
                         getReduced=FALSE){
  annot<-annotation
  idx<-annot$START>annot$END
  annot$START[idx]<-annotation$END
  annot$END[idx]<-annotation$START
  annotation<-annot
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  chrs<-chrs[chrs%in%unique(annotation$CHROM)]
  annotRange<-sapply(chrs,function(chr){
    annot<-annotation[annotation$CHROM==chr,c("ID","START","END")]
    start<-as.integer(annot$START-maxgap)
    start[start<0]<-0
    end<-as.integer(annot$END+maxgap)
    subject<-IRanges(start, end, names=annot$ID)
    if(getReduced)subject<-reduce(subject)
    if(getTree)subject<-IntervalTree(subject)
    subject@metadata<-list(mappedAnnotations=annot$ID)
    subject
  })
  annotRange
}

##------------------------------------------------------------------------
##get markers and genes from computed evse
##return a named character vector with the RiskAssociatedSNPs mapped in the VSE analysis
##..names indicate the LD cluster (i.e. the RiskSNP assignment) 
getMappedClusters<-function(object){
  MarkerIDs<-NULL
  junk<-lapply(object@variantSet,function(vset){
    tp<-lapply(vset@metadata$blocks,names)
    tp<-unlist(tp)
    names(tp)<-vset@metadata$markers[vset@metadata$index]
    MarkerIDs<<-c(MarkerIDs,tp)
    NULL
  })
  mappedIds<-getMappedMarkers(object)
  MarkerIDs<-MarkerIDs[names(MarkerIDs)%in%mappedIds$RiskSNP]
  MarkerIDs<-MarkerIDs[MarkerIDs%in%mappedIds$RiskAssociatedSNP]
  return(MarkerIDs)
}
##return a list with RiskSNPs and RiskAssociatedSNPs mapped in the VSE analysis
getMappedMarkers<-function(object){
  RiskSNP<-NULL
  RiskAssociatedSNP<-NULL
  junk<-sapply(names(object@results$evse),function(reg){
    eqtls<-object@results$evse[[reg]]$eqtls
    RiskSNP<<-c(RiskSNP,eqtls$RiskSNP)
    RiskAssociatedSNP<<-c(RiskAssociatedSNP,eqtls$RiskAssociatedSNP)
    NULL
  })
  rsid<-object@validatedMarkers$rsid
  RiskSNP<-rsid[rsid%in%RiskSNP]
  RiskAssociatedSNP<-unique(RiskAssociatedSNP)
  return(list(RiskSNP=RiskSNP,RiskAssociatedSNP=RiskAssociatedSNP))
}
getEvseMatrix<-function(object,what="probs"){
  mappedIds<-getMappedMarkers(object)
  RiskSNP<-mappedIds$RiskSNP
  evsemtx<-NULL
  if(what=="probs"){
    vl=1;cl="Pr(>F)";efun=min
  } else if(what=="fstat"){
    vl=0;cl="F";efun=max
  } else if(what=="coef") {
    vl=0;cl=".R2";efun=function(x)x[which.max(abs(x))]
  }
  junk<-sapply(names(object@results$evse),function(reg){
    eqtls<-object@results$evse[[reg]]$eqtls
    rvec<-rep(vl,length(RiskSNP))
    names(rvec)<-RiskSNP
    junk<-sapply(RiskSNP,function(rs){
      idx<-which(eqtls$RiskSNP==rs)
      if(length(idx)>0){
        rvec[rs]<<-efun(eqtls[idx,cl])
      }
      NULL
    })
    evsemtx<<-rbind(evsemtx,rvec)
    NULL
  })
  rownames(evsemtx)<-names(object@results$evse)
  evsemtx
}


###########################################################################
## VSE analysis
###########################################################################
##-------------------------------------------------------------------------
vsea<-function(vSet,rSet,annot,verbose=TRUE){
  #compute vse
  resvset<-get.avsdist(vSet=vSet,annot=annot)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    #snow::clusterEvalQ(cl, expr=library("IRanges"))
    snow::clusterExport(cl, list("get.avsdist","IRanges","overlapsAny","ff"),
                        envir=environment())
    resrset<-parSapply(cl, 1:length(rSet), function(i) {
      sum(get.avsdist(rSet[[i]],annot))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet)) 
      sum(get.avsdist(rSet[[i]],annot))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=resvset,nulldist=resrset,nclusters=length(resvset)))
}

##-------------------------------------------------------------------------
##get avs overlap for a variant set
get.avsdist<-function(vSet,annot){
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]   
    if(!is.null(subject)){
      res<-rep(FALSE,length(query@metadata$markers))
      ov<-query@metadata$index[overlapsAny(query,subject)]
      res[ov]<-TRUE
    } else {
      res<-rep(FALSE,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    res
  })
  resvset<-unlist(clusterMapping)
  return(resvset)
}

##-------------------------------------------------------------------------
##get annotation overlap
##annot should be named!
getAnnotOverlap<-function(vSet,annot){
  annotMapping<-lapply(1:length(annot),function(i){
    chr<-names(annot[i])
    query<-annot[[i]]
    subject<-vSet[[chr]]  
    if(!is.null(subject)){
      res<-overlapsAny(query,subject)
    } else {
      res<-rep(FALSE,length(query))
    }
    names(res)<-names(query)
    res
  })
  annotMapping<-unlist(annotMapping)
  return(annotMapping)
}

###########################################################################
## EVSE analysis
###########################################################################

##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
evsea<-function(vSet, rSet, annot, gxdata, snpdata, pValueCutoff=0.01,verbose=TRUE){
  #compute evse
  mtally<-get.eqtldist(vSet, annot, gxdata, snpdata, pValueCutoff)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    #snow::clusterEvalQ(cl, library("ff"))
    #snow::clusterEvalQ(cl, library("IRanges"))
    snow::clusterExport(cl, list("get.eqtldist","eqtlTest","IRanges","overlapsAny","findOverlaps","ff"),
                        envir=environment())
    resrset<-parSapply(cl, 1:length(rSet), function(i) {
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
    if(verbose)close(pb)
  }
  return(list(mtally=mtally,nulldist=resrset,nclusters=length(mtally)))
}

##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
evseaproxy<-function(rSet, annot, gxdata, snpdata, pValueCutoff=0.01,verbose=TRUE){
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    #snow::clusterEvalQ(cl, library("ff"))
    #snow::clusterEvalQ(cl, library("IRanges"))
    snow::clusterExport(cl, list("get.eqtldist","eqtlTest","IRanges","overlapsAny","findOverlaps","ff"),
                        envir=environment())
    nulldist<-parSapply(cl, 1:length(rSet), function(i) {
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    nulldist<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
    if(verbose)close(pb)
  }
  return(nulldist)
}

##-------------------------------------------------------------------------
##get avs/eqtl dist for a variant set
get.eqtldist<-function(vSet,annot,gxdata,snpdata,pValueCutoff=0.01){
  # mapping tally
  clusterMapping<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-metadata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        snpList<-metadata(query)$mappedMarkers[[j]]
        ov<-overlaps@queryHits%in%which(query@metadata$index==j)
        if(any(ov) && length(snpList)>0){
          ov<-unique(overlaps@subjectHits[ov])
          gList<-geneList[ov]
          if(length(gList)>0){
            res<-eqtlTest(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          } else {
            res<-1.0
          }
        } else {
          res<-1.0
        }
        return(res)
      })
    } else {
      res<-rep(1.0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  clusterMapping<-unlist(clusterMapping)
  clusterMapping<pValueCutoff
}

##-------------------------------------------------------------------------
##eqtl test
##run two-way manova with multiple additive factors (gx as response varible)
eqtlTest<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  # set names for formulae
  colnames(gxdata)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(snpdata)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  # run lm
  fm1<-paste(colnames(snpdata),collapse="+")
  fmla <- formula( paste("gxdata ~",fm1, collapse=" ") )
  resfit<-lm(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
  if(ncol(gxdata)>1){
    resf<-summary(manova(resfit))
    resf<-resf$stats[,"Pr(>F)"]
  } else {
    resf<-anova(resfit) 
    resf<-resf[["Pr(>F)"]]
  }
  if(sum(is.na(resf))==length(resf)){
    res=1
  } else {
    res<-min(resf,na.rm=TRUE)
  }
  return(res)
}

#-------------------------------------------------------------------------
#extract eqtls, return consolidated results
# eqtlExtract<-function(vSet,annot,gxdata,snpdata){ 
#   # mapping tally
#   clusterMapping<-lapply(1:length(vSet),function(i){
#     chr<-names(vSet[i])
#     query<-vSet[[i]]
#     subject<-annot[[chr]]
#     if(!is.null(subject)){
#       overlaps<-findOverlaps(query,subject)
#       geneList<-metadata(subject)$mappedAnnotations
#       res<-lapply(1:length(query@metadata$markers),function(j){
#         snpList<-metadata(query)$mappedMarkers[[j]]
#         ov<-overlaps@queryHits%in%which(query@metadata$index==j)
#         ov<-unique(overlaps@subjectHits[ov])
#         gList<-geneList[ov]
#         if(length(gList)>0 && length(snpList)>0){
#           res<-eqtlTestDetailed(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
#           rownames(res)<-rownames(snpdata)[snpList]
#         } else {
#           res<-matrix(NA,nrow=length(snpList),ncol=2)
#           rownames(res)<-snpList;colnames(res)<-c("F","Pr(>F)")
#         }
#         return(res)
#       })
#       names(res)<-query@metadata$markers
#     } else {
#       res<-NA
#     }
#     return(res)
#   })
#   #---simplify list
#   res<-list()
#   lapply(1:length(clusterMapping),function(i){
#     tp<-clusterMapping[[i]]
#     if(is.list(tp)){
#       lapply(names(tp),function(j){
#         tpp<-tp[[j]]
#         tpp<-tpp[!is.na(tpp[,2]),]
#         if(ncol(tpp)>2 && nrow(tpp)>0){
#           res[[j]]<<-tpp
#         }
#       })
#     }
#     NULL
#   })
#   clusterMapping<-res
#   #---get summary
#   summ<-data.frame(NULL,stringsAsFactors=FALSE)
#   lapply(names(clusterMapping),function(riskSNP){
#     tp<-clusterMapping[[riskSNP]]
#     Fstat<-tp[,1,drop=FALSE]
#     Pstat<-tp[,2,drop=FALSE]
#     sapply(rownames(Fstat),function(associatedSNP){
#       tpp<-data.frame(riskSNP,associatedSNP,Fstat[associatedSNP,],Pstat[associatedSNP,],stringsAsFactors=FALSE)
#       summ<<-rbind(summ,tpp)
#       NULL
#     })
#   })
#   if(nrow(summ)>0){
#     colnames(summ)<-c("RiskSNP","RiskAssociatedSNP","F","Pr(>F)")
#   }
#   summ
# }

#-------------------------------------------------------------------------
#other table formats to extract eqtls, return consolidated results
eqtlExtract<-function(vSet,annot,gxdata,snpdata,pValueCutoff){
  eqtls<-eqtlExtractFull(vSet,annot,gxdata,snpdata)
  eqtls<-eqtls[eqtls$'Pr(>F)'<pValueCutoff,]
  rownames(eqtls)<-NULL
  eqtls
}
eqtlExtractFull<-function(vSet,annot,gxdata,snpdata){ 
  # mapping tally
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-metadata(subject)$mappedAnnotations
      res<-lapply(1:length(query@metadata$markers),function(j){
        snpList<-metadata(query)$mappedMarkers[[j]]
        ov<-overlaps@queryHits%in%which(query@metadata$index==j)
        ov<-unique(overlaps@subjectHits[ov])
        gList<-geneList[ov]
        if(length(gList)>0 && length(snpList)>0){
          res<-eqtlTestDetailed(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          rownames(res)<-rownames(snpdata)[snpList]
        } else {
          res<-matrix(NA,nrow=length(snpList),ncol=2)
          rownames(res)<-snpList;colnames(res)<-c("F","Pr(>F)")
        }
        return(res)
      })
      names(res)<-query@metadata$markers
    } else {
      res<-NA
    }
    return(res)
  })
  #---simplify list
  res<-list()
  lapply(1:length(clusterMapping),function(i){
    tp<-clusterMapping[[i]]
    if(is.list(tp)){
      lapply(names(tp),function(j){
        tpp<-tp[[j]]
        tpp<-tpp[!is.na(tpp[,2]),]
        if(ncol(tpp)>2 && nrow(tpp)>0){
          res[[j]]<<-tpp
        }
      })
    }
    NULL
  })
  clusterMapping<-res
  #---get summary
  summ<-data.frame(NULL,stringsAsFactors=FALSE)
  lapply(names(clusterMapping),function(riskSNP){
    tp<-clusterMapping[[riskSNP]]
    Fstat<-tp[,1,drop=FALSE]
    Pstat<-tp[,2,drop=FALSE]
    geneid<-colnames(tp)[-c(1,2)]
    sapply(rownames(Fstat),function(associatedSNP){
      tpp<-data.frame(riskSNP,associatedSNP,geneid,Fstat[associatedSNP,],Pstat[associatedSNP,],stringsAsFactors=FALSE)
      summ<<-rbind(summ,tpp)
      NULL
    })
  })
  if(nrow(summ)>0){
    colnames(summ)<-c("RiskSNP","RiskAssociatedSNP","GeneID","F","Pr(>F)")
  }
  summ
}


#-------------------------------------------------------------------------
##eqtl test for 'eqtlExtract' functions
##run two-way manova with multiple additive factors (gx as response varible)
eqtlTestDetailed<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  #set names for formulae
  names(snpList)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  names(geneList)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(gxdata)<-names(geneList)
  colnames(snpdata)<-names(snpList)
  #run lm
  fm1<-paste(colnames(snpdata),collapse="+")
  fmla <- formula( paste("gxdata ~",fm1, collapse=" ") )
  resfit<-lm(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
  if(ncol(gxdata)>1){
    resf<-summary(manova(resfit))
    resf<-as.data.frame(resf$stats)
    resf<-resf[1:(nrow(resf)-1),]
    resf<-resf[,c("approx F","Pr(>F)")]
  } else {
    resf<-anova(resfit)
    resf<-resf[1:(nrow(resf)-1),]
    resf<-as.data.frame(resf)
    resf<-resf[,c("F value","Pr(>F)")]
  }
  #library('gplots')
  #plotmeans(G2 ~ S6, data=cbind(gxdata,snpdata), col="red", barcol="blue",connect=FALSE,pch=15, las=1)
  #get SNP stats
  colnames(resf)<-c("F","Pr(>F)")
  #get gene stats
  resaov<-summary(aov(resfit))
  resp<-sapply(1:length(resaov),function(i){
    tp<-resaov[[i]][["Pr(>F)"]]
    tp[-length(tp)]
  })
  resp<-matrix(resp,ncol=length(geneList))
  colnames(resp)<-geneList
  #combine results
  res<-cbind(resf,resp)
  res<-res[names(snpList),]
  rownames(res)<-snpList
  return(res)
}

#-------------------------------------------------------------------------
eqtlExtractAnova<-function(vSet,annot,gxdata,snpdata){ 
  # mapping tally
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-metadata(subject)$mappedAnnotations
      resfit<-NULL
      junk<-lapply(1:length(query@metadata$markers),function(j){
        snpList<-metadata(query)$mappedMarkers[[j]]
        ov<-overlaps@queryHits%in%which(query@metadata$index==j)
        ov<-unique(overlaps@subjectHits[ov])
        gList<-geneList[ov]
        if(length(gList)>0 && length(snpList)>0){
          res<-eqtlTestDetailedAnova(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          res$RiskAssociatedSNP<-rownames(snpdata)[res$RiskAssociatedSNP]
          res<-data.frame(RiskSNP=query@metadata$markers[j],res,stringsAsFactors=FALSE)
          resfit<<-rbind(resfit,res)
        }
        NULL
      })
    } else {
      resfit<-NA
    }
    return(resfit)
  })
  #---simplify list
  summ<-NULL
  lapply(1:length(clusterMapping),function(i){
    summ<<-rbind(summ,clusterMapping[[i]])
    NULL
  })
  if(nrow(summ)>0){
    colnames(summ)<-c(".RiskSNP",".RiskAssociatedSNP",".GeneID",".F",".Pr(>F)",".Coef",".R2")
  }
  summ
}

#-------------------------------------------------------------------------
eqtlTestDetailedAnova<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  #set names for formulae
  names(snpList)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  names(geneList)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(gxdata)<-names(geneList)
  colnames(snpdata)<-names(snpList)
  #run lm
  resfit<-NULL
  if(ncol(snpdata)>1 && ncol(gxdata)>1){
    junk<-sapply(colnames(snpdata),function(S){
      fmla <- formula( paste("gxdata ~",S, collapse=" ") )
      raov <- aov(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
      cf<-raov$coefficients[S,]
      sf<-sapply(summary(raov),function(rv){
        tp<-as.data.frame(rv)
        tp<-tp[1:(nrow(tp)-1),4:5]
        as.numeric(tp)
      })
      rownames(sf)<-c("F","Prob")
      colnames(sf)<-names(cf)
      sf<-t(sf)
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=rownames(sf),sf,Coef=cf,stringsAsFactors=FALSE)
      resfit<<-rbind(resfit,raov)
    })
  } else {
    junk<-sapply(colnames(snpdata),function(S){
      fmla <- formula( paste("gxdata ~",S, collapse=" ") )
      raov <- aov(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
      cf<-raov$coefficients[2]
      sf<-as.data.frame(summary(raov)[[1]])
      sf<-sf[1,4:5]
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=colnames(gxdata),F=sf[,1],Prob=sf[,2],Coef=cf,stringsAsFactors=FALSE)
      resfit<<-rbind(resfit,raov)
    })
  }
  rownames(resfit)<-NULL
  resfit$RiskAssociatedSNP<-snpList[resfit$RiskAssociatedSNP]
  resfit$GeneID<-geneList[resfit$GeneID]
  #--
  #run cor
  R2<-cor(gxdata,snpdata)
  colnames(R2)<-snpList
  rownames(R2)<-geneList
  summ<-NULL
  sapply(colnames(R2),function(i){
    tp<-data.frame(RiskAssociatedSNP=i,GeneID=rownames(R2),R2=R2[,i],stringsAsFactors=FALSE)
    summ<<-rbind(summ,tp)
    NULL
  })
  rownames(summ)<-NULL
  #---
  resfit<-cbind(resfit,R2=summ$R2)
  return(resfit)
}

#-------------------------------------------------------------------------
getUniverseCounts1<-function(vSet,annotation,maxgap){
  #count markers in hapmap
  clusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE, getReduced=FALSE)
  # mapping tally
  geneCounts<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      geneList<-metadata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        gList<-geneList[overlapsAny(subject,query[query@metadata$index==j])]
        length(gList)
      })
    } else {
      res<-rep(0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  geneCounts<-unlist(geneCounts)
  #merge counts
  counts<-t(rbind(clusterCounts,geneCounts))
  colnames(counts)<-c("markers","annotation")
  counts
}

#-------------------------------------------------------------------------
getUniverseCounts2<-function(vSet,annotation,maxgap){
  #count markers in hapmap
  clusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count markers mapped to the genotype data
  clusterCounts<-rbind(clusterCounts,unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$mappedMarkers,length))
  })))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation=annotation,maxgap=maxgap,getTree=FALSE, getReduced=FALSE)
  # mapping tally
  geneCounts<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      geneList<-metadata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        gList<-geneList[overlapsAny(subject,query[query@metadata$index==j])]
        length(gList)
      })
    } else {
      res<-rep(0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  geneCounts<-unlist(geneCounts)
  #merge counts
  counts<-t(rbind(clusterCounts,geneCounts))
  colnames(counts)<-c("totalMarkers","markers","annotation")
  counts
}

#-------------------------------------------------------------------------
vseformat<-function(resavs, pValueCutoff=0.01, boxcox=TRUE){
  groups<-rep(1,length(resavs)) #'groups' nao ativo!
  ntests<-length(resavs)
  #get mtally
  mtally<-sapply(names(resavs),function(i){
    resavs[[i]]$mtally
  })
  #get mtally (transformed/normalized if possible)
  nulldist<-sapply(names(resavs),function(i){
    tp<-c(resavs[[i]]$nulldist,sum(resavs[[i]]$mtally))
    if(median(tp)>0){
      tp<-(tp-median(tp))/sd(tp)
    }
    tp
  })
  escore<-nulldist[nrow(nulldist),]
  nulldist<-nulldist[-nrow(nulldist),,drop=FALSE]
  #----get stats
  pvals <- pnorm(escore, lower.tail=FALSE)
  if(is.null(ntests))ntests<-length(pvals)
  ci <- qnorm(1-(pValueCutoff/ntests))
  #----powerTransform
  if(boxcox){
    ptdist<-sapply(1:ncol(nulldist),function(i){
      ndat<-nulldist[,i]
      nnull<-length(ndat)
      edat<-escore[i]
      nd<-as.integer(nnull*.05)
      nd<-max(nd,10)
      qt<-quantile(ndat,probs=seq(0,1,length.out=nd),names=FALSE)
      checkd<-!all(ndat>=0) && min(qt)<max(qt)
      if(checkd && shapiro.test(qt)$p.value<0.05){
        minval<-min(c(nulldist[,i],escore[i]))
        minval<-ifelse(minval<=0,abs(minval)+1,minval)
        ndat<-ndat+minval
        edat<-edat+minval
        l<-coef(powerTransform(c(ndat,edat)), round=TRUE)
        ptdat<-bcPower(c(ndat,edat),l)
        ptdat<-(ptdat-median(ptdat))/sd(ptdat)
        return(ptdat)
      } else {
        return(c(ndat,edat))
      }
    })
    colnames(ptdist)<-names(escore)
    escore<-ptdist[nrow(ptdist),]
    nulldist<-ptdist[-nrow(ptdist),,drop=FALSE]
    pvals <- pnorm(escore, lower.tail=FALSE)
  }
  #----order
  if(length(groups)>1){
    gs<-unique(groups)
    ord<-lapply(1:length(gs),function(i){
      idx<-sort.list(escore[groups==i],decreasing=TRUE)
      labs<-names(escore)[groups==i]
      labs[idx]
    })
    ord<-unlist(ord)
    escore<-escore[ord]
    pvals<-pvals[ord]
    mtally<-mtally[,ord]
    nulldist<-nulldist[,ord]
  }
  return(list(mtally=mtally,nulldist=nulldist,escore=escore,pvalue=pvals,ci=ci))
}

##-------------------------------------------------------------------------
##check if snow cluster is loaded
isParallel<-function(){
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
  all(c(b1,b2))
}

##-------------------------------------------------------------------------
## get ld data from RTNdata
getldata<-function(chrom="chr22",ldfilter="DprimeLOD"){
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  if(chrom==chrs[1]){data('ldatachr1', envir=environment());ldata<-get("ldatachr1")
  } else if(chrom==chrs[2]){data('ldatachr2', envir=environment());ldata<-get("ldatachr2")
  } else if(chrom==chrs[3]){data('ldatachr3', envir=environment());ldata<-get("ldatachr3")
  } else if(chrom==chrs[4]){data('ldatachr4', envir=environment());ldata<-get("ldatachr4")
  } else if(chrom==chrs[5]){data('ldatachr5', envir=environment());ldata<-get("ldatachr5")
  } else if(chrom==chrs[6]){data('ldatachr6', envir=environment());ldata<-get("ldatachr6")
  } else if(chrom==chrs[7]){data('ldatachr7', envir=environment());ldata<-get("ldatachr7")
  } else if(chrom==chrs[8]){data('ldatachr8', envir=environment());ldata<-get("ldatachr8")
  } else if(chrom==chrs[9]){data('ldatachr9', envir=environment());ldata<-get("ldatachr9")
  } else if(chrom==chrs[10]){data('ldatachr10', envir=environment());ldata<-get("ldatachr10")
  } else if(chrom==chrs[11]){data('ldatachr11', envir=environment());ldata<-get("ldatachr11")
  } else if(chrom==chrs[12]){data('ldatachr12', envir=environment());ldata<-get("ldatachr12")
  } else if(chrom==chrs[13]){data('ldatachr13', envir=environment());ldata<-get("ldatachr13")
  } else if(chrom==chrs[14]){data('ldatachr14', envir=environment());ldata<-get("ldatachr14")
  } else if(chrom==chrs[15]){data('ldatachr15', envir=environment());ldata<-get("ldatachr15")
  } else if(chrom==chrs[16]){data('ldatachr16', envir=environment());ldata<-get("ldatachr16")
  } else if(chrom==chrs[17]){data('ldatachr17', envir=environment());ldata<-get("ldatachr17")
  } else if(chrom==chrs[18]){data('ldatachr18', envir=environment());ldata<-get("ldatachr18")
  } else if(chrom==chrs[19]){data('ldatachr19', envir=environment());ldata<-get("ldatachr19")
  } else if(chrom==chrs[20]){data('ldatachr20', envir=environment());ldata<-get("ldatachr20")
  } else if(chrom==chrs[21]){data('ldatachr21', envir=environment());ldata<-get("ldatachr21")
  } else if(chrom==chrs[22]){data('ldatachr22', envir=environment());ldata<-get("ldatachr22")
  } else if(chrom==chrs[23]){data('ldatachrX', envir=environment());ldata<-get("ldatachrX")
  }
  if(ldfilter=="DprimeLOD"){
    ldata<-ldata[ldata$DPLOD>0,]
  } else {
    ldata<-ldata[ldata$Rsquare>0,]
  }
  return(ldata)
}