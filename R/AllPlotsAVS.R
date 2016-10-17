
###########################################################################
###AVS plots (VSE/EVSE methods)
###########################################################################

#--------------------------------------------------------------------------
avs.plot1<-function(object, what="vse", fname=what, ylab="genomic annotation",
                    xlab="Number of clusters mapping to genomic annotation",
                    breaks="Sturges", maxy=200, pValueCutoff=1e-2, width=8, height=3){
  #checks
  if(class(object)!="AVS"){
    stop("not an 'AVS' object!")
  }
  tnai.checks(name="avs.plot.what",para=what)
  tnai.checks(name="fname",para=fname)
  tnai.checks(name="ylab",para=ylab)
  tnai.checks(name="xlab",para=xlab)
  tnai.checks(name="maxy",para=maxy)
  tnai.checks(name="pValueCutoff",para=pValueCutoff)
  tnai.checks(name="width",para=width)
  tnai.checks(name="height",para=height)
  
  if( !is.null(object@results[[what]]) && object@status[toupper(what)]=="[x]" ){
    for(ix in names(object@results[[what]])){
      resavs<-object@results[[what]][[ix]]
      avsplot1(resavs$mtally,resavs$nulldist,fname=paste(fname,"_",ix,sep=""),
               ylab=ylab, xlab=xlab, breaks=breaks, maxy=maxy, 
               pValueCutoff=pValueCutoff, width=width,height=height)
    }
  } else {
    warning(paste("NOTE:",toupper(what),"results not available!"),call. = FALSE)
  }
}

#--------------------------------------------------------------------------
avs.plot2<-function(object, what="evse", fname=what, width=14, height=2.5,
                    rmargin=1, bxseq=seq(-4,8,2), decreasing=TRUE, ylab="Annotation",
                    xlab="Clusters of risk-associated and linked SNPs", tfs=NULL){
  if(class(object)!="AVS"){
    stop("not an 'AVS' object!")
  }
  tnai.checks(name="avs.plot.what",para=what)
  tnai.checks(name="tfs",para=tfs)
  tnai.checks(name="fname",para=fname)
  tnai.checks(name="width",para=width)
  tnai.checks(name="height",para=height)
  tnai.checks(name="rmargin",para=rmargin)
  tnai.checks(name="bxseq",para=bxseq)
  tnai.checks(name="decreasing",para=decreasing)
  tnai.checks(name="ylab",para=ylab)
  tnai.checks(name="xlab",para=xlab)
  
  if( !is.null(object@results[[what]]) && object@status[toupper(what)]=="[x]" ){
    stats<-object@results$stats[[what]]
    if(any(names(stats)=="escore")){
      stats$score<-stats$escore #just to correct label!
    }
    if(!is.null(tfs)){
      stats<-gettfs(stats,tfs)
    }
    ucounts<-object@results$counts[[what]]
    avsplot2(stats=stats,ucounts=ucounts,fname=fname, height=height, width=width, rmargin=rmargin, 
             bxseq=bxseq, ylab=ylab, xlab=xlab, label=toupper(what),decreasing=decreasing)
  } else {
    warning(paste("NOTE:",toupper(what),"results not available!"),call. = FALSE)
  }
}
gettfs<-function(stats,tfs){
  if(!all(tfs %in% colnames(stats$mtally))){
    stop("All 'tfs' should be available in the AVS object!")
  }
  stats$mtally<-stats$mtally[,tfs]
  stats$nulldist<-stats$nulldist[,tfs]
  stats$score<-stats$score[tfs]
  stats$pvalue<-stats$pvalue[tfs]
  stats
}

###########################################################################
###########################################################################

#-------------------------------------------------------------------------
avsplot1<-function(mtally,nulldist,fname,ylab, xlab, breaks, maxy, 
                   pValueCutoff,width,height){
  nclusters<-length(mtally)
  mtally<-sum(mtally)
  pdf(file=paste(fname,".pdf",sep=""),width=width,height=height)
  xlim=c(0,nclusters)
  ylim=c(0,maxy)
  layout(matrix(c(1,2), 2, 1),heights=c(2,1.2))
  #---plot 1
  par(mar=c(0,4.2,2,1),mgp=c(1.75,0.5,0),tcl=-0.2)
  plot(0,ylim=ylim,xlim=xlim,type="n",axes=FALSE, ylab=paste(ylab,"\n(frequency)",sep=""), xlab="")
  axis(2,yaxp=c(0,ylim[2],4),las=2,cex.axis=0.8,pos=-1.5)
  lines(y=c(0,0),x=c(0,xlim[2]), lty=5,cex=3)
  hist(nulldist,axes=FALSE,main="",ylab="",add=TRUE,col="grey90",breaks=breaks)
  #---plot 2
  par(mar=c(4,4.2,0,1),mgp=c(1.75,0.5,0),tcl=-0.2)
  plot(0,ylim=c(0,1.5),xlim=xlim,type="n",axes=FALSE, ylab="", xlab=xlab)
  boxplot(nulldist,horizontal=TRUE,axes=FALSE,add=TRUE, boxwex=2, outwex=1, range=1.5, 
          pch="|",cex=0.4,lty=1,ldy=0.5, outline=FALSE)
  #axis(1,xaxp=c(0,xlim[2],xlim[2]),cex.axis=0.8, pos=-0.5, padj=-0.25)
  axis(side=1, at=seq(0,xlim[2],2), labels=FALSE, cex.axis=0.65, pos=-0.5, padj=-0.25)
  mtext(seq(0,xlim[2],2), side=1, at=seq(0,xlim[2],2), line = 0.45, cex=0.8)
  #----stat
  nper<-length(nulldist)
  pvalue <- 1-( sum(mtally>nulldist)/nper )
  pcol<-ifelse(pvalue<pValueCutoff,"red","black")
  points(x=mtally,y=1,col=pcol,pch=18,cex=1.5)
  pval<-as.character(ifelse(pvalue<(1/nper),"P < 0.001",signif(pvalue,3)))
  if(pvalue<pValueCutoff)mtext(pval,side=3,at=nclusters-4,line=-0.75,cex=0.8)
  dev.off()
}

#-------------------------------------------------------------------------
avsplot2<-function(stats,ucounts,fname="vseplot",height=2, width=14, rmargin=1, 
                   bxseq=seq(-4,8,2),ylab=ylab, xlab="Clusters of risk-associated and linked SNPs",
                   label="EVSE",decreasing=TRUE){
  #---shortcut to set left margin
  lmargin=0.3
  #---universeCounts
  clustersz<-ucounts[,"markers"]
  annotsz<-ucounts[,"annotation"]
  names(clustersz)<-rownames(ucounts)
  names(annotsz)<-rownames(ucounts)
  if(ncol(ucounts)==3){
    ovlab<-"genes"
  } else {
    ovlab<-"overlaps"
  }
  #---mtally
  mtally<-stats$mtally
  labs<-rownames(mtally)
  clustersz<-clustersz[labs]
  annotsz<-annotsz[labs]
  nulldist<-stats$nulldist
  mtl<-colSums(mtally)
  #----stats
  nper<-nrow(nulldist)
  pvalue<-stats$pvalue
  score<-stats$score
  ci<-stats$ci
  #----sort by score (untie with mtally)
  tp<-colSums(mtally);tp<-tp/sum(tp)
  idx<-sort.list(score+tp,decreasing=!decreasing)
  score<-score[idx]
  mtally<-mtally[,idx,drop=FALSE]
  nulldist<-nulldist[,idx,drop=FALSE]
  mtl<-mtl[idx]
  pvalue<-pvalue[idx]
  #----set labels on 1 or 2 cols
  cnames<-colnames(nulldist)
  ccnames<-strsplit(cnames, ".", fixed=TRUE)
  if(all(unlist(lapply(ccnames,length))==2) && length(ccnames)==length(cnames)){
    cnames<-as.data.frame(ccnames,stringsAsFactors=FALSE)
    ylab<-unlist(strsplit(ylab, ".", fixed=TRUE))
    cctype<-TRUE
  } else {
    cctype<-FALSE
  }
  #---start pdf
  pdf(file=paste(fname,".pdf",sep=""),width=width,height=height)
  nc<-ncol(nulldist)
  ylim<-c(0.35,nc+0.65)
  at.y<-seq(1,nc,1)
  xlim=c(min(bxseq),max(bxseq))
  layout(matrix(c(1,2), 1, 2),widths=c(1+(lmargin/3),3))
  #---plot 1
  par(mai=c(0.4,2+lmargin,1.5,0.1),mgp=c(2,0.5,0), tcl=-0.2)
  plot.new()
  par(usr=c(xlim,ylim))
  abline(v=ci, lmitre=5, col="gray75", lwd=1.0)
  boxplot(nulldist, horizontal=TRUE, axes=FALSE, add=TRUE, boxwex=0.7, range=1.5, pch="|", 
          cex=0.6, lty=1, lwd=0.75, at=at.y, outline=FALSE)
  axis(side=3,at=bxseq[-1], labels=bxseq[-1], cex.axis=0.7, padj=0.5, hadj=0.5, las=1, lwd=1.0)
  tx<-paste("Enrichment\nscore (",label,")",sep="")
  mtext(tx, side=3, line=2, cex=0.75, las=1, adj=0, at=-0.5)
  mtext("Bonferroni", side=3, line=1.2, cex=0.75, las=1, col="gray75", adj=0, at=-0.5)
  #----add stats
  pcol<-ifelse(score>ci,"red","black")
  points(x=score,y=at.y,col=pcol,pch=18,cex=1)
  pval<- round(-log10(pvalue),2)
  mtext(format(pval), side=2, at=at.y, line=0.45, cex=0.75, las=2)
  if(cctype){
    mtext(ylab[1], side=2, line=10, cex=0.75, las=2, adj=0, padj=0, at=nc+1)
    mtext(cnames[1,], side=2, at=at.y, line=10, cex=0.75, las=2, adj=0)
    if(length(ylab)>1)mtext(ylab[2], side=2, line=6, cex=0.75, las=2, adj=0, padj=0, at=nc+1)
    mtext(cnames[2,], side=2, at=at.y, line=6, cex=0.75, las=2, adj=0)
  } else {
    mtext(ylab[1], side=2, line=8, cex=0.75, las=2, adj=0, padj=0, at=nc+1)
    mtext(cnames, side=2, at=at.y, line=8, cex=0.75, las=2, adj=0)
  }
  mtext("P value\n(-log10)", side=2, line=1, cex=0.75, las=2, adj=0.5, padj=0, at=nc+1)
  #---plot 2
  binsp <- ( 85 - length(clustersz) ) * 0.115 * rmargin
  par( mai = c( 0.4, 0.0, 1.5, 0.7 + binsp ) )
  mat.tally<-mtally
  mat.tally[,]<-as.numeric(mtally)
  mat.tally[,score>ci]<-2
  mat.tally[!mtally]<-NA
  nc<-(1/nc)
  #---set cols for 1 or 2 levels
  colgrid<-c("grey75","grey50")
  lv<-levels(as.factor(mat.tally))
  if(length(lv)==1){
    if(lv==2){
      colgrid<-"grey50"
    } else {
      colgrid<-"grey75"
    }
  }
  #---
  sp1<-1/(length(labs)-1)
  sp2<-1/(ncol(mtally)-1)
  image(mat.tally,axes=FALSE,col=colgrid,ylim=c(-nc,1+nc),cex=2)
  mtext(xlab, side=3, line=5, cex=0.8, las=1, adj=0.5)
  mtext(labs, side=3, at=seq(0,1,sp1), line=1.2, cex=0.73, las=2, adj=0)
  mtext(clustersz, side=3, at=seq(0,1,sp1), line=1, cex=0.7, las=2, adj=1)
  mtext("markers", side=3, at=1+sp1, line=0.8, cex=0.7, las=1, adj=0, padj=1, font=3)
  mtext(mtl, side=4, at=seq(0,1,sp2), line=0.5, cex=0.7, las=2, adj=0.5)
  mtext("mapping tally", side=3, at=1+sp1, line=1.2, cex=0.7, las=2, adj=0, font=3)
  if(ovlab=="genes"){
    mtext(annotsz, side=1, at=seq(0,1,sp1), line=0.3, cex=0.7, las=2, adj=1)
    mtext(ovlab, side=1, at=1+sp1, line=-0.5, cex=0.7, las=1, adj=0, padj=1, font=3)
  }
  dev.off()
}
