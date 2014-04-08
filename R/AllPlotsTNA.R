################################################################################
##########################          PLOTS           ############################
################################################################################

##------------------------------------------------------------------------------
##Plot enrichment analysis from TNA objects.
tna.plot.gsea1<-function(object, regulon.order="size", ntop=NULL, tfs=NULL, filepath=".", 
                        ylimPanels=c(0.0,3.5,0.0,0.8), heightPanels=c(1,1,3), width=6, height=5, 
                        ylabPanels=c("Phenotype","Regulon","Enrichment score"), 
                        xlab="Position in the ranked list of genes", labPheno="tna_test",
                        alpha=0.5, sparsity=10, splitcor=FALSE, autoformat=TRUE, ...) {
  #checks
  if(class(object)!="TNA" || object@status$analysis["GSEA1"]!="[x]"){
    cat("-invalid 'GSEA1' status! \n")
    stop("NOTE: gsea plot requires results from 'tna.gsea1' analysis!")
  }
  regulon.order<-tnai.checks(name="regulon.order",regulon.order)
  tnai.checks(name="ntop",ntop)
  tnai.checks(name="tfs",tfs)
  tnai.checks(name="filepath",filepath)
  tnai.checks(name="ylimPanels",ylimPanels)
  tnai.checks(name="heightPanels",heightPanels)
  tnai.checks(name="width",width)
  tnai.checks(name="height",height)
  tnai.checks(name="ylabPanels",ylabPanels)
  tnai.checks(name="xlab",xlab)
  tnai.checks(name="labPheno",labPheno)
  tnai.checks(name="alpha",alpha)
  tnai.checks(name="autoformat",autoformat)
  tnai.checks(name="splitcor",splitcor)
  if(!is.null(tfs)){
    resgsea<-tna.get(object, what="gsea1", reportNames=TRUE)
    idx<-(rownames(resgsea)%in%tfs+resgsea$Regulon%in%tfs)>0
    if(all(!idx)){
      stop("one or more input 'tfs' not found in the 'gsea1' results!")
    }
    resgsea<-resgsea[idx,]
  } else {
    resgsea<-tna.get(object, what="gsea1", ntop=ntop, reportNames=TRUE)
  }
  if(!is.null(resgsea) && nrow(resgsea)>0){
    if(regulon.order!='none'){
      decreasing<-ifelse(regulon.order=='Observed.Score',TRUE,FALSE)
      resgsea<-resgsea[sort.list(resgsea[,regulon.order],decreasing=decreasing),]
    }
    gs.names<-rownames(resgsea)
    gs.labels<-resgsea$Regulon
  } else {
    stop("gsea1 slot is empty or null!")
  }
  ##-----get gene sets used in the current gsea1 analysis
  if(object@para$gsea1$tnet=="cdt"){
    rgcs<-object@listOfModulators
    if(ylabPanels[2]=="Regulon")ylabPanels[2]<-"Modulators"
  } else if(object@para$gsea1$tnet=="ref"){
    rgcs<-tna.get(object,what="refregulons.and.mode")
  } else {
    rgcs<-tna.get(object,what="regulons.and.mode")
  }
  ##-----get ordered phenotype
  phenotype<-object@phenotype
  if(object@para$gsea1$orderAbsValue)phenotype<-abs(phenotype)
  phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
  ##----get stat resolution
  pvresolu<-signif(1/(object@para$gsea1$nPermutations+1), digits=5)
  pvcutoff<-paste("< ",as.character(format(pvresolu,scientific=TRUE,digits=2)),collapse="",sep="")
  ##-----get merged data
  tests<-get.merged.data1(gs.names,gs.labels,phenotype,rgcs,resgsea,object@para$gsea1$exponent)
  ##-----fix pvalue report for the resolution
  idx<-tests$pv==pvresolu
  tests$adjpv[idx]<-pvcutoff
  ##-----check format
  if(autoformat)ylimPanels<-check.format1(tests)
  ##-----make plot
  if(!splitcor){
    make.plot1(tests,filepath,labPheno,heightPanels,ylimPanels,ylabPanels,xlab,width, height,alpha,sparsity, ...=...)
  } else {
    make.plot1.1(tests,filepath,labPheno,heightPanels,ylimPanels,ylabPanels,xlab,width, height,alpha,sparsity, ...=...)
  }
}
##------------------------------------------------------------------------------
##Plot 2-tailed enrichment analysis from TNA objects.
tna.plot.gsea2<-function(object, regulon.order="size", ntop=NULL, tfs=NULL, filepath=".", 
                        ylimPanels=c(-1.0,3.0,-0.5,0.5), heightPanels=c(1.5,0.7,5.0), width=5, height=4, 
                        ylabPanels=c("Phenotype","Regulon","Enrichment score"), 
                        xlab="Position in the ranked list of genes", labPheno="tna_test",
                        alpha=1.0, sparsity=10, autoformat=TRUE, ...) {
  #checks
  if(class(object)!="TNA" || object@status$analysis["GSEA2"]!="[x]"){
    cat("-invalid 'GSEA2' status! \n")
    stop("NOTE: gsea plot requires results from 'tna.gsea2' analysis!")
  }
  regulon.order<-tnai.checks(name="regulon.order",regulon.order)
  tnai.checks(name="ntop",ntop)
  tnai.checks(name="tfs",tfs)
  tnai.checks(name="filepath",filepath)
  tnai.checks(name="ylimPanels",ylimPanels)
  tnai.checks(name="heightPanels",heightPanels)
  tnai.checks(name="width",width)
  tnai.checks(name="height",height)
  tnai.checks(name="ylabPanels",ylabPanels)
  tnai.checks(name="xlab",xlab)
  tnai.checks(name="labPheno",labPheno)
  tnai.checks(name="alpha",alpha)
  tnai.checks(name="autoformat",autoformat)
  if(!is.null(tfs)){
    resgsea<-tna.get(object, what="gsea2", ntop=-1, reportNames=TRUE)
    idx<-(rownames(resgsea$differential)%in%tfs+resgsea$differential$Regulon%in%tfs)>0
    if(all(!idx)){
      stop("one or more input 'tfs' not found in the 'gsea2' results!")
    }
    resgsea$differential<-resgsea$differential[idx,]
    resgsea$positive<-resgsea$positive[idx,]
    resgsea$negative<-resgsea$negative[idx,]
  } else {
    resgsea<-tna.get(object, what="gsea2", ntop=ntop, reportNames=TRUE)
  }
  if(!is.null(resgsea) && nrow(resgsea$differential)>0){
    if(regulon.order!='none'){
      decreasing<-ifelse(regulon.order=='Observed.Score',TRUE,FALSE)
      idx<-sort.list(resgsea$differential[,regulon.order],decreasing=decreasing)
      resgsea$differential<-resgsea$differential[idx,]
      resgsea$positive<-resgsea$positive[idx,]
      resgsea$negative<-resgsea$negative[idx,]
    }
    gs.names<-rownames(resgsea$differential)
    gs.labels<-resgsea$differential$Regulon
  } else {
    stop("gsea2 is empty or null!")
  }
  ##-----get gene sets used in the current gsea analysis
  if(object@para$gsea2$tnet=="cdt"){
    rgcs<-object@listOfModulators
    if(ylabPanels[2]=="Regulon")ylabPanels[2]<-"Modulators"
  } else if(object@para$gsea2$tnet=="ref"){
    rgcs<-tna.get(object,what="refregulons.and.mode")
  } else {
    rgcs<-tna.get(object,what="regulons.and.mode")
  }
  ##-----get ordered phenotype
  phenotype<-object@phenotype
  phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
  ##----get stat resolution
  pvresolu<-signif(1/(object@para$gsea2$nPermutations+1), digits=5)
  pvcutoff<-paste("< ",as.character(format(pvresolu,scientific=TRUE,digits=2)),collapse="",sep="")
  ##----plot
  for(i in 1:length(gs.names)){
    ##-----get merged data
    tests<-get.merged.data2(gs.names[i],gs.labels[i],phenotype,rgcs,resgsea$differential[i,,drop=FALSE],
                            object@para$gsea2$exponent)
    tests$pv[["up"]]<-resgsea$positive[i,"Pvalue"]
    tests$adjpv[["up"]]<-resgsea$positive[i,"Adjusted.Pvalue"]
    tests$pv[["down"]]<-resgsea$negative[i,"Pvalue"]
    tests$adjpv[["down"]]<-resgsea$negative[i,"Adjusted.Pvalue"]
    tests$adjpv[]<-paste("= ",as.character(format(tests$adjpv,scientific=TRUE,digits=2)),sep="")
    ##-----fix pvalue report for the resolution
    idx<-tests$pv==pvresolu
    tests$adjpv[idx]<-pvcutoff
    ##-----check format
    if(autoformat)ylimPanels<-check.format2(tests)
    ##-----make plot
    make.plot2(tests,filepath,labPheno,heightPanels,ylimPanels,ylabPanels,xlab,width,height,alpha,sparsity, ...=...)
  }
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
get.merged.data1<-function(gs.names,gs.labels,phenotype,rgcs,resgsea,exponent){
  res<-list()
  for(gs.name in gs.names){
    test<-gseaScores4RTN(geneList=phenotype, geneSet=names(rgcs[[gs.name]]), 
                         exponent=exponent, mode="graph",geneSetMode=rgcs[[gs.name]])
    res$enrichmentScores[[gs.name]]<-test$enrichmentScore
    res$runningScores[[gs.name]]<-test$runningScore
    res$positions[[gs.name]]<-test$positions
    res$pvals[[gs.name]]<-resgsea[gs.name,][["Pvalue"]]
    res$adjpvals[[gs.name]]<-resgsea[gs.name,][["Adjusted.Pvalue"]]
  }
  tests<-list()
  tests[["enrichmentScores"]]<-res$enrichmentScores 
  tests[["runningScores"]]<-as.data.frame(res$runningScores,stringsAsFactors=FALSE)
  tests[["positions"]]<-as.data.frame(res$positions,stringsAsFactors=FALSE)
  tests[["pv"]]<-res$pvals
  tests[["adjpv"]]<-paste("= ",as.character(format(res$adjpvals,scientific=TRUE,digits=2)),sep="")
  tests[["geneList"]]<-phenotype
  tests[["labels"]]<-gs.labels
  tests
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea2
get.merged.data2<-function(gs.name,gs.label,phenotype,rgcs,resgsea,exponent){
  res<-list()
  gs<-rgcs[[gs.name]]
  test<-gseaScores4RTN(geneList=phenotype, geneSet=names(gs), 
                       exponent=exponent, mode="graph",geneSetMode=gs)
  res$positions[[gs.name]]<-test$positions
  res$pvals[[gs.name]]<-resgsea[gs.name,][["Pvalue"]]
  res$adjpvals[[gs.name]]<-resgsea[gs.name,][["Adjusted.Pvalue"]]
  testup<-gseaScores4RTN(geneList=phenotype, geneSet=names(gs[gs>0]), 
                         exponent=exponent, mode="graph",geneSetMode=gs[gs>0])
  testdown<-gseaScores4RTN(geneList=phenotype, geneSet=names(gs[gs<0]), 
                           exponent=exponent, mode="graph",geneSetMode=gs[gs<0])
  res$testup$enrichmentScores[[gs.name]]<-testup$enrichmentScore
  res$testup$runningScores[[gs.name]]<-testup$runningScore
  res$testdown$enrichmentScores[[gs.name]]<-testdown$enrichmentScore
  res$testdown$runningScores[[gs.name]]<-testdown$runningScore
  tests<-list()
  tests$testup[["enrichmentScores"]]<-res$testup$enrichmentScores
  tests$testup[["runningScores"]]<-as.data.frame(res$testup$runningScores,stringsAsFactors=FALSE)
  tests$testdown[["enrichmentScores"]]<-res$testdown$enrichmentScores
  tests$testdown[["runningScores"]]<-as.data.frame(res$testdown$runningScores,stringsAsFactors=FALSE)
  tests[["positions"]]<-as.data.frame(res$positions,stringsAsFactors=FALSE)  
  tests[["geneList"]]<-phenotype
  tests[["label"]]<-gs.label
  tests$pv[["pv"]]<-res$pvals
  tests$adjpv[["pv"]]<-res$adjpvals
  tests
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
check.format1<-function(tests){
  ylimPanels<-rep(0,4)
  tp<-c(min(min(tests$geneList)),max(tests$geneList))
  tpp<-as.integer(tp)
  if(tp[1]<tpp[1])tpp[1]=tpp[1]-1
  if(tp[2]>tpp[2])tpp[2]=tpp[2]+1
  ylimPanels[1:2]<-tpp
  tp<-sapply(1:length(tests$labels),function(i){
    c(min(tests$runningScores[,i]),max(tests$runningScores[,i]))
  })
  tp<-c(min(tp),max(tp))
  if(min(tests$geneList)>=0 && tp[1]>=0)tp[1]=0
  tpp<-round(tp,digits=1)
  if(tp[1]<tpp[1])tpp[1]=tpp[1]-0.1
  if(tp[2]>tpp[2])tpp[2]=tpp[2]+0.1
  ylimPanels[3:4]<-tpp
  ylimPanels
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea2
check.format2<-function(tests){
  ylimPanels<-rep(0,4)
  tp<-c(min(min(tests$geneList)),max(tests$geneList))
  tpp<-as.integer(tp)
  if(tp[1]<tpp[1])tpp[1]=tpp[1]-1
  if(tp[2]>tpp[2])tpp[2]=tpp[2]+1
  ylimPanels[1:2]<-tpp
  tp<-sapply(1:length(tests$label),function(i){
    tp1<-min(c(tests$testup$runningScores[,i],tests$testdown$runningScores[,i]))
    tp2<-max(c(tests$testup$runningScores[,i],tests$testdown$runningScores[,i]))
    c(tp1,tp2)
  })
  tp<-c(min(tp),max(tp))
  if(min(tests$geneList)>=0 && tp[1]>=0)tp[1]=0
  tpp<-round(tp,digits=1)
  if(tp[1]<tpp[1])tpp[1]=tpp[1]-0.1
  if(tp[2]>tpp[2])tpp[2]=tpp[2]+0.1
  tpp[2]<-ifelse(tpp[2]>=0.3,tpp[2],0.3)
  ylimPanels[3:4]<-tpp
  ylimPanels
}


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
make.plot1 <- function(tests, filepath, labPheno, heightPanels, ylimPanels, ylabPanels,
                       xlab, width, height, alpha, sparsity, ...) {
  pdf(file=file.path(filepath, paste(labPheno,".pdf", sep="")), width=width, height=height)
  gsea.plot1(tests$runningScore, tests$enrichmentScore, tests$positions, tests$adjpv, tests$geneList, 
             tests$labels, heightPanels, ylimPanels, ylabPanels, xlab, labPheno, alpha, sparsity, ...=... )
  dev.off()
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
gsea.plot1 <- function(runningScore, enrichmentScore, positions, adjpv, 
                       geneList, labels, heightPanels, ylimPanels, ylabPanels, xlab, 
                       labPheno, alpha, sparsity, ...) {
  #-------------------------------------------------
  #set text size levels
  cexlev=c(1.3,1.2,1.1,0.9,0.8)
  #set colors
  rsc<-positions
  ng<-ncol(rsc)
  get.alpha<-function (colour,alpha=1.0) {
    col <- col2rgb(colour, TRUE)/255
    alpha <- rep(alpha, length.out = length(colour))
    rgb(col[1, ], col[2, ], col[3, ], alpha)
  }
  rsc.colors<-get.alpha(palette(), alpha)
  if(ng>length(rsc.colors))rsc.colors<-get.alpha(colorRampPalette(rsc.colors)(ng),alpha)
  #-------------------------------------------------
  #set hits
  for(i in 1:ng){
    idx<-rsc[,i]!=0
    rsc[idx,i]<-i
  }      	
  rsc<-as.matrix(rsc)
  rsc.vec<-cbind(rep(1:nrow(rsc),ng),as.vector(rsc))
  rsc.vec<-rsc.vec[rsc.vec[,2]!=0,]
  #-------------------------------------------------
  # layout
  np<-sum(heightPanels>0)
  ht<-heightPanels
  ht<-ht[ht>0]
  layout(matrix(1:np, np, 1, byrow=TRUE), heights=ht)
  # plot1
  par(family="Times")
  if(heightPanels[1]>0){
    par(mar=c(0.5, 6.5, 1.5, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    plot(x=c(1,max(rsc.vec[,1])),y=c(min(geneList),max(geneList)), type="n", 
         axes= FALSE,xlab="", ylab=ylabPanels[1], cex.lab=cexlev[1], ylim=ylimPanels[1:2], ...=...)
    if(min(geneList)<0)abline(h=0,lwd=0.6)
    sq<-c(1:length(geneList))%%sparsity;sq[sq>1]<-0
    sq<-as.logical(sq)
    lines(x=c(1:length(geneList))[sq],y=geneList[sq],col="grey75",lwd=1.2)
    nn=ifelse(min(geneList)<0,4,3)
    pp<-pretty(c(geneList,ylimPanels[1:2]),n=nn)
    axis(2,line=0, cex.axis=cexlev[2], las=2, at=pp, labels=pp,...=...)
    #box(lwd=0.8)
    if(!is.null(labPheno)){
      legend("topright", legend=labPheno, col="grey75", pch="---", 
             bty="n",cex=cexlev[3], pt.cex=1.0, ...=...)   	
    }
  }
  #-------------------------------------------------
  # plot2
  if(heightPanels[2]>0){
    par(mar=c(0.0, 6.5, 0.0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    plot(x=c(1,max(rsc.vec[,1])),y=c(0,max(rsc.vec[,2])), type="n", axes=FALSE, xlab="",
         ylab=ylabPanels[2],cex.lab=cexlev[1], ...=...)
    for(i in 1:ng){
      idx<-rsc.vec[,2]==i
      xx<-rsc.vec[idx,1]
      yy<-rsc.vec[idx,2]
      segments(xx,yy-0.9,xx, yy-0.1, col=rsc.colors[i],lwd=0.2)
    }
    axis(2,las=2, at=c(1:ng)-0.5,labels=labels, line=0, cex.axis=cexlev[5], ...=...)
    #box(lwd=0.8)
  }
  #-------------------------------------------------
  # plot3
  if(heightPanels[3]>0){
    rsc.colors<-get.alpha(rsc.colors,1.0)
    par(mar=c(5, 6.5, 0.0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    cc<-as.matrix(runningScore)
    plot(x=c(1,nrow(cc)),y=c(min(cc),max(cc)), type="n", axes=FALSE, xlab="",
         ylim=ylimPanels[c(3,4)],ylab=ylabPanels[3], cex.lab=cexlev[1], ...=...)
    par(mgp=c(3.0,0.5,0))
    title(xlab=xlab, cex.lab=cexlev[1], ...=...)
    if(min(cc)<0)abline(h=0,lwd=0.6)
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScore[[i]])
      segments(xx,0, xx, yy[xx], col=rsc.colors[i],lwd=1, lty=3)
    }
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      sq<-c(1:length(xx))%%sparsity;sq[sq>1]<-0
      sq<-as.logical(sq)
      lines(x=xx[sq],y=yy[sq],col=rsc.colors[i],lwd=0.7)
    }
    axis(1,cex.axis=cexlev[2], ...=...)
    axis(2,las=2,cex.axis=cexlev[2], ...=...)
    #box(lwd=0.8)
    labels<-paste(labels," (adj.p ",format(adjpv,scientific=TRUE,digits=2),")",sep="")
    #labels=sub("=","<",labels)
    legend("topright", legend=labels, col=rsc.colors, pch="---", bty="n",cex=cexlev[4], 
           pt.cex=1.0, title=ylabPanels[2], title.adj = 0, ...=...)
  }
}


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
make.plot1.1 <- function(tests, filepath, labPheno, heightPanels, ylimPanels, ylabPanels,
                       xlab, width, height, alpha, sparsity, ...) {
  pdf(file=file.path(filepath, paste(labPheno,".pdf", sep="")), width=width, height=height)
  gsea.plot1.1(tests$runningScore, tests$enrichmentScore, tests$positions, tests$adjpv, tests$geneList, 
             tests$labels, heightPanels, ylimPanels, ylabPanels, xlab, labPheno, alpha, sparsity, ...=... )
  dev.off()
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea1
gsea.plot1.1 <- function(runningScore, enrichmentScore, positions, adjpv, geneList, 
                             labels, heightPanels, ylimPanels, ylabPanels, xlab, labPheno, 
                             alpha, sparsity, ...) {
  #-------------------------------------------------
  positions<-as.matrix(positions)
  positions[positions==1]=2
  positions[positions==-1]=1
  #set text size levels
  cexlev=c(1.3,1.2,1.1,0.9,0.8)
  #set colors
  ng<-ncol(positions)
  get.alpha<-function (colour,alpha=1.0) {
    col <- col2rgb(colour, TRUE)/255
    alpha <- rep(alpha, length.out = length(colour))
    rgb(col[1, ], col[2, ], col[3, ], alpha)
  }
  rsc.colors1<-get.alpha(palette(), alpha)
  if(ng>length(rsc.colors1))rsc.colors1<-get.alpha(colorRampPalette(rsc.colors1)(ng),alpha)
  rsc.colors2<-get.alpha(c("#96D1FF","#FF8E91"), alpha)
  #-------------------------------------------------
  #set hits
  rsc1<-rsc2<-positions
  for(i in 1:ng){
    idx<-rsc1[,i]!=0
    rsc1[idx,i]<-i
  }
  rsc.vec<-cbind(rep(1:nrow(rsc1),ng),as.vector(rsc1),as.vector(rsc2))
  rsc.vec<-rsc.vec[rsc.vec[,2]!=0,]  
  #-------------------------------------------------
  #-------------------------------------------------
  # layout
  np<-sum(heightPanels>0)
  ht<-heightPanels
  ht<-ht[ht>0]
  layout(matrix(1:np, np, 1, byrow=TRUE), heights=ht)
  # plot1
  par(family="Times")
  if(heightPanels[1]>0){
    par(mar=c(0.5, 6.5, 1.5, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2,family="Times")
    plot(x=c(1,max(rsc.vec[,1])),y=c(min(geneList),max(geneList)), type="n", 
         axes= FALSE,xlab="", ylab=ylabPanels[1], cex.lab=cexlev[1], ylim=ylimPanels[1:2], ...=...)
    if(min(geneList)<0)abline(h=0,lwd=0.6)    #segments(0, 0, length(geneList), 0,col="grey70")
    sq<-c(1:length(geneList))%%sparsity;sq[sq>1]<-0
    sq<-as.logical(sq)
    lines(x=c(1:length(geneList))[sq],y=geneList[sq],col="grey75",lwd=1.2)
    nn=ifelse(min(geneList)<0,4,3)
    pp<-pretty(c(geneList,ylimPanels[1:2]),n=nn)
    axis(2,line=0, cex.axis=cexlev[2], las=2, at=pp, labels=pp,...=...)
    #box()
    if(!is.null(labPheno)){
      legend("topright", legend=labPheno, col="grey75", pch="---", 
             bty="n",cex=cexlev[3], pt.cex=1.0, ...=...)     
    }
#     legend("bottom", legend=c("negative interaction","positive interaction"), col=rsc.colors2, pch="---", bty="n",cex=cexlev[4],horiz=TRUE, 
#            pt.cex=1.4, title=NULL, title.adj = 0,inset=-0.05, ...=...)
  }
  #-------------------------------------------------
  # plot2
  if(heightPanels[2]>0){
    par(mar=c(0.0, 6.5, 1.0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    plot(x=c(1,max(rsc.vec[,1])),y=c(0,max(rsc.vec[,2])), type="n", axes=FALSE, xlab="",
         ylab=ylabPanels[2],cex.lab=cexlev[1], ...=...)
    for(i in 1:ng){
      idx<-rsc.vec[,2]==i
      xx<-rsc.vec[idx,1]
      yy<-rsc.vec[idx,2]
      cc<-rsc.vec[idx,3]
      segments(xx,yy-0.9,xx, yy-0.1, col=rsc.colors2[cc],lwd=0.2)
    }
    axis(2,las=2, at=c(1:ng)-0.5,labels=labels, line=0, cex.axis=cexlev[5], ...=...)
    par(xpd=TRUE,mar=c(0.0, 6.5, 0.0, 1.0))
    legend("top", legend=c("negative interaction","positive interaction"), col=rsc.colors2, pch="---", bty="n",cex=cexlev[4],horiz=TRUE, 
           pt.cex=1.4, title=NULL, title.adj = 0,inset=-0.18, ...=...)
    par(xpd=FALSE,mar=c(0.0, 6.5, 1.0, 1.0))
  }
  #-------------------------------------------------
  # plot3
  if(heightPanels[3]>0){
    rsc.colors1<-get.alpha(rsc.colors1,1.0)
    par(mar=c(5, 6.5, 0.0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    cc<-as.matrix(runningScore)
    plot(x=c(1,nrow(cc)),y=c(min(cc),max(cc)), type="n", axes=FALSE, xlab="",
         ylim=ylimPanels[c(3,4)],ylab=ylabPanels[3], cex.lab=cexlev[1], ...=...)
    par(mgp=c(3.0,0.5,0))
    title(xlab=xlab, cex.lab=cexlev[1], ...=...)
    if(min(cc)<0)abline(h=0,lwd=0.6)   #segments(0, 0, nrow(cc), 0,col="grey70")
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScore[[i]])
      segments(xx,0, xx, yy[xx], col=rsc.colors1[i],lwd=1, lty=3)
    }
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      sq<-c(1:length(xx))%%sparsity;sq[sq>1]<-0
      sq<-as.logical(sq)
      lines(x=xx[sq],y=yy[sq],col=rsc.colors1[i],lwd=0.7)
    }
    axis(1,cex.axis=cexlev[2], ...=...)
    axis(2,las=2,cex.axis=cexlev[2], ...=...)
    #box()
    labels<-paste(labels," (adj.p ",format(adjpv,scientific=TRUE,digits=2),")",sep="")
    #labels=sub("=","<",labels)
    legend("topright", legend=labels, col=rsc.colors1, pch="---", bty="n",cex=cexlev[4], 
           pt.cex=1.0, title=ylabPanels[2], title.adj = 0, ...=...)
  }
}


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea2
make.plot2 <- function(tests, filepath, labPheno, heightPanels, ylimPanels, ylabPanels,
                       xlab, width, height, alpha, sparsity, ...) {
  pdf(file=file.path(filepath, paste(labPheno,"_",tests$label,".pdf", sep="")), width=width, height=height)
  gsea.plot2(tests$testup$runningScore, tests$testup$enrichmentScore, 
             tests$testdown$runningScore, tests$testdown$enrichmentScore, 
             tests$positions, tests$adjpv, tests$geneList, tests$label, heightPanels, 
             ylimPanels, ylabPanels, xlab, labPheno, alpha, sparsity, ...=... )
  dev.off()
}
#-------------------------------------------------------------------------------------
#--subfunction for tna.plot.gsea2
gsea.plot2 <- function(runningScoreUp, enrichmentScoreUp, runningScoreDown, enrichmentScoreDown,
                       positions, adjpv, geneList, label, heightPanels, ylimPanels, ylabPanels, 
                       xlab, labPheno, alpha, sparsity, ...) {
  #-------------------------------------------------
  positions<-as.matrix(positions)
  positions[positions==1]=2
  positions[positions==-1]=1
  #set text size levels
  cexlev=c(1.3,1.2,1.1,0.9,0.8)
  #set colors
  ng<-ncol(positions)
  get.alpha<-function (colour,alpha=1.0) {
    col <- col2rgb(colour, TRUE)/255
    alpha <- rep(alpha, length.out = length(colour))
    rgb(col[1, ], col[2, ], col[3, ], alpha)
  }
  rsc.colors<-get.alpha(c("#96D1FF","#FF8E91"), alpha)
  #-------------------------------------------------
  #set hits
  rsc1<-rsc2<-positions
  for(i in 1:ng){
    idx<-rsc1[,i]!=0
    rsc1[idx,i]<-i
  }
  rsc.vec<-cbind(rep(1:nrow(rsc1),ng),as.vector(rsc1),as.vector(rsc2))
  rsc.vec<-rsc.vec[rsc.vec[,2]!=0,]  
  #-------------------------------------------------
  # layout
  np<-sum(heightPanels>0)
  ht<-heightPanels
  ht<-ht[ht>0]
  layout(matrix(1:np, np, 1, byrow=TRUE), heights=ht)
  # plot1
  par(family="Times")
  if(heightPanels[1]>0){
    xlim=c(0,length(geneList))
    par(mar=c(0.0, 6.5, 1.5, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2,family="Times")
    plot(x=c(1,max(rsc.vec[,1])),y=c(min(geneList),max(geneList)), type="n", 
         axes= FALSE,xlab="", ylab=ylabPanels[1], cex.lab=cexlev[1], ylim=ylimPanels[1:2],xlim=xlim, ...=...)
    if(min(geneList)<0)abline(h=0,lwd=0.6)    #segments(0, 0, length(geneList), 0,col="grey70")
    sq<-c(1:length(geneList))%%sparsity;sq[sq>1]<-0
    sq<-as.logical(sq)
    lines(x=c(1:length(geneList))[sq],y=geneList[sq],col="grey75",lwd=1.2)
    nn=ifelse(min(geneList)<0,4,3)
    pp<-pretty(c(geneList,ylimPanels[1:2]),n=nn)
    axis(2,line=0, cex.axis=cexlev[2], las=2, at=pp, labels=pp,...=...)
    if(!is.null(labPheno)){
      legend("topright", legend=labPheno, col="grey75", pch="---", 
             bty="n",cex=cexlev[3], pt.cex=1.0, ...=...)     
    }
  }
  #-------------------------------------------------
  # plot2
  if(heightPanels[2]>0){
    par(mar=c(0.0, 6.5, 0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    plot(x=c(1,max(rsc.vec[,1])),y=c(0,max(rsc.vec[,2])+1), type="n", axes=FALSE, xlab="",
         ylab="",cex.lab=cexlev[1], ...=...)
    for(i in 1:ng){
      idx<-rsc.vec[,2]==i
      xx<-rsc.vec[idx,1]
      yy<-rsc.vec[idx,2]
      cc<-rsc.vec[idx,3]
      segments(xx,yy-0.9,xx, yy-0.1, col=rsc.colors[cc],lwd=0.2)
    }
    axis(2,las=2, at=c(1:ng)-0.5,labels=ylabPanels[2], line=0, cex.axis=cexlev[3], ...=...)
    legend("top", legend=c("negative interaction","positive interaction"), col=rsc.colors, 
           pch="---", bty="n",cex=cexlev[4],horiz=TRUE, pt.cex=1.4, title=NULL, title.adj = 0,
           inset=0, ...=...)
  }
  #-------------------------------------------------
  # plot3
  if(heightPanels[3]>0){
    rsc.colors<-get.alpha(rsc.colors,1.0)
    par(mar=c(5, 6.5, 0.0, 1.0),mgp=c(4.5,0.5,0),tcl=-0.2)
    cc<-as.matrix(runningScoreDown)
    plot(x=c(1,nrow(cc)),y=c(min(cc),max(cc)), type="n", axes=FALSE, xlab="",
         ylim=ylimPanels[c(3,4)],ylab=ylabPanels[3], cex.lab=cexlev[1], ...=...)
    par(mgp=c(3.0,0.5,0))
    title(xlab=xlab, cex.lab=cexlev[1], ...=...)
    if(min(cc)<0)abline(h=0,lwd=0.6)
    #---
    cc<-as.matrix(runningScoreDown)
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScoreDown[[i]])
      segments(xx,0, xx, yy[xx], col=rsc.colors[1],lwd=1.2, lty=3)
    }
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      sq<-c(1:length(xx))%%sparsity;sq[sq>1]<-0
      sq<-as.logical(sq)
      lines(x=xx[sq],y=yy[sq],col=rsc.colors[1],lwd=1.0)
    }
    #---
    cc<-as.matrix(runningScoreUp)
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScoreUp[[i]])
      segments(xx,0, xx, yy[xx], col=rsc.colors[2],lwd=1.2, lty=3)
    }
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      sq<-c(1:length(xx))%%sparsity;sq[sq>1]<-0
      sq<-as.logical(sq)
      lines(x=xx[sq],y=yy[sq],col=rsc.colors[2],lwd=1.0)
    }
    cc<-as.matrix(runningScoreDown)
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScoreDown[[i]])
      points(xx, yy[xx], bg=rsc.colors[1],col=rsc.colors[1], lwd=1, cex=1, pch=21)
    }
    cc<-as.matrix(runningScoreUp)
    for(i in 1:ng){
      yy<-cc[,i]
      xx<-c(1:nrow(cc))
      xx<-which(yy==enrichmentScoreUp[[i]])
      points(xx, yy[xx], bg=rsc.colors[2],col=rsc.colors[2], lwd=1, cex=1, pch=21)
    }
    #---
    axis(1,cex.axis=cexlev[2], ...=...)
    axis(2,las=2,cex.axis=cexlev[2], ...=...)
    adjpv<-c(adjpv["up"],adjpv["down"],adjpv["pv"])
    lbstat<-paste(c("Positive","Negative","Differential")," (adj.p ",adjpv,")",sep="")
    legend("topright", legend=lbstat, col=c(rsc.colors[2],rsc.colors[1],NA), pch=20, bty="n",cex=cexlev[4], 
           pt.cex=1.2, title=ylabPanels[2], title.adj = 0,  ...=...)
    legend("bottomleft", legend=label, col=NA, pch=NA, bty="n",cex=cexlev[1]*1.3, pt.cex=1.2, title=NULL,  ...=...)
  }
}


