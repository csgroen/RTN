# Unit tests fot TNI-class methods
test_tni <- function(){
  data(dt4rtn)
  tfs4test<-dt4rtn$tfs[c("PTTG1","FOXM1")]
  rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=tfs4test)
  rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
  #tni.permutation
  rtni<-tni.permutation(rtni,nPermutations=10)
  res<-tni.get(rtni,what="refnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.bootstrap
  rtni<-tni.bootstrap(rtni,nBootstraps=10)
  res<-tni.get(rtni,what="refnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.dpi.filter
  rtni<-tni.dpi.filter(rtni)
  res<-tni.get(rtni,what="tnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.conditional
  annot <- res<-tni.get(rtni,what="annotation")
  idx <- annot$SYMBOL%in%c("FGF2","ERBB2")
  mod4test<-annot$PROBEID[idx]
  rtni<-tni.conditional(rtni, modulators=mod4test, minRegulonSize=1)
  res<-tni.get(rtni,what="cdt")
  checkTrue(is.list(res))
  #tni.graph
  # res<-tni.graph(rtni, gtype="rmap", tfs=tfs4test)
  # checkTrue(is.igraph(res))
}
# Unit tests fot TNA-class methods
test_tna <- function(){
  data(dt4rtn)
  tfs4test<-dt4rtn$tfs[c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")]
  rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=tfs4test)
  rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
  rtni<-tni.permutation(rtni,nPermutations=10)
  rtni<-tni.bootstrap(rtni,nBootstraps=10)
  rtni<-tni.dpi.filter(rtni)
  rtna<-tni2tna.preprocess(rtni, phenotype=dt4rtn$pheno, hits=dt4rtn$hits, phenoIDs=dt4rtn$phenoIDs)
  #tna.mra
  rtna <- tna.mra(rtna, minRegulonSize=1)
  res<-tna.get(rtna,what="mra")
  checkTrue(is.data.frame(res) && ncol(res)==8)
  #tna.overlap
  rtna <- tna.overlap(rtna, minRegulonSize=1)
  res<-tna.get(rtna,what="overlap")
  checkTrue(is.data.frame(res) && ncol(res)==9)
  #tna.gsea1
  rtna <- tna.gsea1(rtna, nPermutations=10, minRegulonSize=1, stepFilter=FALSE)
  res <- tna.get(rtna,what="gsea1")
  checkTrue(is.data.frame(res) && ncol(res)==5)
  #tna.gsea2
  rtna <- tna.gsea2(rtna, nPermutations=10, minRegulonSize=1, stepFilter=FALSE)
  res<-tna.get(rtna,what="gsea2")
  checkTrue(is.list(res) && length(res)==3)
  #tna.synergy
  rtna <- tna.synergy(rtna,nPermutations=10, minRegulonSize=1, stepFilter=FALSE)
  res<-tna.get(rtna,what="synergy")
  checkTrue(is.data.frame(res) && ncol(res)==9)
  #tna.shadow
  rtna <- tna.shadow(rtna,nPermutations=10, minRegulonSize=1, stepFilter=FALSE)
  res<-tna.get(rtna,what="shadow")
  checkTrue(is.data.frame(res) && ncol(res)==9)
  #tna.graph
  # res<-tna.graph(rtna, tnet="ref", gtype="amap", tfs=names(tfs4test))
  # checkTrue(is.igraph(res))
}
