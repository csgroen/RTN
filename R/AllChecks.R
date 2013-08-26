##This function is used for argument checking
tnai.checks <- function(name, para) {
  if(name=="alpha") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'alpha' should be an integer or numeric value <=1 and >=0!\n")
  }
  else if(name=="tna.what"){
    opts<-c("tnet","tfs","pheno","regulons","refregulons","para","mra","gsea1","gsea2","overlap","synergy","shadow","summary","status",
            "regulons.and.pheno","refregulons.and.pheno","regulons.and.mode","refregulons.and.mode")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ) )
  }
  else if(name=="tni.what"){
    opts<-c("gexp","tfs","para","refnet","tnet","refregulons","regulons","cdt", "summary","status","regulons.and.mode","refregulons.and.mode")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ) )
  }
  else if(name=="tni.gtype"){
    opts<-c("rmap","amap","mmap")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options: \n", paste(opts,collapse = ", ") ) )
  }
  else if(name=="tna.gtype"){
    opts<-c("rmap","amap")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options: \n", paste(opts,collapse = ", ") ) )
  }
  else if(name=="order") {
    if(!is.logical(para) || length(para)!=1)
      stop("'order' should be a logical value!\n")
  }
  else if(name=="regulon.order") {
    opts<-c('name','score','size', 'pvalue','adj.pvalue','none')
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'regulon.order' should be any one of the options: \n", paste(opts,collapse = ", ") ) )
    para<-switch(which(opts%in%para),'Regulon','Observed.Score','Regulon.Size','Pvalue','Adjusted.Pvalue','none')
    return(para)
  }
  else if(name=="ntop") {
    if(!is.null(para) && ( !(is.numeric(para) || is.integer(para)) || length(para)!=1 || round(para,0)!=para) )
      stop("'ntop' should be an integer value!\n")
  }
  else if(name=="estimator"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("pearson", "kendall", "spearman")))
      stop("'estimator' should be any one of 'pearson', 'kendall' and 'spearman'!\n")
  }
  else if(name=="reportNames"){
    if(!is.logical(para) || length(para)!=1)
      stop("'reportNames' should be a logical value!\n")
  }
  else if(name=="autoformat"){
    if(!is.logical(para) || length(para)!=1)
      stop("'autoformat' should be a logical value!\n")
  }
  else if(name=="splitcor"){
    if(!is.logical(para) || length(para)!=1)
      stop("'splitcor' should be a logical value!\n")
  }
  else if(name=="labpair") {
    if( (!is.matrix(para) || ncol(para)!=2) )
      stop(" 'labpair' should be a two-column matrix with regulon pair labels!\n ")
  }
  else if(name=="tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!\n")
  }  
  else if(name=="gsea.tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi","ref","cdt")))
      stop("'tnet' should be any one of 'dpi', 'ref' and 'cdt'!\n")
  }
  else if(name=="TNI"){
    if(!class(para)=="TNI")
      stop("'object' should be an object of class 'TNI'!\n")
    if( sum(names(para@results)%in%c("tn.ref", "tn.dpi") )<2 )
      stop("'object' should contain results in the slot 'result'!\n")
    if(!is.character(para@transcriptionFactors) || 
      any(is.na(para@transcriptionFactors)) || any(para@transcriptionFactors==""))
      stop(" 'object' should contain TF names in the slot 'transcriptionFactors'!\n ")
  }
  else if(name=="transcriptionalNetwork") {
    if( !is.matrix(para) )
      stop(" 'transcriptionalNetwork' should be a matrix with targets on rows and TFs on cols!\n ")
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("'transcriptionalNetwork' should be a matrix of named rows and cols (unique names)!\n")  
    if(sum(is.na(para))>0)
      stop("'transcriptionalNetwork' matrix should have no NAs!\n")
  }
  else if(name=="referenceNetwork") {
    if( !is.matrix(para) )
      stop(" 'referenceNetwork' should be a matrix with targets on rows and TFs on cols!\n ")
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("'referenceNetwork' should be a matrix of named rows and cols (unique names)!\n")  
    if(sum(is.na(para))>0)
      stop("'referenceNetwork' matrix should have no NAs!\n")
  }  
  else if(name=="gexp") {
    if( !is.matrix(para) || !is.numeric(para[1,]))
      stop(" 'gexp' should be a numeric matrix with genes on rows and samples on cols!\n ")
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("the 'gexp' matrix should be named on rows and cols (unique names)!\n")
  }
  else if(name=="transcriptionFactors") {
    if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
      stop(" 'transcriptionFactors' should be a character vector, without NA or empty names!\n ")
    if(length(unique(para)) < length(para))
      stop(" 'transcriptionFactors' should have unique identifiers!\n ")
  }
  else if(name=="tfs") {
    if(!is.null(para)){
    if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
      stop(" 'tfs' should be a character vector, without NA or empty names!\n ")
    }
  }
  else if(name=="modulators") {
    if(!is.null(para)){
      if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
        stop(" 'modulators' should be a character vector, without NA or empty names!\n ") 
    }
  }
  else if(name=="sampling") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1 || round(para,0)!=para )
      stop("'sampling' should be an integer value >= 1 and <= 100!\n")
  }
  else if(name=="amapCutoff") {
    if(!is.null(para) && (!is.numeric(para) || length(para)!=1 || para>1 || para<0))
      stop("'amapCutoff' should be a numeric value >=0 and <=1!\n")
  }
  else if(name=="amapFilter"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper","quantile","custom")))
      stop("'amapFilter' should be any one of 'quantile', 'phyper' and 'custom'!\n")
  }
  else if(name=="rgc") {
    if(!is.list(para))
      stop("'listOfRegulons' should be a list of gene sets!\n")
    if(  is.null(names(para)) || length(unique(names(para)))<length(names(para))  )
      stop("'listOfRegulons' should be a list of named gene sets (unique names)!\n")
    if(any(unlist(lapply(para,length))==0))
      stop("Empty regulon(s) found: see 'listOfRegulons' derived from the transcriptional network!\n")
  }
  else if(name=="rgcs") {
    if(!is.list(para))
      stop("'listOfRegulonPairs' should be a list of gene set collections!\n")
    if(is.null(names(para)))
      stop("'listOfRegulonPairs' should be a list of named gene set collections!\n")
    if(!all(unlist(lapply(para,is.list))))
      stop("Each gene set collection in 'listOfRegulonPairs' should be a list of gene sets!\n")
    if(any(unlist(lapply(para,length))==0))
      stop("Empty gene set collection(s) in 'listOfRegulonPairs'!\n")
  }
  else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
      !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!\n")
  }
  else if(name=="coverage") {
    if(!is.character(para) || length(para)!=1 || 
         !(para %in% c("tnet", "all")))
      stop("'coverage' should be any one of 'tnet' and 'all'!\n")
  }
  else if(name=="pValueCutoff") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'pValueCutoff' should be an integer or numeric value <=1 and >=0!\n")
  }
  else if(name=="miThreshold") {
    if(!is.null(para)){
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0)
      stop("'miThreshold' should be an integer or numeric value >=0!\n")
    }
  } 
  else if(name=="consensus") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1)
      stop("'consensus' should be an integer or numeric value <=100 and >=1!\n")
  }
  else if(name=="pooledNullDistribution") {
    if(!is.logical(para) || length(para)!=1)
      stop("'pooledNullDistribution' should be a logical value!\n")
  }
  else if(name=="nPermutations") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nPermutations' should be an integer >=1 !\n'")
  }
  else if(name=="parChunks") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'parChunks' should be an integer >=1 !\n'")
  }
  else if(name=="nBootstraps") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nBootstraps' should be an integer >=1 !\n'")
  }
  else if(name=="minRegulonSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'minRegulonSize' should be an integer >=1 !\n'") 
  }
  else if(name=="minIntersectSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || para>100)
      stop("'minIntersectSize' should be an integer >=1 and <=100!\n'")
  }
  else if(name=="statFilter") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper", "pchisq")))
      stop("'statFilter' should be any one of 'phyper' and 'pchisq'!\n")
  }
  else if(name=="statUniverse") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("all", "tnet")))
      stop("'statUniverse' should be any one of 'all' and 'tnet'!\n")
  }
  else if(name=="eps") {
    if(!is.numeric(para) || length(para)!=1 || para<0)
      stop("'eps' should be an numeric value >=0 !\n'")   
  }
  else if(name=="exponent") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1)
      stop("'exponent' should be an integer or numeric value >=1 !\n")
  }  
  if(name=="cvfilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'cvfilter' should be a logical value!\n")
  }
  else if(name=="verbose") {
    if(!is.logical(para) || length(para)!=1)
      stop("'verbose' should be a logical value!\n")
  }
  else if(name=="tnet") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!\n")
  }
  else if(name=="stepFilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'stepFilter' should be a logical value!\n")
  } 
  else if(name=="orderAbsValue") {
    if(!is.logical(para) || length(para)!=1)
      stop("'orderAbsValue' should be a logical value!\n")
  }
  else if(name=="phenotype") {
    if(!is.null(para)){
      if( !( is.numeric(para) || is.integer(para) ) || length(para)==0 )
        stop("'phenotype' should be a named numeric or integer vector with length > 0!\n")
      if(is.null(names(para)) || any(is.na(names(para))) || any(names(para)==""))
        stop("'phenotype' should be a named vector, without NA or empty names!\n")      
    }
  }
  else if(name=="hits") {
    if(!is.null(para)){
      if(!is.vector(para) || length(para)==0 )
        stop("'hits' should be a character vector with length > 0!\n")      
    }
  }
  else if(name=="phenoIDs") {
    if(is.null(para)){
      return(para)
    } else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2 ){
        stop("'phenoIDs' should be a matrix of characters with length >= 2!\n")
      }
      para[,1]<-as.character(para[,1])
      para[,2]<-as.character(para[,2])
      rownames(para)<-as.character(para[,1])
      if( (any(is.na(para[,1])) || any(para[,1]=="") ) ){
        stop("Col 1 in 'phenoIDs' matrix should have no NA or empty names!\n")      
      }
      return(as.matrix(para)) 
    }
  }
  else if(name=="gexpIDs"){
    if(is.null(para)){
      return(para)
    } else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2){
        stop("'gexpIDs' should be a data frame (or a matrix of characters) with ncol >= 2!\n")
      }
      sapply(1:ncol(para),function(i){
        para[,i]<<-as.character(para[,i])
      })
      if( (any(is.na(para[,1])) || any(para[,1]=="") ) ){
        stop("Col 1 in 'gexpIDs' matrix should have no empty names or NAs!\n")      
      }
      rownames(para)<-para[,1]
      para<-data.frame(para,stringsAsFactors=FALSE)
      #check colnames
      if(is.null(colnames(para))){
        colnames(para)<-paste("C",1:ncol(para),sep="")
      }
      return(para)
    }
  }
  else if(name=="duplicateRemoverMethod") {
    if(!is.character(para) || length(para) != 1 || !(para %in% c("max","min","average")))
      stop("'duplicateRemoverMethod' should be only one of the following character strings: 'max', 'min' and 'average'")
  }
  else if(name=="gs") {
    if(!is.character(para) || length(para)==0 || any(is.na(para)) || any(para==""))
      stop("'geneSet/regulon' should be a character vector with length > 0, without NA or empty names!\n")
  } 
  else if(name=="universe") {
    if(!is.character(para) || any(is.na(para)) || any(para=="")) 
      stop("'universe' should be a character vector without any NA or empty values!\n")
  }  
  else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
      !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!\n")
  }
  else if(name=="globalAdjustment") {
    if(!is.logical(para) || length(para)!=1)
      stop("'globalAdjustment' should be a logical value!\n")
  }
  else if(name=="tna.para") {
    if(missing(para))
      stop("'para' should be provided as a list!\n")
    ##default parameters
    para.default<-list(pValueCutoff = 0.05,pAdjustMethod = "BH", nPermutations = 1000, minRegulonSize = 15,exponent = 1)
    ##check if input parameters are supported
    if(length(setdiff(names(para),names(para.default)))>0) 
      stop("Some parameters in 'para' are not supported. Check the right format of para!\n")
    ##fill out default parameters for non-specified ones
    para.unspecified<-setdiff(names(para.default),names(para))
    if(length(para.unspecified)>0)
      for(i in 1:length(para.unspecified)) {
        para[[para.unspecified[i]]]<-para.default[[para.unspecified[i]]]
      }
    ##check data type in para
    if(!(is.integer(para$pValueCutoff) || is.numeric(para$pValueCutoff)) || length(para$pValueCutoff)!=1 || para$pValueCutoff>1)
      stop("'pValueCutoff' should be an integer or numeric value <=1!\n")
    if(!is.character(para$pAdjustMethod) || length(para$pAdjustMethod)!=1 || 
      !(para$pAdjustMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!\n")
    if(!(is.integer(para$nPermutations) || is.numeric(para$nPermutations)) || length(para$nPermutations)!=1 || para$nPermutations<1)
      stop("'nPermutations' should be an integer >=1 !\n'")
    if(!(is.integer(para$minRegulonSize) || is.numeric(para$minRegulonSize)) || length(para$minRegulonSize)!=1 || para$minRegulonSize<1)
      stop("'minRegulonSize' should be an integer >=1 !\n'")
    if(!(is.integer(para$exponent) || is.numeric(para$exponent)) || length(para$pValueCutoff)!=1 || para$exponent<1)
      stop("'exponent' should be an integer or numeric value >=1 !\n")
    return(para)
  }
  else if(name=="tni.para") {
    if(missing(para))stop("'para' should be provided as a list!\n")
    ##default parameters
    para.default<-list(pValueCutoff=0.001, pAdjustMethod="BH", globalAdjustment=TRUE, estimator="pearson",
                       nPermutations=1000, nBootstraps=100, pooledNullDistribution=TRUE, consensus=95, eps=0)
    ##check if input parameters are supported
    if(length(setdiff(names(para),names(para.default)))>0) 
      stop("Some parameters in 'para' are not supported. Check the right format of para!\n")
    ##fill out default parameters for non-specified ones
    para.unspecified<-setdiff(names(para.default),names(para))
    if(length(para.unspecified)>0)
      for(i in 1:length(para.unspecified)) {
        para[[para.unspecified[i]]]<-para.default[[para.unspecified[i]]]
      }
    ##check data type in para
    if(!(is.integer(para$pValueCutoff) || is.numeric(para$pValueCutoff)) || 
      length(para$pValueCutoff)!=1 || para$pValueCutoff>1)
      stop("'pValueCutoff' should be an integer or numeric value <=1!\n")
    if(!is.character(para$pAdjustMethod) || length(para$pAdjustMethod)!=1 || 
      !(para$pAdjustMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!\n")
    if(!is.logical(para$globalAdjustment) || length(para$globalAdjustment)!=1)
      stop("'globalAdjustment' should be a logical value!\n")
    if(!is.character(para$estimator) || length(para$estimator)!=1 || 
      !(para$estimator %in% c("pearson", "kendall", "spearman")))
      stop("'estimator' should be any one of 'pearson', 'kendall' and 'spearman'")
    if(!(is.integer(para$nPermutations) || is.numeric(para$nPermutations)) || length(para$nPermutations)!=1 || para$nPermutations<1)
      stop("'nPermutations' should be an integer >=1 !\n'")
    if(!(is.integer(para$nBootstraps) || is.numeric(para$nBootstraps)) || length(para$nBootstraps)!=1 || para$nBootstraps<1)
      stop("'nBootstraps' should be an integer >=1 !\n'")
    if(!(is.integer(para$consensus) || is.numeric(para$consensus)) || length(para$consensus)!=1 || 
      para$consensus>100 || para$consensus<1)
      stop("'consensus' should be an integer or numeric value <=100 and >=1!\n")   
    if(!is.logical(para$pooledNullDistribution) || length(para$pooledNullDistribution)!=1)
      stop("'pooledNullDistribution' should be a logical value!\n")
    if(!is.numeric(para$eps) || length(para$eps)!=1 || para$eps<0)
      stop("'eps' should be an numeric value >=0 !\n'")    
    return(para)
  }
  else if(name=="filepath"){
    if(!is.character(para) || length(para)!=1)
      stop("'filepath' should be a single character!\n")
  }
  else if(name=="ylimPanels") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)!=4 )
      stop("'ylimPanels' should be a numeric vector with length=4 !\n")
  }
  else if(name=="heightPanels") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)!=3 )
      stop("'heightPanels' should be a numeric vector with length=3 !\n")
  }
  else if(name=="width") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 )
      stop("'width' should be an integer or numeric value!\n")
  }  
  else if(name=="height") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 )
      stop("'height' should be an integer or numeric value!\n")
  }
  else if(name=="ylabPanels") {
    if(!is.character(para) || length(para)!=3)
      stop(" 'ylabPanels' should be a character vector with length=3 !\n ")
  }
  else if(name=="xlab") {
    if(!is.character(para) || length(para)!=1)
      stop(" 'xlab' should be a character with length=1 !\n ")
  }  
  else if(name=="labPheno") {
    if(!is.character(para) || length(para)!=1)
      stop(" 'labPheno' should be a character with length=1 !\n ")
  }
}
