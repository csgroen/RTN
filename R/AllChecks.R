##This function is used for argument checking
tnai.checks <- function(name, para) {
  if(name=="alpha") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'alpha' should be an integer or numeric value <=1 and >=0 !\n",call.=FALSE)
  }
  else if(name=="maxgap") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 || as.integer(para)<para)
      stop("'maxgap' should be an integer value >=0!\n",call.=FALSE)
  }
  else if(name=="tna.what"){
    opts<-c("tnet","tfs","pheno","regulons","refregulons","para","mra","gsea1","gsea2","overlap","synergy","shadow","summary","status",
            "regulons.and.pheno","refregulons.and.pheno","regulons.and.mode","refregulons.and.mode")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
  else if(name=="tni.what"){
    opts<-c("gexp","tfs","para","refnet","tnet","refregulons","regulons","cdt","cdtrev", "summary","status",
            "regulons.and.mode","refregulons.and.mode")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ) ,call.=FALSE )
  }
  else if(name=="avs.what"){
    opts<-c("markers","validatedMarkers","variantSet","randomSet","summary",
            "status","linkedMarkers","randomMarkers","vse","evse")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
  else if(name=="avs.plot.what"){
    opts<-c("vse","evse")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'avs.plot.what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
  else if(name=="reldata"){
    opts<-c("rel27CEU-NCBIB36")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("available 'reldata' option:", paste(opts,collapse = ", ") ),call.=FALSE)
  }
  else if(name=="ldfilter"){
    opts<-c("DprimeLOD","Rsquare")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("available 'ldfilter' option:", paste(opts,collapse = ", ") ),call.=FALSE)
  }
  else if(name=="snpop"){
    if(!is.null(para)){
      if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
        stop("'snpop' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    }
  } 
  else if(name=="markers") {
    if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
      stop("'markers' should be a character vector, without 'NA' or empty names!",call.=FALSE)
  }
  else if(name=="tni.gtype"){
    opts<-c("rmap","amap","mmap","mmapDetailed")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
  }
  else if(name=="tna.gtype"){
    opts<-c("rmap","amap")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
  }
  else if(name=="order") {
    if(!is.logical(para) || length(para)!=1)
      stop("'order' should be a logical value!",call.=FALSE)
  }
  else if(name=="regulon.order") {
    opts<-c('name','score','size', 'pvalue','adj.pvalue','none')
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'regulon.order' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
    para<-switch(which(opts%in%para),'Regulon','Observed.Score','Regulon.Size','Pvalue','Adjusted.Pvalue','none')
    return(para)
  }
  else if(name=="ntop") {
    if(!is.null(para) && ( !(is.numeric(para) || is.integer(para)) || length(para)!=1 || round(para,0)!=para) )
      stop("'ntop' should be an integer value!",call.=FALSE)
  }
  else if(name=="nsplit") {
    if(!is.null(para) && ( !(is.numeric(para) || is.integer(para)) || length(para)!=1 || round(para,0)!=para) )
      stop("'nsplit' should be an integer value!",call.=FALSE)
  }
  else if(name=="estimator"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("pearson", "kendall", "spearman")))
      stop("'estimator' should be any one of 'pearson', 'kendall' and 'spearman'!",call.=FALSE)
  }
  else if(name=="pwtransform"){
    if(!is.logical(para) || length(para)!=1)
      stop("'pwtransform' should be a logical value!",call.=FALSE)
  }
  else if(name=="reportNames"){
    if(!is.logical(para) || length(para)!=1)
      stop("'reportNames' should be a logical value!",call.=FALSE)
  }
  else if(name=="autoformat"){
    if(!is.logical(para) || length(para)!=1)
      stop("'autoformat' should be a logical value!",call.=FALSE)
  }
  else if(name=="splitcor"){
    if(!is.logical(para) || length(para)!=1)
      stop("'splitcor' should be a logical value!",call.=FALSE)
  }
  else if(name=="labpair") {
    if( (!is.matrix(para) || ncol(para)!=2) )
      stop("'labpair' should be a two-column matrix with regulon pair labels!",call.=FALSE)
  }
  else if(name=="tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }  
  else if(name=="gsea.tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi","ref","cdt")))
      stop("'tnet' should be any one of 'dpi', 'ref' and 'cdt'!",call.=FALSE)
  }
  else if(name=="TNI"){
    if(!class(para)=="TNI")
      stop("'object' should be an object of class 'TNI'!",call.=FALSE)
    if( sum(names(para@results)%in%c("tn.ref", "tn.dpi") )<2 )
      stop("'object' should contain results in the slot 'result'!",call.=FALSE)
    if(!is.character(para@transcriptionFactors) || 
      any(is.na(para@transcriptionFactors)) || any(para@transcriptionFactors==""))
      stop("'object' should contain TF names in the slot 'transcriptionFactors'!",call.=FALSE)
  }
  else if(name=="transcriptionalNetwork") {
    if( !is.matrix(para) )
      stop(" 'transcriptionalNetwork' should be a matrix with targets on rows and TFs on cols!",call.=FALSE)
    if( is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("'transcriptionalNetwork' should be a matrix of named rows and cols (unique names)!",call.=FALSE)  
    if(sum(is.na(para))>0)
      stop("'transcriptionalNetwork' matrix should have no NAs!",call.=FALSE)
  }
  else if(name=="referenceNetwork") {
    if( !is.matrix(para) )
      stop("'referenceNetwork' should be a matrix with targets on rows and TFs on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("'referenceNetwork' should be a matrix of named rows and cols (unique names)!",call.=FALSE)  
    if(sum(is.na(para))>0)
      stop("'referenceNetwork' matrix should have no NAs!",call.=FALSE)
  }  
  else if(name=="gexp") {
    if( !is.matrix(para) || !is.numeric(para[1,]))
      stop("'gexp' should be a numeric matrix with genes on rows and samples on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
      length(unique(rownames(para))) < length(rownames(para)) ||
      length(unique(colnames(para))) < length(colnames(para)))
      stop("the 'gexp' matrix should be named on rows and cols (unique names)!",call.=FALSE)
  }
  else if(name=="gxdata") {
    if( !is.matrix(para) || !is.numeric(para[1,]))
      stop("'gxdata' should be a numeric matrix with genes on rows and samples on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
           length(unique(rownames(para))) < length(rownames(para)) ||
           length(unique(colnames(para))) < length(colnames(para)))
      stop("the 'gxdata' matrix should be named on rows and cols (unique names)!",call.=FALSE)
  }
  else if(name=="transcriptionFactors") {
    if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
      stop("'transcriptionFactors' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    if(length(unique(para)) < length(para))
      stop("'transcriptionFactors' should have unique identifiers!",call.=FALSE)
  }
  else if(name=="tfs") {
    if(!is.null(para)){
    if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
      stop("'tfs' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    }
  }
  else if(name=="modulators") {
    if(!is.null(para)){
      if( !(is.character(para) || is.numeric(para)) || any(is.na(para)) || any(para==""))
        stop("'modulators' should be a character vector, without 'NA' or empty names!",call.=FALSE) 
    }
  }
  else if(name=="sampling") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1 || round(para,0)!=para )
      stop("'sampling' should be an integer value >= 1 and <= 100 !",call.=FALSE)
  }
  else if(name=="probs") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'probs' should be a numeric value >= 0 and <= 1 !",call.=FALSE)
  }
  else if(name=="prob") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'prob' should be a numeric value >= 0 and <= 1 !",call.=FALSE)
  }
  else if(name=="amapCutoff") {
    if(!is.null(para) && (!is.numeric(para) || length(para)!=1 || para>1 || para<0))
      stop("'amapCutoff' should be a numeric value >=0 and <=1!",call.=FALSE)
  }
  else if(name=="amapFilter"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper","quantile","custom")))
      stop("'amapFilter' should be any one of 'quantile', 'phyper' and 'custom'!",call.=FALSE)
  }
  else if(name=="listOfRegulons") {
    if(!is.list(para))
      stop("'listOfRegulons' should be a list of gene sets!\n")
    if(  is.null(names(para)) || length(unique(names(para)))<length(names(para))  )
      stop("'listOfRegulons' should be a named list (unique names)!",call.=FALSE)
    #if(any(unlist(lapply(para,length))==0))
    #  stop("empty regulon(s) found!",call.=FALSE)
  }
  else if(name=="listOfRegulonPairs") {
    if(!is.list(para))
      stop("'listOfRegulonPairs' should be a list of gene set collections!",call.=FALSE)
    if(is.null(names(para)))
      stop("'listOfRegulonPairs' should be a named list!",call.=FALSE)
    if(!all(unlist(lapply(para,is.list))))
      stop("each gene set collection in 'listOfRegulonPairs' should be a list of gene sets!",call.=FALSE)
    #if(any(unlist(lapply(para,length))==0))
    #  stop("empty gene set collection(s) in 'listOfRegulonPairs'!",call.=FALSE)
  }
  else if(name=="glist") {
    if(!is.null(para)){
      if(!is.list(para))
        stop("'glist' should be a list of gene sets or regulons!\n")
      if(  is.null(names(para)) || length(unique(names(para)))<length(names(para))  )
        stop("'glist' should be a named list (unique names)!",call.=FALSE)
      if(any(!unlist(lapply(para,class))=="character")){
        para<-lapply(para,function(gl){
          as.character(gl)
        })
      }
      #if(any(unlist(lapply(para,length))==0))
      #  stop("empty gene set(s) found in 'glist'!",call.=FALSE)
      return(para)
    }
  }
  else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
      !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm','hochberg','hommel','bonferroni','BH','BY','fdr' and 'none'!",call.=FALSE)
  }
  else if(name=="coverage") {
    if(!is.character(para) || length(para)!=1 || 
         !(para %in% c("tnet", "all")))
      stop("'coverage' should be any one of 'tnet' and 'all'!",call.=FALSE)
  }
  else if(name=="pValueCutoff") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'pValueCutoff' should be an integer or numeric value <=1 and >=0 !",call.=FALSE)
  }
  else if(name=="miThreshold") {
    if(is.character(para)){
      if(!is.character(para) || length(para)!=1 || !(para %in% c("md","md.tf","md.tf.tar")))
        stop("'miThreshold' should be any one of 'md', 'md.tf' and 'md.tf.tar'!",call.=FALSE)
    } else {
      if( !(is.integer(para) || is.numeric(para)) || any(para<0) )
        stop("custom 'miThreshold' should have integer or numeric values >=0 !",call.=FALSE)
    }
  }
  else if(name=="consensus") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1)
      stop("'consensus' should be an integer or numeric value <=100 and >=1 !",call.=FALSE)
  }
  else if(name=="pooledNullDistribution") {
    if(!is.logical(para) || length(para)!=1)
      stop("'pooledNullDistribution' should be a logical value!",call.=FALSE)
  }
  else if(name=="nPermutations") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nPermutations' should be an integer >=1 !",call.=FALSE)
  }
  else if(name=="nrand") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nrand' should be an integer >=1 !",call.=FALSE)
  }
  else if(name=="parChunks") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'parChunks' should be an integer >=1 !",call.=FALSE)
  }
  else if(name=="nBootstraps") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nBootstraps' should be an integer >=1 !",call.=FALSE)
  }
  else if(name=="minRegulonSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'minRegulonSize' should be an integer >=1 !",call.=FALSE) 
  }
  else if(name=="minSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'minSize' should be an integer >=1 !",call.=FALSE) 
  }
  else if(name=="evse.minSize") {
    if(!(is.integer(para) || is.numeric(para)) || any(para<1) || any(round(para,0)!=para) )
      stop("'minSize' should be integer >=1 !",call.=FALSE)
    if( length(para)<1 || length(para)>2)
      stop("'minSize' should either be a single value or a vector of length = 2 !",call.=FALSE)
    if(length(para)==1)para<-c(para,para)
    return(para)
  }
  else if(name=="minIntersectSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || para>100)
      stop("'minIntersectSize' should be an integer >=1 and <=100 !",call.=FALSE)
  }
  else if(name=="minAgreement") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || para>100)
      stop("'minAgreement' should be an integer >=1 and <=100 !",call.=FALSE)
  }
  else if(name=="statFilter") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper","pchisq")))
      stop("'statFilter' should be any one of 'phyper' and 'pchisq'!",call.=FALSE)
  }
  else if(name=="statUniverse") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("all","tnet","reg")))
      stop("'statUniverse' should be any one of 'all' and 'tnet'!",call.=FALSE)
  }
  else if(name=="eps") {
    if(!is.numeric(para) || length(para)!=1 || para<0)
      stop("'eps' should be an numeric value >=0!",call.=FALSE)   
  }
  else if(name=="exponent") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1)
      stop("'exponent' should be an integer or numeric value >=1 !",call.=FALSE)
  }  
  else if(name=="cvfilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'cvfilter' should be a logical value!",call.=FALSE)
  }
  else if(name=="verbose") {
    if(!is.logical(para) || length(para)!=1)
      stop("'verbose' should be a logical value!",call.=FALSE)
  }
  else if(name=="mask") {
    if(!is.logical(para) || length(para)!=1)
      stop("'mask' should be a logical value!",call.=FALSE)
  }
  else if(name=="decreasing") {
    if(!is.logical(para) || length(para)!=1)
      stop("'decreasing' should be a logical value!",call.=FALSE)
  }
  else if(name=="fineMapping") {
    if(!is.logical(para) || length(para)!=1)
      stop("'fineMapping' should be a logical value!",call.=FALSE)
  }
  else if(name=="boxcox") {
    if(!is.logical(para) || length(para)!=1)
      stop("'boxcox' should be a logical value!",call.=FALSE)
  }
  else if(name=="mergeColinked") {
    if(!is.logical(para) || length(para)!=1)
      stop("'mergeColinked' should be a logical value!",call.=FALSE)
  }
  else if(name=="tnet") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }
  else if(name=="stepFilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'stepFilter' should be a logical value!",call.=FALSE)
  } 
  else if(name=="orderAbsValue") {
    if(!is.logical(para) || length(para)!=1)
      stop("'orderAbsValue' should be a logical value!",call.=FALSE)
  }
  else if(name=="phenotype") {
    if(!is.null(para)){
      if( !( is.numeric(para) || is.integer(para) ) || length(para)==0 )
        stop("'phenotype' should be a named numeric or integer vector with length >0 !",call.=FALSE)
      if(is.null(names(para)) || any(is.na(names(para))) || any(names(para)==""))
        stop("'phenotype' should be a named vector, without 'NA' or empty names!",call.=FALSE)      
    }
  }
  else if(name=="hits") {
    if(!is.null(para)){
      if(!is.vector(para) || length(para)==0 )
        stop("'hits' should be a character vector with length >0 !",call.=FALSE)      
    }
  }
  else if(name=="phenoIDs") {
    if(is.null(para)){
      return(para)
    } else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2 ){
        stop("'phenoIDs' should be a dataframe (or a matrix of characters) with ncol >=2 !",call.=FALSE)
      }
      junk<-sapply(1:ncol(para),function(i){
        tp<-para[,i]
        if(is.list(tp))tp<-unlist(tp)
        if(i<=2 || is.factor(tp))tp<-as.character(tp)
        para[,i]<<-tp
      })
      rownames(para)<-para[,1]
      if( (any(is.na(para[,1])) || any(para[,1]=="") ) ){
        stop("Col 1 in 'phenoIDs' matrix should have no 'NA' or empty values!",call.=FALSE)      
      }
      para<-as.data.frame(para,stringsAsFactors=FALSE)
      return(para)
    }
  }
  else if(name=="annotation"){
    if( !is.data.frame(para) || ncol(para)!=4 ){
      stop("'annotation' should be a dataframe with ncol=4 !",call.=FALSE)
    }
    colnames(para)<-toupper(colnames(para))
    cnames<-c("CHROM","START","END")
    if( !all(cnames%in%colnames(para)) ){
      tp1<-"\nCol names in the 'annotation' data do not match valid labels!\n"
      tp2<-"Please, revise the 'annotation' format:\n"
      tp3<-"Name col1: <CHROM> chromosome name\n"
      tp4<-"Name col2: <START> start position\n"
      tp5<-"Name col3: <END> end position\n"
      tp6<-"Name col4: <ID> any genomic id or name of the line\n"
      stop(tp1,tp2,tp3,tp4,tp5,tp6,call.=FALSE)
    }
    idx<-c(which(!colnames(para)%in%cnames),match(cnames,colnames(para)))
    para<-para[,idx,drop=FALSE]
    #---
    idlab<-colnames(para)[1]
    if(any(duplicated(para[,1]))){
      stop(paste("Col '",idlab,"' in the 'annotation' data should have non-duplicated IDs!",sep=""),call.=FALSE)
    }
    colnames(para)[1]<-"ID"
    #---
    sapply(1:ncol(para),function(i){
      tp<-para[,i]
      if( any( is.na(tp) || any(tp=="") ) ){
        stop("'annotation' matrix should have no 'NA' or empty values!",call.=FALSE)      
      }
      if(is.list(tp))para[,i]<<-unlist(tp)
    })
    #---
    para$ID<-as.character(para$ID)
    para$CHROM<-as.character(para$CHROM)
    chrs<-c(paste("chr",1:22,sep=""),"chrX","chrY")
    chrChecks<-!para$CHROM%in%chrs
    if( any(chrChecks) ){
      n<-sum(chrChecks)/length(chrChecks)
      tp1<-paste("chromosome values in 'annotation' should be listed in ", sep="")
      tp2<-"[chr1, chr2, chr3, ..., chr22, chrX]!"
      if(n>0.95){
        stop(tp1,tp2,call.=FALSE)
      } else {
        nonvalid<-paste(unique(para$CHROM[chrChecks]),collapse=", ")
        tp3<-"\n...the following values were removed: "
        warning(tp1,tp2,tp3,nonvalid,call.=FALSE)
        para<-para[!chrChecks,]
      }
    }
    #---
    b1<-is.integer(para$START) || is.numeric(para$START)
    b2<-is.integer(para$END) || is.numeric(para$END)
    if( !(b1 && b2) ){
      stop("Chromosome start/end positions in 'annotation' should be integer or numeric vectors!",call.=FALSE)
    }
    rownames(para)<-para$ID
    return(para)
  }
  else if(name=="gexpIDs"){
    if(is.null(para)){
      return(para)
    } else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2){
        stop("'gexpIDs' should be a data frame (or a matrix of characters) with ncol >=2 !",call.=FALSE)
      }
      junk<-sapply(1:ncol(para),function(i){
        tp<-para[,i]
        if(is.list(tp))tp<-unlist(tp)
        tp<-as.character(tp)
        para[,i]<<-tp
      })
      if( (any(is.na(para[,1])) || any(para[,1]=="") ) ){
        stop("Col 1 in 'gexpIDs' should have no NA or empty value!",call.=FALSE)      
      }
      if( length(unique(para[,1])) < length(para[,1]) )
        stop("Col 1 in 'gexpIDs' matrix should have unique ids!",call.=FALSE)
      rownames(para)<-para[,1]
      para<-data.frame(para,stringsAsFactors=FALSE)
      #check colnames
      if(is.null(colnames(para))){
        colnames(para)<-paste("C",1:ncol(para),sep="")
      }
      colnames(para)<-toupper(colnames(para))
      if(any(duplicated(colnames(para)))){
        stop("'gexpIDs' matrix should have unique col names (not case sensitive)!")
      }
      return(para)
    }
  }
  else if(name=="duplicateRemoverMethod") {
    if(!is.character(para) || length(para) != 1 || !(para %in% c("max","min","average")))
      stop("'duplicateRemoverMethod' should be one of the following character strings: 'max', 'min' and 'average'",call.=FALSE)
  }
  else if(name=="gs") {
    if(!is.character(para) || length(para)==0 || any(is.na(para)) || any(para==""))
      stop("'geneSet/regulon' should be a character vector with length > 0, without 'NA' or empty names!",call.=FALSE)
  } 
  else if(name=="universe") {
    if(!is.character(para) || any(is.na(para)) || any(para=="")) 
      stop("'universe' should be a character vector without 'NA' or empty values!",call.=FALSE)
  }  
  else if(name=="globalAdjustment") {
    if(!is.logical(para) || length(para)!=1)
      stop("'globalAdjustment' should be a logical value!",call.=FALSE)
  }
  else if(name=="tna.para") {
    if(missing(para))
      stop("'para' should be provided as a list!",call.=FALSE)
    ##default parameters
    para.default<-list(pValueCutoff = 0.05,pAdjustMethod = "BH", nPermutations = 1000, minRegulonSize = 15,exponent = 1)
    ##check if input parameters are supported
    if(length(setdiff(names(para),names(para.default)))>0) 
      stop("Some parameters in 'para' are not supported. Check the right format of para!",call.=FALSE)
    ##fill out default parameters for non-specified ones
    para.unspecified<-setdiff(names(para.default),names(para))
    if(length(para.unspecified)>0)
      for(i in 1:length(para.unspecified)) {
        para[[para.unspecified[i]]]<-para.default[[para.unspecified[i]]]
      }
    ##check data type in para
    if(!(is.integer(para$pValueCutoff) || is.numeric(para$pValueCutoff)) || length(para$pValueCutoff)!=1 || para$pValueCutoff>1)
      stop("'pValueCutoff' should be an integer or numeric value <=1!",call.=FALSE)
    if(!is.character(para$pAdjustMethod) || length(para$pAdjustMethod)!=1 || 
      !(para$pAdjustMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!",call.=FALSE)
    if(!(is.integer(para$nPermutations) || is.numeric(para$nPermutations)) || length(para$nPermutations)!=1 || para$nPermutations<1)
      stop("'nPermutations' should be an integer >=1!",call.=FALSE)
    if(!(is.integer(para$minRegulonSize) || is.numeric(para$minRegulonSize)) || length(para$minRegulonSize)!=1 || para$minRegulonSize<1)
      stop("'minRegulonSize' should be an integer >=1!",call.=FALSE)
    if(!(is.integer(para$exponent) || is.numeric(para$exponent)) || length(para$pValueCutoff)!=1 || para$exponent<1)
      stop("'exponent' should be an integer or numeric value >=1 !",call.=FALSE)
    return(para)
  }
  else if(name=="tni.para") {
    if(missing(para))stop("'para' should be provided as a list!",call.=FALSE)
    ##default parameters
    para.default<-list(pValueCutoff=0.001, pAdjustMethod="BH", globalAdjustment=TRUE, estimator="pearson",
                       nPermutations=1000, nBootstraps=100, pooledNullDistribution=TRUE, consensus=95, eps=0)
    ##check if input parameters are supported
    if(length(setdiff(names(para),names(para.default)))>0) 
      stop("Some parameters in 'para' are not supported. Check the right format of para!",call.=FALSE)
    ##fill out default parameters for non-specified ones
    para.unspecified<-setdiff(names(para.default),names(para))
    if(length(para.unspecified)>0)
      for(i in 1:length(para.unspecified)) {
        para[[para.unspecified[i]]]<-para.default[[para.unspecified[i]]]
      }
    ##check data type in para
    if(!(is.integer(para$pValueCutoff) || is.numeric(para$pValueCutoff)) || 
      length(para$pValueCutoff)!=1 || para$pValueCutoff>1)
      stop("'pValueCutoff' should be an integer or numeric value <=1 !",call.=FALSE)
    if(!is.character(para$pAdjustMethod) || length(para$pAdjustMethod)!=1 || 
      !(para$pAdjustMethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'!",call.=FALSE)
    if(!is.logical(para$globalAdjustment) || length(para$globalAdjustment)!=1)
      stop("'globalAdjustment' should be a logical value!",call.=FALSE)
    if(!is.character(para$estimator) || length(para$estimator)!=1 || 
      !(para$estimator %in% c("pearson", "kendall", "spearman")))
      stop("'estimator' should be any one of 'pearson', 'kendall' and 'spearman'",call.=FALSE)
    if(!(is.integer(para$nPermutations) || is.numeric(para$nPermutations)) || length(para$nPermutations)!=1 || para$nPermutations<1)
      stop("'nPermutations' should be an integer >=1!",call.=FALSE)
    if(!(is.integer(para$nBootstraps) || is.numeric(para$nBootstraps)) || length(para$nBootstraps)!=1 || para$nBootstraps<1)
      stop("'nBootstraps' should be an integer >=1!",call.=FALSE)
    if(!(is.integer(para$consensus) || is.numeric(para$consensus)) || length(para$consensus)!=1 || 
      para$consensus>100 || para$consensus<1)
      stop("'consensus' should be an integer or numeric value <=100 and >=1 !",call.=FALSE)   
    if(!is.logical(para$pooledNullDistribution) || length(para$pooledNullDistribution)!=1)
      stop("'pooledNullDistribution' should be a logical value!",call.=FALSE)
    if(!is.numeric(para$eps) || length(para$eps)!=1 || para$eps<0)
      stop("'eps' should be an numeric value >=0!",call.=FALSE)    
    return(para)
  }
  else if(name=="filepath"){
    if(!is.character(para) || length(para)!=1)
      stop("'filepath' should be a single character!",call.=FALSE)
  }
  else if(name=="ylimPanels") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)!=4 )
      stop("'ylimPanels' should be a numeric vector with length =4 !",call.=FALSE)
  }
  else if(name=="heightPanels") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)!=3 )
      stop("'heightPanels' should be a numeric vector with length =3 !",call.=FALSE)
  }
  else if(name=="bxseq") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)<2 )
      stop("'bxseq' should be a numeric vector with length >=2 !",call.=FALSE)
  }
  else if(name=="width") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'width' should be an integer or numeric value >0 !",call.=FALSE)
  }  
  else if(name=="height") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'height' should be an integer or numeric value >0 !",call.=FALSE)
  }
  else if(name=="maxy") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0)
      stop("'maxy' should be an integer or numeric value >0 !",call.=FALSE)
  }
  else if(name=="rmargin") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'rmargin' should be an integer or numeric value >0 !",call.=FALSE)
  }
  else if(name=="bxsp") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0)
      stop("'bxsp' should be an integer or numeric value >0 !",call.=FALSE)
  }
  else if(name=="ylabPanels"){
    if(!is.character(para) || length(para)!=3)
      stop("'ylabPanels' should be a character vector with length=3 !",call.=FALSE)
  }
  else if(name=="xlab"){
    if(!is.character(para) || length(para)!=1)
      stop("'xlab' should be a character value with length=1 !",call.=FALSE)
  }
  else if(name=="ylab"){
    if(!is.character(para) || length(para)!=1)
      stop("'ylab' should be a character value with length=1 !",call.=FALSE)
  }
  else if(name=="lab"){
    if(!is.character(para) || length(para)!=1)
      stop("'lab' should be a character value with length=1 !",call.=FALSE)
  }
  else if(name=="tlt"){
    if(!is.character(para) || length(para)!=1)
      stop("'tlt' should be a character value with length=1 !",call.=FALSE)
  }
  else if(name=="fname") {
    if(!is.character(para) || length(para)!=1)
      stop("'fname' should be a character value with length=1!",call.=FALSE)
  } 
  else if(name=="lab"){
    if(!is.character(para) || length(para)!=1)
      stop("'lab' should be a character value with length=1 !",call.=FALSE)
  }
  else if(name=="labPheno") {
    if(!is.character(para) || length(para)!=1)
      stop("'labPheno' should be a character value with length=1 !",call.=FALSE)
  } 
  else if(name=="idkey") {
    if(!is.null(para)){
      if(!is.character(para) || length(para)!=1)
        stop("'idkey' should be a character value with length=1 !",call.=FALSE)
    }
  }
  else {
    stop("...<",name,"> check is missing!!",call. = FALSE)
  }
}
