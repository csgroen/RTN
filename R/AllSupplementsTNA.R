

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##      OPTIMIZED GSEA FUNCTIONS FOR REGULONS AND TNETS
##       -- IMPLEMENTS AN INTERFACE FOR HTSanalyzeR --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##This function computes observed and permutation-based scores associated 
##with a gene set enrichment analysis for a collection of regulons.
gsea1tna <- function(listOfRegulons, phenotype, exponent=1, 
                     nPermutations=1000,verbose=TRUE) {	 
  #OBS1:names provided as integer values!
  #OBS2:max/min sizes checked in the previous functions!
  pheno.names <- as.integer(names(phenotype))
  nRegulons <- length(listOfRegulons)
  if(nRegulons > 0) {
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
    if(isParallel() && nRegulons>1) {
      if(verbose)cat("-Performing GSEA (parallel version - ProgressBar disabled)...\n")
      if(verbose)cat("--For", length(listOfRegulons), "regulons...\n")      
      scores <- gseaScoresBatchParallel4RTN(geneList=phenotype, geneNames.perm=perm.gL,
                                        collectionOfGeneSets=listOfRegulons,exponent=exponent,
                                        nPermutations=nPermutations)
      sapply(1:nRegulons, function(i){
        scoresperm[i,]<<-unlist(scores["scoresperm",i])
        scoresObserved[i]<<-unlist(scores["scoresObserved",i])
        NULL
      })
    } else {
      if(verbose) cat("-Performing gene set enrichment analysis...\n")
      if(verbose) cat("--For", length(listOfRegulons), "regulons...\n")
      if(verbose) pb <- txtProgressBar(style=3)
      for(i in 1:nRegulons) {
        scores <- gseaScoresBatch4RTN(geneList=phenotype, geneNames.perm=perm.gL, 
                                  geneSet=listOfRegulons[[i]],exponent=exponent,
                                  nPermutations=nPermutations)
        scoresObserved[i] <- scores$scoresObserved
        scoresperm[i,] <- scores$scoresperm
        if(verbose) setTxtProgressBar(pb, i/nRegulons)
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
gsea2tna <- function(listOfRegulons, phenotype, exponent=1, nPermutations=1000, 
                     verbose1=TRUE, verbose2=TRUE) {   
  #OBS1:names provided as integer values!
  #OBS2:max/min sizes checked in the previous functions!
  pheno.names <- as.integer(names(phenotype))
  nRegulons <- length(listOfRegulons)
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
    if(isParallel() && nRegulons>1) {
      if(verbose1 && verbose2)cat("-Performing two-tailed GSEA (parallel version - ProgressBar disabled)...\n")
      if(verbose1 && verbose2)cat("--For", length(listOfRegulons), "regulons...\n")      
      scores <- gseaScoresBatchParallel4RTN(geneList=phenotype, geneNames.perm = perm.gL,
                                            collectionOfGeneSets=listOfRegulons,
                                            exponent=exponent,nPermutations=nPermutations)
      sapply(1:nRegulons, function(i){
        scoresperm[i,]<<-unlist(scores["scoresperm",i])
        scoresObserved[i]<<-unlist(scores["scoresObserved",i])
        NULL
      })
    } else {
      if(verbose1 && verbose2) cat("-Performing two-tailed GSEA analysis...\n")
      if(verbose1 && verbose2) cat("--For", length(listOfRegulons), "regulons...\n")
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

##------------------------------------------------------------------------
##This function computes pairwise ES scores for synergy analysis
##..here synergy is inferred from the 'effectsize' by
##..comparing the enrichment score of the intesect related to the union!
##..i.e. the null distribution is computed by resampling the union!
pairwiseSynergy <- function(collectionsOfPairs, phenotype, exponent=1, 
                      nPermutations=1000, minIntersectSize=1, verbose=TRUE) {
  ##min intersect
  minter<-as.integer(minIntersectSize)*0.01
  if(isParallel() && length(collectionsOfPairs)>1) {
    if(verbose)cat("-Performing synergy analysis (parallel version - ProgressBar disabled)...\n")
    if(verbose)cat("--For", length(collectionsOfPairs), "regulon pairs...\n")
    #get permutation scores
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("gseaScores4RTN"),envir=environment())
    permScores<-parSapply(cl,1:length(collectionsOfPairs), function(i){
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
    if(verbose)cat("-Performing synergy analysis...\n")
    if(verbose)cat("--For", length(collectionsOfPairs), "regulon pairs...\n")
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
  ##min intersect
  minter<-as.integer(minIntersectSize)*0.01
  if( isParallel() && length(collectionsOfPairs)>1) {
    #get permutation scores
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("gseaScores4RTN"),envir=environment())
    permScores<-parSapply(cl,1:length(collectionsOfPairs), function(i){
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

##------------------------------------------------------------------------
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
  results <- matrix(NA, nrow=0, ncol=8)
  colnames(results) <- c("Regulon1","Regulon2","Universe.Size", "R1.Size", "R2.Size", 
                         "Expected.Overlap", "Observed.Overlap", "Pvalue")
  for(i in 1:length(regulon.filtered)){
      if(verbose) setTxtProgressBar(pb, i/length(regulon.filtered))
      r1<-regulon.filtered[i]
      r.filter<-which(regulon.filtered > r1)
      res <- t(sapply(
        r.filter, function(r2) {tna.hyper(listOfRegulons[r1], universe, listOfRegulons[r2])}
        ))
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
    results <- matrix(NA, nrow=0, ncol=7)
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
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = F)
  ex <- (n/N)*m
  if(m == 0 | n == 0) HGTresults <- 1
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe.Size", "R1.Size", "R2.Size", 
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(hyp.vec)
}

#--------------------------------------------------------------------
#The next 2 functions were extracted from the old version of HTSanalyzeR
#due to a compatibility issues.
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
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = F)
  ex <- (n/N)*m
  if(m == 0) HGTresults <- NA
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe Size", "Gene Set Size", "Total Hits", 
                      "Expected Hits", "Observed Hits", "Pvalue")
  return(hyp.vec)
}

#--------------------------------------------------------------------
#The next 3 functions have being extracted from HTSanalyzeR
#due to compatibility issues!
#--------------------------------------------------------------------
##This function computes enrichment scores for GSEA, running score and 
##position of hits for a gene set.
gseaScores4RTN <- function(geneList, geneSet, exponent=1, mode="score") {
  ##geneSet can either be a character vector with gene ids, or a named integer
  ##vector with additional information (e.g. pos/neg targets)
  if( is.character(geneSet) || is.null(names(geneSet)) ){
    geneSetType<-rep(1,length(geneSet))
    names(geneSetType)<-geneSet
  } else {
    geneSetType<-geneSet
    geneSet<-names(geneSet)
  }
  ##the geneSet should be a subset of the gene universe, i.e. we keep 
  ##only those element of the gene set that appear in the geneList		
  geneSet<-intersect(names(geneList), geneSet)
  ##compute the size of the gene set and of the genelist	
  nh <- length(geneSet)
  N <- length(geneList)
  ##initialize the ES, runningES and the Phit and Pmiss by position 
  ##(the actual values of Phit and Pmiss are actually cumulative sums 
  ##of these 'by position' values)	
  ES <- 0
  Phit <- rep(0, N)
  Pmiss <- rep(0, N)
  runningES <- rep(0, N)
  hits <- rep(FALSE, N)
  hitsType <- rep(0, N)
  ##Stop if the geneSet is larger than the gene universe	
  if(nh > N) {
    stop("NOTE: gene set is larger than gene list!")
  } else {
    ##compute the positions of the hits in the geneList (0 if there 
    ##is no match, 1 if there is a match)	
    hits[which(!is.na(match(names(geneList), geneSet)))] <- TRUE
    ##same as hits, but including mode information is available
    hitsType[which(!is.na(match(names(geneList), names(geneSetType[geneSetType>0]))))] <-  1
    hitsType[which(!is.na(match(names(geneList), names(geneSetType[geneSetType<0]))))] <- -1
    ##if sum(hits)=0 then there is no match between geneList and 
    ##geneSet, and all scores stay at 0.		
    if(sum(hits)!=0){
      ##fill the Phit by position		
      Phit[which(hits)]<-abs(geneList[which(hits)])^exponent
      NR=sum(Phit)
      ##fill the Pmiss by positions			
      Pmiss[which(!hits)]<-1/(N-nh)
      ##do the cumulative sums	and compute the runningES		
      Phit=cumsum(Phit/NR)
      Pmiss=cumsum(Pmiss)
      runningES<-Phit-Pmiss
      ##compute the maximal (positive) and minimal (or maximal 
      ##negative) values of the ES, and choose which one is kept			
      ESmax<-max(runningES)
      ESmin<-min(runningES)
      ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
    }
  }
  ##return the relevant information according to mode  	
  if(mode=="score"){
    return(ES)
  } else if(mode=="graph"){
    return(list("enrichmentScore"=ES, "runningScore"=runningES, "positions"=hitsType))
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
  scores <- parSapply(getOption("cluster"), 1:length(collectionOfGeneSets), function(i) {
    gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm, geneSet = collectionOfGeneSets[[i]], 
                         exponent = exponent, nPermutations = nPermutations)
  }
  )
  return(scores)
}
