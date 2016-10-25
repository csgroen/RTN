
#######################################################
#Functions for partial correlation analysis
#######################################################

fdrtool2format <- function(fdrout, parcor) {
    
    TFindex <- unique(fdrout$index_TF)
    a <- fdrout
    
    Ngenes <- length(rownames(parcor))
    NTFs <- length(colnames(parcor))
    TFtest <- colnames(parcor)
    NAMES <- rownames(parcor)
    m <- matrix(NA, ncol=NTFs, nrow=Ngenes)
    colnames(m) <- TFtest; rownames(m) <- NAMES
    pval <- qval <- prob <- m
    
    for (i in 1:length(TFtest)) {
        TFname <- TFtest[i]
        TFres1 <- fdrout[fdrout[,"TF"] == TFtest[i],]
        TFres1 <- data.frame("target"=as.character(TFres1[,"target"]),"pval"=TFres1[,"pval"],"qval"=TFres1[,"qval"],"prob"=TFres1[,"prob"])
        TFres2 <- fdrout[fdrout[,"target"] == TFtest[i],]
        TFres2 <- data.frame("target"=as.character(TFres2[,"TF"]),"pval"=TFres2[,"pval"],"qval"=TFres2[,"qval"],"prob"=TFres2[,"prob"])
        res <- data.frame(rbind(TFres1,TFres2))
        res$target <- as.character(res$target)
        res <- rbind(res, c(TFname, 0,0,1))
        
        pval[,i] <-  as.numeric(res[match(NAMES, res[,"target"]),"pval"])
        qval[,i] <-  as.numeric(res[match(NAMES, res[,"target"]),"qval"])
        prob[,i] <-  as.numeric(res[match(NAMES, res[,"target"]),"prob"])
    
    }
    
    return(list("tn"=parcor, "qval"=qval, "prob"=prob, "pval"=pval))
}

takevector <- function(pmim) {
    
    x <- meltX(pmim)
    
    #remove duplicated entries
    #(n^2 - n)/2
    MIN <- pmin(x$index_target, x$index_TF)
    MAX <- pmax(x$index_target, x$index_TF)
    indexes_order <- paste(MIN,MAX,sep=".")
    x <- x [ !duplicated(indexes_order),]
    
    #remove diagonal (n)
    if (sum(x[,"value"] == 999) != ncol(pmim)) {stop("there are more values in the diagonal than expected")}
    x <- x[x[,"value"] != 999,]
    
    #check size
    sizeX <- nrow(x)
    expected <- ncol(pmim)*nrow(pmim) - ( (ncol(pmim)*ncol(pmim) - ncol(pmim)) / 2 ) - ncol(pmim)
    if(sizeX != expected){stop("inconsistent with expected")}
    
    return(x) 
}

#based on analytical solution
computePvaluesPcor <-function(pmim, verbose=TRUE){
    if (verbose) {cat("-taking vector\n")}
    parcor_vect <- takevector(pmim) #takes the vector of TF.target without any duplicated entries (TF.target == target.TF)
    if (verbose) {cat("-applying fdrtool\n")}
    fdr.out <- fdrtool(parcor_vect[,"value"],statistic="correlation", plot=F, verbose=F)
    pval = fdr.out$pval
    qval = fdr.out$qval
    prob = 1 - fdr.out$lfdr #posterior probabilities are 1 - local FDR
    fdrout = cbind(parcor_vect, pval, qval, prob)
    res <- fdrtool2format(fdrout, pmim)
    
    #return results
    return(list(tn=pmim,qval=res[["qval"]], pval=res[["pval"]], prob=res[["prob"]]))  
}

meltX <- function(pmim) {

    #getrid of the diagonal
    diag(pmim[colnames(pmim),]) <- 999
     
    #melt X
    x <- melt(pmim, value.name = "value")
    colnames(x) <- c("target","TF","value")
    x$target <- as.character(x$target)
    x$TF <- as.character(x$TF)
    x$value <- as.numeric(as.character(x$value))
    index_target <- match(x$target, rownames(pmim))
    index_TF <- match(x$TF, rownames(pmim))
    
    x <- data.frame(x[,1:3],"index_target"=index_target, "index_TF"=index_TF)
    
    return(x)
}

computeZscores <- function(tfmode, pval, alt="two.sided") {
    #compute the Z score for each network
    if (any(pval==0)) {
    minPVAL <- min(pval[pval !=0])
    if(minPVAL > 2.2204461e-16) {minPVAL <- 2.220446e-16} 
    pval[pval==0] <- minPVAL #avoid Inf for Zscore
    }
    if (any(pval==1)) {
        pval[pval==1] <- 0.9999999 #avoid -Inf for Zscore
    }
    
    if (alt=="two.sided") {
        z <- qnorm(pval/2, lower.tail=F) * tfmode
    }
    
    if (alt=="one.sided") {
        z <- qnorm(pval,lower.tail=F) * tfmode
    }

    return(z)
}



applycorpcor <- function(object, verbose=TRUE) {
    gexp= t(object@gexp)
    if (verbose) {cat("-applying pcor.shrink\n")}
    pmim <- pcor.shrink(gexp, verbose=F)
    tfs <- object@transcriptionFactors
    pmim <- pmim[,match(tfs,colnames(pmim))]
    net <- computePvaluesPcor(pmim,  verbose=verbose)
    net$tfmode <- sign(net$tn)
    net$z.score <- computeZscores(net$tfmode, net$pval, alt="two.sided")
    object@results.pcor <- net
    return(object)
}

##############################################################################
#Network Fusion
##############################################################################

combine.pv <- function(p, w, sig, alt="two.sided") {
    
    if (alt == "two.sided") {
        z <- computeZscores(sig, p, alt="two.sided")
        z.combined <- sum(z*w) / sqrt(sum(w^2))
        p.combined <- 2*pnorm(abs(z.combined), lower.tail=F)
        return(list("p"=p.combined, "z"=z.combined))
    }
    if (alt == "one.sided") {
        z <- computeZscores(sig, p, alt="one.sided")
        z.combined <- sum(z*w) / sqrt(sum(w^2))
        p.combined <- pnorm(z.combined, lower.tail=F)
        return(list("p"=p.combined, "z"=z.combined))
    }
}

combineZscores <- function(w, z) {
    z.combined <- sum(z*w) / sqrt(sum(w^2))
    p.combined <- 2*pnorm(abs(z.combined), lower.tail=F)
    return(list("p"=p.combined, "z"=z.combined))
}


fuseNetwork <- function(objectlist) {

    #get lists
    sampleSize <- c()
    pvalNets <- list()
    parcorNets <- list()
    probNets <- list()
    zNets <- list()
    for (object in objectlist) { 
        net <- object@results.pcor
        gexp <- object@gexp
        sampleSize <- c(sampleSize, ncol(gexp)) 
        pvalNets <- c(pvalNets, net['pval'])
        parcorNets <- c(parcorNets, net['tn']) 
        probNets <- c(probNets, net['prob'])
        zNets <- c(zNets, net['z.score'])
    }
    names(sampleSize) <- names(pvalNets) <- names(parcorNets) <- names(probNets) <- names(zNets) <- names(objectlist)
    
    #intersect all common TFs and genes in all networks
    TFs <- Reduce(intersect,lapply(pvalNets, colnames))
    genes <- Reduce(intersect,lapply(pvalNets, rownames))
    
    #combine Pvalues
    emptydframe <- matrix(NA, ncol=length(TFs), nrow=length(genes))
    rownames(emptydframe) <- genes
    colnames(emptydframe) <- TFs
    
    tfmode <- emptydframe
    combinedZ <- emptydframe
    combinedPvals <- emptydframe
    
    #will create combinedZ; tfmode; combinedPvals
    cat("-combining p-values")
    pb <- txtProgressBar(min=0,max=length(TFs),style=3)
    for (i in 1:length(TFs)) {
        tf <- TFs[i]
        z <- lapply(zNets, function (x) {x[match(genes, rownames(x)),match(tf,colnames(x))]})
        z <- do.call("cbind", z)
        w <- sampleSize[colnames(z)]
        comb <- lapply(1:nrow(z), function(j) { combineZscores(w, z[j,])})
        comb <- do.call("rbind", comb)
        #comb <- apply(ps, 1, function(p) combine.test(p=p,weight=w, method="z.transform"))
        combinedPvals[,tf] <- unlist(comb[,"p"])
        z.score <- unlist(comb[,"z"])
        combinedZ[,tf] <- z.score
        if (any(is.infinite(z.score))) {stop('infinite Z scores, unable to determine sign')}
        tfmode[,tf] <- sign(z.score)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    
    #perform pvalue adjustment
    #From combinedPvals will calculate qval and prob
    print("-performing pvalue adjustment with fdrtool")
    pcor_vect <- takevector(combinedPvals) #takes the vector of TF.target without any duplicated entries (TF.target == target.TF)
    fdr.out <- fdrtool(pcor_vect[,"value"],statistic="pvalue", plot=T, verbose=F)
    pval = fdr.out$pval
    qval = fdr.out$qval
    prob = 1 - fdr.out$lfdr #posterior probabilities is 1 - local FDR
    fdrout = cbind(pcor_vect, pval, qval, prob)
    res <- fdrtool2format(fdrout, combinedPvals)
    combinedQvals <- res[["qval"]]
    combinedProb <- res[["prob"]]            
    
    #make sure is the same order
    tfmode <- tfmode[match(genes, rownames(tfmode)), match(TFs, colnames(tfmode))]
    combinedZ <- combinedZ[match(genes, rownames(combinedZ)), match(TFs, colnames(combinedZ))]
    combinedPvals <- combinedPvals[match(genes, rownames(combinedPvals)), match(TFs, colnames(combinedPvals))]
    combinedQvals <- combinedQvals[match(genes, rownames(combinedQvals)), match(TFs, colnames(combinedQvals))]
    combinedProb <- combinedProb[match(genes, rownames(combinedProb)), match(TFs, colnames(combinedProb))]
    
    net <- list("fusedNet"=list("qval"=combinedQvals, "pval"=combinedPvals, "prob"=combinedProb, "tfmode"=tfmode, "z.score"=combinedZ))
    return(net)
}



##############################################################################
#Convert 'tni' object to 'Regulon' (to run VIPER)
##############################################################################

TFuse2regulon <- function(object, pvalueCutoff= 0.05, filterTFs=NULL, returnclassregulon=T) {
    
    if (class(object) == "TNIpcor") {
        net <- object@results.pcor
        tfmode <- net$tfmode
        prob <- net$prob
        qval <- net$qval
    }
    if (class(object) == "list")    {
        net <- object
        if (!("tfmode" %in% names(net))) {stop("tfmode not found")}
        if (!("prob" %in% names(net))) {stop("prob not found")}
        if (!("qval" %in% names(net))) {stop("qval not found")}
        tfmode <- net$tfmode
        prob <- net$prob
        qval <- net$qval
    }
    
    if (is.null(filterTFs)) {TFs <- colnames(tfmode)}
    if (!is.null(filterTFs)) {TFs <- filterTFs; TFs <- intersect(TFs, colnames(tfmode))}
    scores <- prob[,TFs]
    mode1 <- tfmode[,TFs]
    
    aracne <- list()
    for (tf in TFs) {
        reg <- qval[,tf]
        which.reg <- reg <= pvalueCutoff
        likelihood <- scores[which.reg,match(tf, colnames(scores))]
        tfmode <- mode1[which.reg,match(tf, colnames(mode1))]
        aracne[[tf]] <- list("tfmode"=tfmode, "likelihood"=likelihood)
    }
    
    # removing missing data from the aracne regulon
    aracne <- aracne[names(aracne) != "NA"]
    aracne <- lapply(aracne, function(x) {
        filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
        x$tfmode <- x$tfmode[filtro]
        x$likelihood <- x$likelihood[filtro]
        return(x)
    })
    regul <- aracne[sapply(aracne, function(x) length(names(x$tfmode)))>0]
    if (returnclassregulon) class(regul) <- "regulon"
    return(regul)
}
