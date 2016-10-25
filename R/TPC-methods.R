##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "TPC",
          function(.Object, tniclassObj) {
          .Object@gexp= tniclassObj@gexp
          .Object@transcriptionFactors= tniclassObj@transcriptionFactors
          .Object@results.pcor= tniclassObj@results
          status = tniclassObj@status
          status = c(status["Preprocess"], "corpcor"="[ ]")
          .Object@status <- status
          .Object
          }
)

#methods
setMethod(
  "tni2tpc.corpcor",
  "TNI",
  function(object, verbose=TRUE){
    object <- new("TPC", tniclassObj=object)
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input 'object' needs preprocessing")
    object <- applycorpcor(object)
    object@status["corpcor"] <- "[x]"
    if(verbose)cat("-Partial correlation analysis complete! \n\n")
    return(object)
  }
)

setMethod(
  "tpc.fuseNetwork",
  "TPC",
  function(objectlist, verbose=TRUE){
    #check if class in objectslist is OK
    fused <- fuseNetwork(objectlist)
    if(verbose)cat("-Network fusion complete! \n\n")
    return(fused)
  }
)

##count significant edges
setMethod(
  "tpc.countedges",
  "TPC",
  function(object, cutoff.method="prob", cutoff.tr = 0.9, verbose=TRUE) {
    
    if (class(object) == "TPC") {net <- object@results.pcor}else{net <- object }
    tfmode <- net$tfmode
    total_edges <- nrow(tfmode) * ncol(tfmode)
    
    if (cutoff.method == "prob") {
      sig_edges <- sum(net$prob > cutoff.tr, na.rm=T)
      perc <- sig_edges * 100 / total_edges
    }
    
    if (cutoff.method == "qval") {
      sig_edges <- sum(net$qval <= cutoff.tr, na.rm=T)
      perc <- sig_edges * 100 / total_edges
    }
    
    if (verbose) {
      cat(paste0("Significant edges: ",sig_edges, " Corresponding to ",round(perc,2), "% of all possible edges\n"))
    }
    return(sig_edges)
  }
  
)
