#' Merge batches after drift correction
#'
#' The output of the within-batch drift correction is a correction object.
#' This function merges peak tables from several batches by extracting information from the correction objects.
#' The user must specify a minimum proportion of qualified batches per feature, i.e. such batches where the QC CV is < the specified limit.
#' There is thus a risk that features with poor quality (in certain batches) are present, but the features are present in high quality in sufficient proportion of batches to anyway warrant inclusion.
#' @param batchList A list of correction objects (after drift correction)
#' @param qualRatio Proportion of
#'
#' @return A list object containing
#' @return `peakTable` A merged peak table
#' @return `batch` Batch identifier (per sample)
#' @return `injection` Injection number (per sample)
#' @export
#'
#' @examples
#' batchList <- list(b1=b1Corr, b2=b2Corr, b3=b3Corr)
#' mergeList <- mergeBatches(batchList)
#' cbind(mergeList$batch,mergeList$injection)
#' dim(mergeList$peakTable)
mergeBatches <- function(batchList, qualRatio=0.5) {
  nBatch <- length(batchList)
  nQual <- ceiling(qualRatio*nBatch)
  if(is.null(names(batchList))) {
    batchNames <- 1:nBatch
  } else {
    batchNames <- names(batchList)
  }
  nSamp <- numeric(nBatch)
  injections <- list()
  qualFeatures <- list()
  peakTables <- list()
  # Extract relevant data from corr objects
  for (batch in 1:nBatch) {
    nSamp[batch] <- length(batchList[[batch]]$TestInjs)
    injections[[batch]] <- batchList[[batch]]$TestInjs
    qualFeatures[[batch]] <- batchList[[batch]]$finalVars
    peakTables[[batch]] <- batchList[[batch]]$TestFeatsCorr
  }
  # Aggregate data
  injections <- do.call(c,injections)
  batches <- rep(batchNames,nSamp)
  qualFeatures <- do.call(c,qualFeatures)
  qualFeatures <- names(which(table(qualFeatures)>=nQual))
  peakTables <- do.call(rbind,peakTables)
  peakTables <- peakTables[,colnames(peakTables)%in%qualFeatures]
  return(list(peakTable=peakTables, batch=batches, injection=injections))
}
