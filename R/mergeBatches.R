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
#' @return `peakTableOrg` A merged peak table of original data
#' @return `peakTableCorr` A merged peak table of drift-corrected data
#' @return `batch` Batch identifier (per sample)
#' @return `injection` Injection number (per sample)
#'
#' @examples \donttest{
#' data('ThreeBatchData') 
#' set.seed(2024)
#' # Get batches
#' batchB <- getBatch(peakTable = PTfill, meta = meta, 
#'                    batch = meta$batch, select = 'B')
#' batchF <- getBatch(peakTable = PTfill, meta = meta, 
#'                    batch = meta$batch, select = 'F')
#' # Drift correction using QCs
#' BCorr <- correctDrift(peakTable = batchB$peakTable, 
#'                       injections = batchB$meta$inj, 
#'                       sampleGroups = batchB$meta$grp, QCID = 'QC', 
#'                       G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'))
#' # More unbiased drift correction using QCs & external reference samples
#' FCorr <- correctDrift(peakTable = batchF$peakTable, 
#'                       injections = batchF$meta$inj,
#'                       sampleGroups = batchF$meta$grp, QCID = 'QC',
#'                       RefID='Ref', G = seq(5,35,by=3), 
#'                       modelNames = c('VVE', 'VEE'))
#' # Merge batches for batch normalization, for example
#' mergedData <- mergeBatches(list(BCorr, FCorr))
#' }
#'
#' @export
mergeBatches <- function(batchList, qualRatio=0.5) {
  # Determining number of batches and batch-ratio needed to surpass to keep a feature
  nBatch <- length(batchList)
  nQual <- ceiling(qualRatio*nBatch)
  #If no batchnames supplied in batchList names batches from 1-nBatches
  if(is.null(names(batchList))) {
    batchNames <- 1:nBatch
  } else {
    batchNames <- names(batchList)
  }
  # Settings up lists and vectors to store batch-specific information in
  nSamp <- numeric(nBatch)
  injections <- list()
  qualFeatures <- list()
  peakTablesOrg <- peakTablesCorr <- list()
  # Extract relevant data from corr objects
  for (batch in 1:nBatch) {
    nSamp[batch] <- length(batchList[[batch]]$TestInjs)
    injections[[batch]] <- batchList[[batch]]$TestInjs
    qualFeatures[[batch]] <- batchList[[batch]]$finalVars
    peakTablesOrg[[batch]] <- batchList[[batch]]$TestFeats
    peakTablesCorr[[batch]] <- batchList[[batch]]$TestFeatsCorr
  }
  # Aggregate data and merge
  injections <- do.call(c,injections)
  batches <- rep(batchNames,nSamp)
  qualFeatures <- do.call(c,qualFeatures)
  qualFeatures <- names(which(table(qualFeatures)>=nQual))
  peakTablesOrg <- do.call(rbind,peakTablesOrg)
  peakTablesOrg <- peakTablesOrg[,colnames(peakTablesOrg)%in%qualFeatures]
  peakTablesCorr <- do.call(rbind,peakTablesCorr)
  peakTablesCorr <- peakTablesCorr[,colnames(peakTablesCorr)%in%qualFeatures]
  return(list(peakTableOrg=peakTablesOrg, peakTableCorr=peakTablesCorr, batch=batches, injection=injections))
}
