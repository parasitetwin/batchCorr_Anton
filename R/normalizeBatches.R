#' Between-batch normalisation
#'
#' Batches are feature-wise normalised by Ref samples if passing heuristic criteria (CV and fold change).
#' Otherwise normalized by population median (still per feature).
#'
#' @param peakTableCorr Multi-batch drift-corrected peak table (samples in rows; features in columns)
#' @param peakTableOrg Multi-batch original peak table (samples in rows; features in columns)
#' @param batches Vector (length=nrow(PeakTab)) of batch identifiers
#' @param sampleGroup Vector (length=nrow(PeakTab)) of sample type (e.g. "sample", "QC", "Ref)
#' @param refGroup Identifier of reference samples in sampleGroups (defaults to "QC")
#' @param population Identifier of population samples in sampleGroups (all (default) or any type of samples present in sampleGroups)
#' @param CVlimit  CV criterion to pass for Ref samples per batch (defaults to 0.3)
#' @param FCLimit Fold-change criterion for intensity (in relation to average intensity FC between batches; defaults to 5)
#' @param medianZero Strategy for substituting median value for population normalization when median is zero (-> Inf). Either 'mean' or 'min' (i.e. lowest non-zero value; default)
#'
#' @return An list object containing:
#' @return peakTable: Normalised multi-batch peak table
#' @return refCorrected: Boolean matrix with info on which batches were normalised by reference samples; others were normalized by population medians
#' @examples
#' # Note that the example data does not include any biological samples, in
#' # which case population = "sample" does not work
#' data("ThreeBatchData")
#' normData <- normalizeBatches(peakTableCorr = PTfill, batches = meta$batch,
#'                              sampleGroup = meta$grp, refGroup = "Ref",
#'                              population = "all")
#' @export
normalizeBatches <- function(peakTableCorr, peakTableOrg, batches, sampleGroup, refGroup='QC', population='all', CVlimit=0.3, FCLimit = 5, medianZero = c('mean', 'min')){
  if (missing(peakTableOrg)) peakTableOrg <- peakTableCorr
  if (!(identical(dim(peakTableCorr), dim(peakTableOrg)) & identical(colnames(peakTableCorr), colnames(peakTableOrg)))) stop('Mismatch between peakTableCorr and peakTableOrg')
  if (missing(medianZero)) medianZero <- 'min'
  if (population=='all') popSample <- rep(TRUE,nrow(peakTableCorr)) else {
    if (!population%in%sampleGroup) (stop('population identifier needs to be present in sampleGroups\nConsider setting population="all"')) else {
      popSample <- sampleGroup==population
    }
  }
  # Extract info and declare/allocate variables
  peakTableNormalized <- peakTableCorr
  uniqBatch <- unique(batches)
  nBatch <- length(uniqBatch)
  nFeat <- ncol(peakTableCorr)
  CVMatrix=matrix(nrow=nBatch,ncol=nFeat) # Ref sample CVs per batch (row) and feature (col)
  rownames(CVMatrix) <- uniqBatch
  colnames(CVMatrix) <- colnames(peakTableCorr)
  RefMeanIntensity <- CVMatrix # Ref sample mean intensity  per batch (row) and feature (col)
  RefNormMatrix <- matrix(FALSE, nrow = nBatch, ncol = nFeat) # Boolean (flag) for which batches (rows) of which features (columns) are normalized by reference samples
  MeanIntensityRatios <- matrix(1, nrow = nBatch, ncol = nBatch) # Mean intensity ratios between batches (all features)
  # Aggregate CV and average intensities for the reference sample type per batch
  for (b in 1:nBatch) {
    batch <- uniqBatch[b]
    peakTableBatch <- peakTableCorr[batches==batch & sampleGroup==refGroup,]
    CVMatrix[b,]=ifelse(cv(peakTableBatch)<=CVlimit,TRUE,FALSE)  # Criterion for QC feature CV
    RefMeanIntensity[b,]=apply(peakTableBatch,2,mean)
  }
  # Find average intensity ratios between batches
  for (b in 1:(nBatch - 1)) {
    for (bb in (b + 1):nBatch) {
      MeanIntensityRatios[bb, b] = mean(RefMeanIntensity[bb, ])/mean(RefMeanIntensity[b,])
      MeanIntensityRatios[b, bb] = 1/MeanIntensityRatios[bb, b]
    }
  }
  # Perform normalisation per feature
  for (feat in 1:nFeat) {
    # Find intensity ratio per feature between batches
    featureIntensityRatios = matrix(1, nrow = nBatch, ncol = nBatch)
    for (b in 1:(nBatch - 1)) {
      for (bb in (b + 1):nBatch) {
        featureIntensityRatios[bb, b] = RefMeanIntensity[bb, feat]/RefMeanIntensity[b,feat]
        featureIntensityRatios[b, bb] = 1/featureIntensityRatios[bb, b]
      }
    }
    # Identify candidates for ref normalization
    candidates = abs(log(featureIntensityRatios/MeanIntensityRatios)) <= log(FCLimit)  # Criterion for intensity ratio
    candidates[is.na(candidates)]=FALSE # Convert missing values -> FALSE (i.e. not a candidate)
    for (b in 1:nBatch) { # Criterion for CV < limit
      if (CVMatrix[b, feat] == FALSE | is.na(CVMatrix[b, feat])) {
        candidates[, b] = candidates[b, ] = FALSE
      }
    }
    # Perform ref normalization
    refBatch = min(which(colSums(candidates) == max(colSums(candidates)))) # Find reference batch
    refCorrFlags = candidates[, refBatch] # Extract flags for reference sample correction
    refCorrIndex = which(refCorrFlags) # Find which batches to normalise to the reference
    refCorrIndex=refCorrIndex[refCorrIndex!=refBatch]
    refIntensity = RefMeanIntensity[refBatch, feat] # Ref batch intensity
    for (b in refCorrIndex) { # Correct batches to reference batch intensity
      correctionFactor = refIntensity/RefMeanIntensity[b, feat]
      peakTableNormalized[batches == uniqBatch[b], feat] = peakTableCorr[batches == uniqBatch[b], feat] * correctionFactor
    }
    RefNormMatrix[, feat] = refCorrFlags
    # Perform population (median) normalisation for the other batches
    if (length(refCorrIndex)+1!=nBatch) {
      refCorrected=peakTableNormalized[batches%in%uniqBatch[c(refBatch,refCorrIndex)] & popSample, feat] # All samples of "population" corrected by reference samples
      refCorrMedian=median(refCorrected) # Extract their median value
      WhichPOPCorr=which(!refCorrFlags) # Which batches to correct by population median instead
      for (n in WhichPOPCorr) {
        populationMedian <- median(peakTableOrg[batches == uniqBatch[n] & popSample, feat])
        if (populationMedian==0) {
          if (medianZero=='min') {
            populationMedian <- peakTableOrg[batches == uniqBatch[n] & popSample, feat]
            populationMedian <- min(populationMedian[populationMedian!=0])
          } else if (medianZero=='mean') {
            populationMedian <- mean(peakTableOrg[batches == uniqBatch[n] & popSample, feat])
          } else stop('Other options not included at present.')
        }
        correctionFactor=refCorrMedian/populationMedian # Calculate correction factor for "population" samples
        peakTableNormalized[batches == uniqBatch[n], feat] = peakTableOrg[batches == uniqBatch[n], feat] * correctionFactor
      }
    }
  }
  return(list(peakTable = peakTableNormalized, refCorrected = RefNormMatrix))
}
