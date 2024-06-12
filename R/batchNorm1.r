#' Extract m/z and rt from peak table
#'
#' Extract features from peak table and report their m/z and rt values
#'
#' @param PT a peak table with variables as columns
#' @param sep character separating mz from rt, e.g. "_"
#' @param start character from which to start the read of peakInfo (from PT colnames)
#' @param timepos Which position carries info about rt (1 for before separator; 2 for after separator)
#' @return a matrix with m/z and rt of features as columns
#'
#' @examples
#' data('ThreeBatchData')
#' # Extract peakinfo (i.e. m/z and rt of features). These column names have 2
#' # leading characters describing LC-MS mode -> start at 3
#' peakIn <- peakInfo(PT = PTnofill, sep = '@', start = 3)
#'
#' @export
peakInfo=function(PT,sep='_',timepos=2,start=1) {
  #Making sure separator is present in colnames of peak table (PT)
  if(!any(grepl(sep, colnames(PT)))){
    message(paste0("Error in 'peakInfo': Separator not present in column names of peak table."))
    stop()
  }
  #Making sure type of start is numeric in nature
  if(!is.numeric(start) || !is.numeric(timepos)){
    whichNotNumeric <- c("'start'", "'timepos'")
    message(paste0("Error in 'peakInfo': ",
                   paste(whichNotNumeric[c(!is.numeric(start),
                                           !is.numeric(timepos))], collapse=" & "),
                   " not numeric."))
    stop()
  }
  #Checking that separator is not longer than 1
  if(length(sep) > 1){
    message(paste0("Error: Only use of one separator possible"))
    stop()
  }
  
  #Making matrix of mz and rt after splitting by separator
  peakInfo=matrix(unlist(strsplit(colnames(PT),sep)),ncol=2,byrow=TRUE)
  peakInfo[,1]=substr(peakInfo[,1],start,max(nchar(peakInfo[,1])))
  peakInfo=matrix(as.numeric(peakInfo),ncol=2)
  #Reverse column order if time is not after separator
  if (timepos!=2) peakInfo=peakInfo[,2:1]
  colnames(peakInfo)=c('mz','rt')
  rownames(peakInfo)=paste('feature',1:nrow(peakInfo),sep='_')
  return(peakInfo)
}

#' Limit variables/features to those common between batches
#'
#' Extract features from peak table and report their m/z and rt values
#' @param batchFeats list with feature names per batch
#' @param batchLimit lower limit of number of batches in which the feature needs to be present to be included in the final list (defaults to length of `batchFeats` list)
#' @return a vector with the features present in at least `batchLimit` batches
#' @noRd
featComb=function(batchFeats,batchLimit) {
  nBatch=length(batchFeats)
  if (missing(batchLimit)) batchLimit=nBatch
  vars=batchFeats[[1]]
  for (b in 2:nBatch) {
    vars=c(vars,batchFeats[[b]])
  }
  varTab=table(vars) #Tabularize number of occurrences of features
  finalVars=t(as.data.frame(varTab[varTab>=batchLimit])) # Bring out the features present >= batchLimit times NB! Not correctly sorted!
  featInfo=peakInfo(finalVars) # Bring out mz and rt to sort them properly
  finalVars=colnames(finalVars)  # take names of features only
  finalVars=finalVars[order(featInfo[,1],featInfo[,2])]  # Sort them according to mz and rt
  return(finalVars)
}


#' Extract features from multiple batch data (deprecated)
#'
#' Extract features present in all batches and combine them into a master peaktable.
#' @param batchObjs a list with batch objects (from within batch drift correction)
#' @param batchLimit lower limit of number of batches in which the feature needs to be present to be included in the final list (defaults to length of `batchFeats` list)
#' @param finalVars a vector with the features to bring out from batch PTs. NB! All features need be present in all batches (defaults to all features present in all batches)
#' @return a combined, but NOT normalised, peaktable - limited to common features
#' @noRd
batchComb=function(batchObjs,batchLimit,finalFeats) {
  warning('Still works but deprecated: Use mergeBatches() instead')
  nBatch=length(batchObjs)
  if (missing(batchLimit)) batchLimit=nBatch
  if (missing(finalFeats)) {
    batchFeats=list()
    for (b in 1:nBatch) {
      batchFeats[[b]]=batchObjs[[b]]$finalVars
    }
    finalFeats=featComb(batchFeats,batchLimit)
  }
  PTComb=subset(batchObjs[[1]]$TestFeatsFinal,select=finalFeats)
  for (b in 2:nBatch) {
    PTComb=rbind(PTComb,subset(batchObjs[[b]]$TestFeatsFinal,select=finalFeats))
  }
  return(PTComb)
}

#' BN: Info on reference samples aggregated on batch level (Deprecated)
#'
#' Reference samples are aggregated on batch level
#' @param PT a multi-batch master peak table
#' @param meta metadata with batch (col1) and sample type (col2)
#' @param grpType a sample type identifier for reference samples
#' @param CVlimit CV criterion to pass for Ref samples per batch
#' @return an object containing:
#' @return CV: boolean per batch & feature in CV<limit
#' @return aveInt: average reference intensity per batch & feature
#' @noRd
refOut=function(PT,meta,grpType='R',CVlimit=0.3) {
  warning('Still works but deprecated: Functionality merged with refCorr')
  batch=meta[,1]
  grp=meta[,2]
  uniqBatch=unique(batch)
  CVMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(CVMat)=uniqBatch
  aveIntMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(aveIntMat)=uniqBatch
  for (b in 1:length(uniqBatch)) {
    bat=uniqBatch[b]
    PTbatch=PT[batch==bat & grp==grpType,]
    CVMat[b,]=ifelse(cv(PTbatch)<=CVlimit,TRUE,FALSE)
    aveIntMat[b,]=apply(PTbatch,2,mean)
  }
  return(list(CV=CVMat,aveInt=aveIntMat))
}

#' BN: Between batch normalisation by Ref samples
#'
#' Batches are normalised by Ref samples if passing heuristic distance criterion.
#'
#' @param PT a multi-batch master peak table
#' @param batch Batch identifier
#' @param grp Sample type identifier (e.g. 'Sample', 'QC', 'Ref' or similar)
#' @param grpType Reference sample identifier (e.g. 'Ref')
#' @param CVlimit  CV criterion to pass for Ref samples per batch
#' @param FCLimit Fold-change criterion for intensity (in relation to average intensity FC between batches)
#'
#' @return an object containing:
#' @return PTRef: Reference sample-normalised multi-batch peak table
#' @return refCorr: Boolean matrix with info on which batches were normalised by reference samples
#' @return CV: boolean per batch & feature in CV<limit
#' @return aveInt: average reference intensity per batch & feature
#' @return PTOrg (indata peaktable)
#' @return batch Batch identifier
#' @return grp Sample type identifier
#' @return grpType Reference sample identifier
#' @importFrom stats median
#' @noRd
refCorr=function (PT, batch, grp, grpType='R', CVlimit=0.3, FCLimit = 5){
  # Extract info and declare/allocate variables
  PTcorr = PT
  uniqBatch=unique(batch)
  nBatch=length(uniqBatch)
  nFeat = ncol(PT)
  CVMat=matrix(nrow=nBatch,ncol=nFeat)
  rownames(CVMat)=uniqBatch
  aveIntMat=CVMat
  refCorrMat = matrix(FALSE, nrow = nBatch, ncol = nFeat)
  meanIntRat = matrix(1, nrow = nBatch, ncol = nBatch)
  # Aggregate CV and average intensities for the reference sample type per batch
  for (b in 1:nBatch) {
    bat=uniqBatch[b]
    PTbatch=PT[batch==bat & grp==grpType,]
    CVMat[b,]=ifelse(cv(PTbatch)<=CVlimit,TRUE,FALSE)  # Criterion for QC feature CV
    aveIntMat[b,]=apply(PTbatch,2,mean)
  }
  # Find average intensity ratios between batches
  for (b in 1:(nBatch - 1)) {
    for (bb in (b + 1):nBatch) {
      meanIntRat[bb, b] = mean(aveIntMat[bb, ])/mean(aveIntMat[b,])
      meanIntRat[b, bb] = 1/meanIntRat[bb, b]
    }
  }
  # Check potential normalisation candidates
  cvFlags = apply(CVMat, 2, function(x) sum(x, na.rm = TRUE))
  whichFeatsCV = which(cvFlags > 1)
  lenCV = length(whichFeatsCV)
  for (feat in 1:nFeat) {
    # feat = whichFeatsCV[lc]
    featIntRat = matrix(1, nrow = nBatch, ncol = nBatch)
    for (b in 1:(nBatch - 1)) { # Find intensity ratio per feature between batches
      for (bb in (b + 1):nBatch) {
        featIntRat[bb, b] = aveIntMat[bb, feat]/aveIntMat[b,feat]
        featIntRat[b, bb] = 1/featIntRat[bb, b]
      }
    }
    featFlags = abs(log(featIntRat/meanIntRat)) <= log(FCLimit)  # Criterion for intensity ratio
    # Cleanup
    featFlags[is.na(featFlags)]=FALSE
    for (b in 1:nBatch) {
      if (CVMat[b, feat] == FALSE | is.na(CVMat[b, feat])) {
        featFlags[, b] = featFlags[b, ] = FALSE
      }
    }
    refBatch = min(which(colSums(featFlags) == max(colSums(featFlags)))) # Find reference batch
    refCorr = featFlags[, refBatch] # Extract flags for reference sample correction
    WhichRefCorr = which(refCorr) # Find which batches to normalise to the reference
    refInt = aveIntMat[refBatch, feat] # Ref batch intensity
    refCorrIndex=WhichRefCorr[WhichRefCorr!=refBatch]
    for (b in refCorrIndex) { # Correct batches to reference batch intensity
      corrFact = refInt/aveIntMat[b, feat]
      PTcorr[batch == uniqBatch[b], feat] = PT[batch == uniqBatch[b], feat] * corrFact
    }
    # Perfrom population (median) normalisation for the other batches
    featCorr=PTcorr[batch %in% uniqBatch[WhichRefCorr], feat] # All samples corrected by reference samples ->
    Median=median(featCorr) # Extract their median value
    WhichPOPCorr=which(!refCorr) # Which batches to correct by population median instead
    for (n in WhichPOPCorr) {
      corrFact=Median/median(PTcorr[batch == uniqBatch[n], feat])
      PTcorr[batch == uniqBatch[n], feat] = PTcorr[batch == uniqBatch[n], feat] * corrFact
    }
    refCorrMat[, feat] = refCorr
  }
  return(list(PTRef = PTcorr, refCorr = refCorrMat, CV = CVMat, aveInt = aveIntMat, PTOrg = PT, batch = batch, grp = grp, grpType = grpType))
}

