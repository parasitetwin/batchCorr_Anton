#' Extract specific batch from peaktable and metadata
#'
#' @param peakTable Multi-batch peak table
#' @param meta Multi-batch metadata (including e.g. batch, injection sequence, sample tye, sample name, ...)
#' @param batch vector (length = nSamples) containing batch information (e.g. A, B, C)
#' @param select which batch to extract
#'
#' @return list object
#' @return $peakTable: Single batch peak table
#' @return $meta: Single batch metadata
#'
#' @examples
#' data('ThreeBatchData') 
#' # Get batches for drift correction
#' batchB <- getBatch(peakTable = PTfill, meta = meta, 
#'                    batch = meta$batch, select = 'B')
#' batchF <- getBatch(peakTable = PTfill, meta = meta, 
#'                    batch = meta$batch, select = 'F')
#'
#' @export
getBatch=function(peakTable,meta,batch,select) {
  #Checking if nrow of peakTable corresponds to nrow/length of meta
  if(is.null(dim(meta))){
    if(length(meta) != nrow(peakTable)){
      # message(paste0("Warning in 'getBatch': Meta length not equivalent to rows in peakTable"))
      stop("Meta length not equivalent to rows in peakTable")
    }
  } else {
    if(nrow(meta) != nrow(peakTable)){
      # message(paste0("Warning in 'getBatch': Meta length not equivalent to rows in peakTable"))
      stop("Meta length not equivalent to rows in peakTable")
    }
  }
  
  #Checking that select is not a vector and throws error if it is
  if(length(select) > 1){
    # message(paste0("Error in 'getBatch': Vector of batches to select is not supported"))
    stop("Number of batches to select > 1")
  }
  
  peakTable=peakTable[batch==select,]
  
  #Checking if nrow in peakTable is larger than 0 and throwing error if not
  if(nrow(peakTable) == 0){
    # message("Error in 'getBatch': No such batch in batch parameter")
    stop("No such batch in batch vector")
  }
  
  #Checking if meta is vector and in that case turning it into a matrix
  if(is.null(dim(meta))){
    meta=matrix(meta,ncol = 1)
  }
  meta=meta[batch==select,]
  
  return(list(peakTable=peakTable,meta=meta))
}

#' Extract specific sample group from peaktable and metadata
#'
#' Matches pattern in `select` to sample group and extracts
#' @param peakTable Multi-batch peak table
#' @param meta Multi-batch metadata (including e.g. batch, injection sequence, sample tye, sample name, ...)
#' @param sampleGroup vector (length = nSamples) containing sample group information (e.g. QC, Sample, Reference)
#' @param select which sample group to extract
#'
#' @return list object
#' @return $peakTable: Sample group peak table
#' @return $meta: Sample group metadata
#' @noRd
#'
#' @examples
#' data(ThreeBatchData)
#' batchB=getBatch(peakTable=PTfill, meta=meta, batch=meta$batch, select='B')
#' batchBQC=.getGroup(peakTable=batchB$peakTable, meta=batchB$meta, sampleGroup=batchB$Meta$grp, select='QC')
#' @noRd
.getGroup=function(peakTable,meta,sampleGroup,select) {
  whichIncl=grep(select,sampleGroup)
  peakTable=peakTable[whichIncl,]
  if(is.null(dim(meta))) meta=matrix(meta,ncol = 1)
  meta=meta[whichIncl,]
  return(list(peakTable=peakTable,meta=meta))
}
