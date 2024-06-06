#' Perform multi-batch alignment to merge features artificially split between batches
#'
#' @param peakInfo matrix with mz and rt in columns 1:2 (see e.g. ?peakInfo)
#' @param PeakTabNoFill Multi-batch peak table including missing values
#' @param PeakTabFilled Multi-batch peak table without missing values
#' @param batches Vector (length=nrow(PeakTab)) of batch identifiers
#' @param sampleGroups Vector (length=nrow(PeakTab)) of sample type (e.g. "sample", "QC", "Ref)
#' @param selectGroup Which sample type to base alignment on (e.g. "sample", "QC" or "Ref")
#' @param mzdiff Tolerance for difference in m/z
#' @param rtdiff Tolerance for difference in retention time
#' @param report Whether to export diagnostic plots into your work directory (defaults to TRUE)
#'
#' @return batchAlign Object
#' @return Aligned peak table available under $PTalign
#' @return Returns NULL if no alignment candidates can be found
#' @export
#'
#' @examples
#' data('ThreeBatchData')
#' # Extract peakinfo (i.e. m/z and rt of features). These column names have 2
#' # leading characters describing LC-MS mode -> start at 3
#' peakIn <- peakInfo(PT = PTnofill, sep = '@', start = 3)  
#' # Perform multi-batch alignment
#' alignBat <- alignBatches(peakInfo = peakIn, PeakTabNoFill = PTnofill,
#'                          PeakTabFilled = PTfill, batches = meta$batch,
#'                          sampleGroups = meta$grp, selectGroup = 'QC')
#' # Extract new peak table
#' PT <- alignBat$PTalign 
alignBatches=function(peakInfo, PeakTabNoFill, PeakTabFilled, batches, sampleGroups, selectGroup, mzdiff=0.002, rtdiff=15, report=TRUE) {
  bFlag <- batchFlag(PTnofill = PeakTabNoFill, batch = batches, sampleGroup = sampleGroups, peakInfo = peakInfo)
  aIQ <- alignIndex(batchflag = bFlag, grpType=selectGroup, mzdiff=mzdiff, rtdiff=rtdiff, report=report)
  if(is.null(aIQ)) return(NULL)
  if (report) plotAlign(batchflag = bFlag, alignindex = aIQ, plotType='pdf')
  bAlign <- batchAlign(batchflag = bFlag, alignindex = aIQ, peaktable_filled = PeakTabFilled, batch = batches)
}
