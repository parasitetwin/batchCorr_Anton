#' Prepare a QC object for drift correction
#'
#' @param peakTable Peak table (without missing values). Samples in rows; Features in columns; colnames as mz@rt
#' @param inj Injection sequence number for QC samples
#'
#' @return QC object
#'
#' @examples QCObject=makeQCObject(QCPeakTable, QCInjections)
#'
#' @noRd
makeQCObject=function (peakTable, inj)
{
  if(length(inj)!=nrow(peakTable)) stop ('mismatch number of samples in peak table and injection sequence')
  QCCV = cv(peakTable)
  QCscale = scale(peakTable, center = FALSE)
  NAs = colSums(is.na(QCscale)) > 0
  QCRawNaRm = peakTable[, !NAs]
  QCFeats = QCscale[, !NAs]
  return(list(inj = inj, Feats = QCFeats, RawFeats = peakTable,
              RawFeatsNaRm = QCRawNaRm, NAs = NAs))
}
