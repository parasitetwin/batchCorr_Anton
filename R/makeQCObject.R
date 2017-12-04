#' Prepare a QC object for drift correction
#'
#' @param QCPeakTab Peak table (without missing values). Samples in rows; Features in columns; colnames as mz@rt
#' @param inj Injection sequence number for QC samples
#'
#' @return QC object
#' @export
#'
#' @examples QCObject=makeQCObject(QCPeakTable, QCInjections)
makeQCObject=function (QCPeakTab, inj) 
{
  QCCV = batchCorr::cv(QCPeakTab)
  QCscale = scale(QCPeakTab, center = FALSE)
  NAs = colSums(is.na(QCscale)) > 0
  QCRawNaRm = QCPeakTab[, !NAs]
  QCFeats = QCscale[, !NAs]
  return(list(inj = inj, Feats = QCFeats, RawFeats = QCPeakTab, 
              RawFeatsNaRm = QCRawNaRm, NAs = NAs))
}
