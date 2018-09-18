library(batchCorr)
data('OneBatchData')
# This will load 3 objects:
# B_PT which is the one-batch dataset without any missing values (after fillPeaks/Imputation)
# B_meta which contains information on batch (B), sample group (QC or Ref) and inj (number in the injection sequence)

##########################
# In the data we have QC samples for drift correction and long-term reference samples
# that can be used for validation of drift correction
#
# We will perform two different types of drift correction
# 1. Drift correction of batch data based ONLY on QC samples
# 2. Drift correction of batch data based on QC samples and validated by reference samples

##########################
# 1. Make drift correction based on QC-samples WITHOUT validation by external reference samples
batchBQC=getGroup(peakTable=B_PT, meta=B_meta, sampleGroup=B_meta$grp, select='QC')
# Make QC, Batch and Ref objects used by batchCorr
BQC=makeQCObject(peakTable = batchBQC$peakTable, inj = batchBQC$meta$inj)
BBatch=makeBatchObject(peakTable = B_PT, inj = B_meta$inj, QCObject = BQC)
# Perform batch correction
# BCorr=driftWrap(QCObject = BQC, BatchObject = BBatch, modelNames = NULL, G = seq(5, 45, by = 10), report = TRUE) # This will take a few minutes...
BCorr=driftWrap(QCObject = BQC, BatchObject = BBatch, modelNames = 'VVE', G = seq(10, 30, by = 2), report = TRUE)

##########################
# 2. Make drift correction based on QC-samples WITH validation by external reference samples
batchBQC=getGroup(peakTable=B_PT, meta=B_meta, sampleGroup=B_meta$grp, select='QC')
batchBRef=getGroup(peakTable=B_PT, meta=B_meta, sampleGroup=B_meta$grp, select='Ref')
# Make QC, Batch and Ref objects used by batchCorr
BQC=makeQCObject(peakTable = batchBQC$peakTable, inj = batchBQC$meta$inj)
BBatch=makeBatchObject(peakTable = B_PT, inj = B_meta$inj, QCObject = BQC)
BRef=makeBatchObject(peakTable = batchBRef$peakTable, inj = batchBRef$meta$inj, QCObject = BQC)
# Perform batch correction
# BCorr=driftWrap(QCObject = BQC, BatchObject = BBatch, modelNames = NULL, G = seq(5, 45, by = 10), RefObject = BRef, report = TRUE) # This will take a few minutes...
BCorr=driftWrap(QCObject = BQC, BatchObject = BBatch, modelNames = 'VVE', G = seq(14, 30, by = 2), RefObject = BRef, report = TRUE)

##########################
# From a coding perspective, it's easy to perform correction with validation if reference samples are present
#      There is virtually no extra computational burden. The main time is consumed during clustering!
# In fact, when reference samples are available (e.g. platform long-term reference) -
#      Validation using reference samples should ALWAYS be performed!
# This reduces the likelihood of introducing bias from "false" correction when random drift is larger than systematic drift
#
# A recommendation is to perform an initial drift correction with modelNames=NULL and G=seq(1,51,by=10)
#      This will giva an indication of which modelNames are relvant to examine in greater detail
#      Then a subselection of modelNames and G can be selected (e.g. modelNames='VVE', G=20:30)

