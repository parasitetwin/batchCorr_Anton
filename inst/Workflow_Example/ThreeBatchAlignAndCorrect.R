library(batchCorr)
data('ThreeBatchData')
# This will load 3 objects:
# PTfill which is the entire 3-batch dataset without any missing values (after fillPeaks/Imputation)
# PTnofill which is the 3-batch dataset with missing values (before fillPeaks/Imputation)
# meta which contains information on batch (B, F or H), sample group (QC or Ref) and inj (number in the injection sequence)


##########################
# We will perform 3 tasks:
# 1. Batch alignment to merge features that were systematically split into different features in different batches
# 2. Perform cluster-based within-batch intensity drift correction
# 3. Between-batch normalization of features (merging back features into one master peak table)




##########################
## Perform batch alignment
# Extract peakinfo (i.e. m/z and rt of features)
peakIn=peakInfo(PT = PTnofill, sep = '@', start = 3) # These column names have 2 leading characters describing LC-MS mode -> start at 3
# Perform multi-batch alignment
alignBat <- alignBatches(peakInfo = peakIn, PeakTabNoFill = PTnofill, PeakTabFilled = PTfill, batches = meta$batch, sampleGroups = meta$grp, selectGroup = 'QC')
# Extract new peak table
PT=alignBat$PTalign
dim(PT)
# [1]    90 11284


##########################
## Perform within-batch intensity drift correction
# Batch B
batchB <- getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'B')
BCorr <- correctDrift(peakTable = batchB$peakTable, injections = batchB$meta$inj, sampleGroups = batchB$meta$grp, QCID = 'QC', G = seq(5,35,by=5), modelNames = c('VVE', 'VEE') )
# Batch F
batchF <- getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'F')
FCorr <- correctDrift(peakTable = batchF$peakTable, injections = batchF$meta$inj, sampleGroups = batchF$meta$grp, QCID = 'QC', G = seq(5,35,by=5), modelNames = c('VVE', 'VEE') )
# Batch H
batchH <- getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'H')
HCorr <- correctDrift(peakTable = batchH$peakTable, injections = batchH$meta$inj, sampleGroups = batchH$meta$grp, QCID = 'QC', G = seq(5,35,by=5),modelNames = c('VVE', 'VEE') )

##########################
## Perform between-batch normalization
mergedData <- mergeBatches(list(BCorr,FCorr,HCorr))
normData <- normalizeBatches(peakTableCorr = mergedData$peakTableCorr, peakTableOrg = mergedData$peakTableOrg, batches = meta$batch, sampleGroup = meta$grp, refGroup = 'Ref', population = 'all')
# normData <- normalizeBatches(peakTable = mergedData$peakTable, batches = meta$batch, sampleGroup = meta$grp, refGroup = 'Ref', population = 'all')
PTnorm <- normData$peakTable

