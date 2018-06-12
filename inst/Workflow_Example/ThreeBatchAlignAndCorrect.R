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
peakIn=peakInfo(PT = PTnofill, sep = '@', start = 3) # column names have 2 leading characters describing LC-MS mode -> start at 3
# Flag presence/missingness on batch level
bF=batchFlag(PTnofill = PTnofill, batch = meta$batch, sampleGroup = meta$grp, peakInfo = peakIn)
# Find possible alignment candidates per sample type
aIQ=alignIndex(batchflag = bF, grpType='QC', mzdiff=0.002, rtdiff=15, report=T, reportName='splits_aIQ')
# Plot achieved alignments (not necessary)
plotAlign(batchflag = bF, alignindex = aIQ, plotType='pdf', reportName='clustPlots_aIQ')
# Perform alignment -> Peaktable
bA=batchAlign(batchflag = bF, alignindex = aIQ, peaktable_filled = PTfill, batch = meta$batch)
# Extract new peak table
PT=bA$PTalign
dim(PT)
# [1]    90 11284


##########################
## Perform within-batch intensity drift correction

# Batch B
# Extract sub-parts of the peak tables and metadata
batchB=getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'B')
batchBQC=getGroup(peakTable=batchB$peakTable, meta=batchB$meta, sampleGroup=batchB$meta$grp, select='QC')
batchBRef=getGroup(peakTable=batchB$peakTable, meta=batchB$meta, sampleGroup=batchB$meta$grp, select='Ref')
# Make QC, Batch and Ref objects used by batchCorr
BQC=makeQCObject(peakTable = batchBQC$peakTable, inj = batchBQC$meta$inj)
BBatch=makeBatchObject(peakTable = batchB$peakTable, inj = batchB$meta$inj, QCObject = BQC)
BRef=makeBatchObject(peakTable = batchBRef$peakTable, inj = batchBRef$meta$inj, QCObject = BQC)
# Perform batch correction
BCorr=driftWrap(QCObject = BQC, BatchObject = BBatch, RefObject = BRef)

# Batch F
# Extract sub-parts of the peak tables and metadata
batchF=getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'F')
batchFQC=getGroup(peakTable=batchF$peakTable, meta=batchF$meta, sampleGroup=batchF$meta$grp, select='QC')
batchFRef=getGroup(peakTable=batchF$peakTable, meta=batchF$meta, sampleGroup=batchF$meta$grp, select='Ref')
# Make QC, Batch and Ref objects used by batchCorr
FQC=makeQCObject(peakTable = batchFQC$peakTable, inj = batchFQC$meta$inj)
FBatch=makeBatchObject(peakTable = batchF$peakTable, inj = batchF$meta$inj, QCObject = FQC)
FRef=makeBatchObject(peakTable = batchFRef$peakTable, inj = batchFRef$meta$inj, QCObject = FQC)
# Perform batch correction
FCorr=driftWrap(QCObject = FQC, BatchObject = FBatch, RefObject = FRef)

# Batch H
# Extract sub-parts of the peak tables and metadata
batchH=getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'H')
batchHQC=getGroup(peakTable=batchH$peakTable, meta=batchH$meta, sampleGroup=batchH$meta$grp, select='QC')
batchHRef=getGroup(peakTable=batchH$peakTable, meta=batchH$meta, sampleGroup=batchH$meta$grp, select='Ref')
# Make QC, Batch and Ref objects used by batchCorr
HQC=makeQCObject(peakTable = batchHQC$peakTable, inj = batchHQC$meta$inj)
HBatch=makeBatchObject(peakTable = batchH$peakTable, inj = batchH$meta$inj, QCObject = HQC)
HRef=makeBatchObject(peakTable = batchHRef$peakTable, inj = batchHRef$meta$inj, QCObject = HQC)
# Perform batch correction
HCorr=driftWrap(QCObject = HQC, BatchObject = HBatch, RefObject = HRef)

