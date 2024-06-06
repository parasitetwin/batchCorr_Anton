

data("OneBatchData")
data("ThreeBatchData")
# Same result with peakTable=PTnofill
batchB=getBatch(peakTable=PTfill, meta=meta, batch=meta$batch, select='B')
identical(batchB$meta$batch, B_meta$batch) # metadata are identical
identical(batchB$peakTable, B_PT) # peakTables are not identical

However, tutorial/vignette:
peakIn <- peakInfo(PT = PTnofill, sep = '@', start = 3)
alignBat <- alignBatches(peakInfo = peakIn, PeakTabNoFill = PTnofill, PeakTabFilled = PTfill, batches = meta$batch, sampleGroups = meta$grp, selectGroup = 'QC')
PT=alignBat$PTalign
batchB <- getBatch(peakTable = PT, meta = meta, batch = meta$batch, select = 'B')
# metadata and peakTable are identical!
identical(batchB$meta$batch, B_meta$batch) # metadata are identical
identical(batchB$peakTable, B_PT) # peakTables are identical

BCorr <- correctDrift(peakTable = batchB$peakTable, injections = batchB$meta$inj, sampleGroups = batchB$meta$grp, QCID = 'QC', G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'))

batchF <- getBatch(peakTable = PTfill, meta = meta, batch = meta$batch, select = 'F')
FCorr <- correctDrift(peakTable = batchF$peakTable, injections = batchF$meta$inj, sampleGroups = batchF$meta$grp, QCID = 'QC', G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'))

mergedData <- mergeBatches(list(BCorr, FCorr, HCorr))

# The example data does not include any biological samples, in which case population = "sample" does not work
normData <- normalizeBatches(peakTableCorr = mergedData$peakTableCorr, 
                             batches = meta$batch, sampleGroup = meta$grp,
                             refGroup = 'Ref', population = 'all')

normData <- normalizeBatches(peakTable = mergedData$peakTable, batches = meta$batch, sampleGroup = meta$grp, refGroup = 'Ref', population = 'sample')



peakInfo: extract m/z and rt of features

alignBatches: between-batch alignment
The PeakTabNoFill argument is used to choose features for alignment
The PeakTabFilled argument is used for batchAlign, the actual aligning of the peakTable

makeBatchObject!
- "Prepare a Batch object for drift correction"
- peakTable contains peakTables from all batches as a matrix/array
- the same for meta, which seems to be a data.frame or smth but is returned as a matrix/array
- returns list of length 2 with injection order and peakTable
- not needed if done using correctDrift()


makeQCObject!
- "Prepare a QC object for drift correction"
- gets coefficient of variance for samples
- scales the peaktable
- checks how many 
- not needed if done using correctDrift()

driftWrap!
- simply calls the low-level functions with data from makeBatchObject and makeQCObject

clust doesn't accept any arguments not present in correctDrift
driftCalc accepts additional smoothFunc and spar arguments
driftCorr accepts a driftCalc list and a table to corrected, which can also be generated from the driftCacl list. driftCorr can low-level be run without reference samples 
cleanVar: doesn't accpet any arguments not present in correctDrift

getBatch: get a part of the dataset (a specific batch)

mergeBatches: merge drift corrected batches to a single dataset
- combines features that after correction pass the criterion of QC CV <30% in at least 50% of the batches


What about between-batch alignment and normalization?
normalizeBatches
Are BatchObjects and QCObjects (not really objects, but lists) used in all functions or just drift correction? Perhaps it is confusing to call the BatchObjects and QCObjects since they aren't really objects. Rename to BatchList and QCList instead?

correct_drift: within-batch drift correction
- data preprocessing
- call to driftWrap which calls clust, driftCalc, driftCorr and cleanVar
- include these function names in the vignette
- but can the low-level functions (clust etc.) be used with some additional functionality or smth than if called by the wrapper?
- call to drift wrap with or without reference samples

normalizeBatches:
- will check if the long-term reference samples pass some criteria
"First, you’ll notice that executing the 2nd function call will generate an error, since no true samples are shipped with the example data and population = 'sample' will not be able to select any samples"
This is a problemo, must be able to execute methinks!


alignBatches.R: top-level, calls batchFlag, alignIndex, plotAlign, batchAlign (all in batchAlign5.R)
- do we want to export the low-level functions?
- batchFlag could, with modification of alignBatches, be passed NAhard, NAsoft(experimental) and quantileSoft(experimental)
- additional arguments accepted in alignIndex: flagtype (experimental)
- additional arguments accepted in plotAlign: clust, mzwidth, rtwidth
- additional arguments accepted in batchAlign: none

- in addition, the report, reporttype and reportname arguments in the functions can be used to plot only some results in a specified manner

batchAlgin5.R: peakTab, grabAlign, align, batchFlag, clusterMatrix, clustSplit, plotClust, alignIndex, plotAlign, aggregateIndex, batchAlign

batchNorm1.R: peakInfo, featComb, batchComb, refOut, refCorr

clustFuncBare3.R: cv, rmsDist, clust, driftCalc, driftCorr, cleanVar, driftWrap (wrapper)

gets.R: getBatch, getGroup

makeBatchObject.R: makeBatchObject

makeQCObject.R: makeQCObject

mergeBatches.R: mergeBatches

normalizeBatches.R: normalizeBatches, smth else?

Refactor: align, driftCalc, driftCorr, cleanVar, batchAlign, refCorr, 



