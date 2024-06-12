test_that("normalizeBatches", {
  
  data("ThreeBatchData", package="batchCorr")
  
  #Normalizing batches previously corrected with corrDrift and merged with mergeBatches
  meta <- meta[-which(meta$batch=="H"),]
  normData <- normalizeBatches(peakTableCorr = mergedData$peakTableCorr,
                               batches = meta$batch,
                               sampleGroup = meta$grp,
                               refGroup = 'Ref',
                               population = 'QC')
  
  #Expecting error when using a population name which doesn't exist
  expect_error(normalizeBatches(peakTableCorr = mergedData$peakTableCorr,
                                batches = meta$batch,
                                sampleGroup = meta$grp,
                                refGroup = 'Ref',
                                population = 'error'))
  
  #Expecting error when using a refGroup which doesn't exist
  expect_error(normalizeBatches(peakTableCorr = mergedData$peakTableCorr,
                                batches = meta$batch,
                                sampleGroup = meta$grp,
                                refGroup = 'error',
                                population = 'QC'))
  
  #Expecting error when using a batch-vector of different length
  expect_error(normalizeBatches(peakTableCorr = mergedData$peakTableCorr,
                                batches = meta$batch[1:5],
                                sampleGroup = meta$grp,
                                refGroup = 'Ref',
                                population = 'QC'))
  
  #Expecting error when using a sampleGroup-vector of different length
  expect_error(normalizeBatches(peakTableCorr = mergedData$peakTableCorr,
                                batches = meta$batch,
                                sampleGroup = meta$grp[1:5],
                                refGroup = 'Ref',
                                population = 'QC'))
  
  
})