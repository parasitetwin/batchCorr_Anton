test_that("mergeBatches",{
  
  data("ThreeBatchData", package="batchCorr")
  
  #Merging datasets previously corrected in corrDrift unit test
  mergedData <- mergeBatches(list(BCorr, FCorr))
  
  #Checking that all samples were included in the merging
  expect_true((nrow(FCorr$TestFeatsFinal) + nrow(BCorr$TestFeatsFinal)) ==
                nrow(mergedData$peakTableCorr))
  
  #Checking error when inputing character in qualRatio
  expect_error(mergeBatches(list(BCorr, FCorr), qualRatio="A"))
  
  #Checking error when only supplying one batch
  expect_error(mergeBatches(FCorr))
})