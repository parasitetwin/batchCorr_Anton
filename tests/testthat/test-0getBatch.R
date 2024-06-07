context("getBatch")
test_that("getBatch",{
  
  data("ThreeBatchData", package="batchCorr")
  
  #Check that batch collected using getBatch is correct
  batchF <- getBatch(PTfill,meta,meta$batch,"F")
  
  #Checking that no other batches were extracted due to mismatch in meta and PT + that number of samples is correct
  batchVector <- c()
  for(i in 1:nrow(batchF$peakTable)){
    batchVector <- c(batchVector,
                     strsplit(rownames(batchF$peakTable)[i], "_")[[1]][3])
  }
  expect_true(length(batchVector) == nrow(batchF$meta) && nrow(batchF$meta) == nrow(PTfill[grepl("BatchF", rownames(PTfill)), ]))
  
  #Expecting false return values when input is a vector - Function causes bugs when this is done
  dimFandH <- dim(getBatch(PTfill, meta, meta$batch, c("F", "H"))$peakTable)[1]
  dimF <- dim(getBatch(PTfill, meta, meta$batch, c("F"))$peakTable)[1]
  dimH <- dim(getBatch(PTfill, meta, meta$batch, c("H"))$peakTable)[1]
  expect_false(dimFandH == dimF && dimFandH == dimH)
  
})
