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
  
  #Expecting that selecting multiple batches should generate error
  expect_error(getBatch(PTfill, meta, meta$batch, select=c("F", "H")))
  
  #Expecting error when trying to collect batch which doesn't exist in meta$batch
  expect_error(getBatch(PTfill, meta, meta$batch, select = "A"))
  
  #Expecting error when using a meta with less samples than there are samples in the peak table
  expect_error(getBatch(PTfill, meta[1:80], meta$batch[1:80], select="F"))
})
