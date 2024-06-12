library(stringr)
library(testthat)

test_that("peakInfo", {
  
  data("ThreeBatchData", package="batchCorr")
  
  #Checking that number of columns is the same as rows in the returned DF
  mz_rt <- peakInfo(PT=PTfill, sep="@", timepos=2, start=3)
  
  expect_true(nrow(mz_rt) == ncol(PTfill))
  
  #Checking that columns correspond to rows correctly
  rtSum <- 0
  mzSum <- 0
  for(i in 1:nrow(mz_rt)){
    if(grepl(mz_rt[i,1],
                    colnames(PTfill)[i])){
      mzSum <- mzSum+1
    }
    if(grepl(mz_rt[i,2],
             colnames(PTfill)[i])){
      rtSum <- rtSum + 1
    }
  }
  
  expect_true(mzSum == ncol(PTfill) && rtSum == ncol(PTfill))
  
  #Checking that wrong separator generates an error
  expect_error(mz_rt <- peakInfo(PT=PTfill, sep="_", timepos=2, start=3))
  
  #Checking that character "start" generates an error
  expect_error(mz_rt <- peakInfo(PT=PTfill, sep="@", timepos=1, start="A"))
  
  #Checking that character "timepos" generates an error
  expect_error(mz_rt <- peakInfo(PT=PTfill, sep="@", timepos="1", start=3))
  
  #Checking that character "timepos" and "start" generates an error
  expect_error(mz_rt <- peakInfo(PT=PTfill, sep="@", timepos="1", start="3"))
  
  
})
