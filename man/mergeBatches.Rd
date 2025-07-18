% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeBatches.R
\name{mergeBatches}
\alias{mergeBatches}
\title{Merge batches after drift correction}
\usage{
mergeBatches(batchList, qualRatio = 0.5)
}
\arguments{
\item{batchList}{A list of correction objects (after drift correction)}

\item{qualRatio}{Proportion of}
}
\value{
A list object containing

`peakTableOrg` A merged peak table of original data

`peakTableCorr` A merged peak table of drift-corrected data

`batch` Batch identifier (per sample)

`injection` Injection number (per sample)
}
\description{
The output of the within-batch drift correction is a correction object.
This function merges peak tables from several batches by extracting information from the correction objects.
The user must specify a minimum proportion of qualified batches per feature, i.e. such batches where the QC CV is < the specified limit.
There is thus a risk that features with poor quality (in certain batches) are present, but the features are present in high quality in sufficient proportion of batches to anyway warrant inclusion.
}
\examples{
\donttest{
data('ThreeBatchData') 
set.seed(2024)
# Get batches
batchB <- getBatch(peakTable = PTfill, meta = meta, 
                   batch = meta$batch, select = 'B')
batchF <- getBatch(peakTable = PTfill, meta = meta, 
                   batch = meta$batch, select = 'F')
# Drift correction using QCs
BCorr <- correctDrift(peakTable = batchB$peakTable, 
                      injections = batchB$meta$inj, 
                      sampleGroups = batchB$meta$grp, QCID = 'QC', 
                      G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'))
# More unbiased drift correction using QCs & external reference samples
FCorr <- correctDrift(peakTable = batchF$peakTable, 
                      injections = batchF$meta$inj,
                      sampleGroups = batchF$meta$grp, QCID = 'QC',
                      RefID='Ref', G = seq(5,35,by=3), 
                      modelNames = c('VVE', 'VEE'))
# Merge batches for batch normalization, for example
mergedData <- mergeBatches(list(BCorr, FCorr))
}

}
