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
B_CorrWithoutRef <- correctDrift(peakTable = B_PT, injections = B_meta$inj, sampleGroups = B_meta$grp, QCID = 'QC', modelNames = c('VVE','VVV'), G = seq(5,25,by=5))

##########################
# 2. Make drift correction based on QC-samples WITH validation by external reference samples
B_CorrWithRef <- correctDrift(peakTable = B_PT, injections = B_meta$inj, sampleGroups = B_meta$grp, QCID = 'QC', RefID = 'Ref', modelNames = c('VVE','VVV'), G = seq(5,25,by=5))

##########################
# From a coding perspective, it's easy to perform correction with validation if reference samples are present
#      There is virtually no extra computational burden. The main time is consumed during clustering!
# In fact, when reference samples are available (e.g. platform long-term reference) -
#      Validation using reference samples should ALWAYS be performed!
# This reduces the likelihood of introducing bias from "false" correction when random drift is larger than systematic drift
