# Prospective changes

All functions excluding the main three, a helper mentioned in the vignette or getBatch, getGroup, makeBatchObject, makeQCObject or mergeBatches could be non-exported functions!? What about normalization sub-functions?
- makeBatchObject and makeQCObject are only needed to manage data for drift correction, don't have to be exported methinks.
- Summa sumarum: the low-level functions are barely needed, except for driftCalc smoothFunc and spar arguments. So perhaps only export correct_drift(), maybe even implement driftWrap in there? Check the notame wrapper on this?
Notame wrappers (between-batch alignment and normalization) use the top-level alignBatches and normaliseBatches functions.
- I suggest that we decide to only do the three main functions, modifying them such that all necessary arguments can be passed to helpers, perhaps excluding "experimental" arguments? + data manipulation functions 
Or decide to include the helpers of the main functions. But that's it.

It is quite possible that the Bioconductor reviewers want to rewrite the package in terms of the data structures and other infrastructure. Do or do not prepare for this in advance?

Is there hierarchical clustering in mclust?

It could be helpful if Carl or I did a round of commenting the code.

The for-loop in normalizeBatches could be modified such that the part identifying candidates for reference sample normalization could be outside the loop

It may be a good idea to SummarizedExperimentize the package, but meibi later? I could do a round of SummarizedExperimenting at the same time as notame is translated to SummarizedExperimentized.

Reminder: Implement multiple types of reference samples in driftCorr? 

Reminder: Do we want to keep experimental arguments: NAsoft, quantileSoft, flagAll, flagSoft, flagType?