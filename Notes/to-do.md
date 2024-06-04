# To-do
**This to-do list is prioritized; prefer tasks at the top of the list. Vilhelm has an overarching sense of the order in which tasks are best performed, and updates this list accordingly**

Write unit tests

Make a filePath argument to all plotting functionality, as Bioconductor does not allow writing to working directory. It could also be useful to return a list of plots if the filePath argument is not supplied.

Make, check and/or fix examples for exported functions; at least the makeBatchObject and makeQCObject examples are not working

Implement BiocParallel parallelization for within-batch drift correction

Remove cat() and print() calls, replace with message, warning or stop where applicable

Document the data in batchCorr.R

Write more code comments for maintainability

Carl? Comment/modify the vignette


