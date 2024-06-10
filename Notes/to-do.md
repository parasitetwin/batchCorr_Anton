# To-do
**This to-do list is prioritized; prefer tasks at the top of the list. Vilhelm has an overarching sense of the order in which tasks are best performed, and updates this list accordingly. You could copy the task name as a commit title.**

Read through the Bioconductor submission guidelines to see if there's something important to add!

R CMD check error: no test files
- Write unit tests, primarily for exported functions

R CMD check warning:
Undocumented code objects:
  ‘B_PT’ ‘B_meta’ ‘PTfill’ ‘PTnofill’ ‘meta’
Undocumented data sets:
  ‘B_PT’ ‘B_meta’ ‘PTfill’ ‘PTnofill’ ‘meta’
R CMD check note: 
no visible global function definition for many code objects
- run R CMD check on the tar.gz, Vilhelm?

BiocCheck::BiocCheck warning: No Bioconductor dependencies detected. Note that some infrastructure packages may not have Bioconductor dependencies. For more information, reach out to the Bioconductor community and/or consider a CRAN submission.

R CMD check note: installed size is 15.9 mb. Subdirectories of 1mb or more: data 15.0 mb
BiocCheck::BiocCheck warning: Data files exceed the 5MB size limit
- One may possibly remove onebatch data, but it's still like 12.7 mb, any suggestions Carl?

BiocCheck::BiocCheck error: Unable to find your email in the Support Site (Anton)

BiocCheck::BiocCheck note: Cannot determine whether maintainer is subscribed to the Bioc-Devel mailing list (requires admin credentials). Subscribe here: https://stat.ethz.ch/mailman/listinfo/bioc-devel
- Anton? 

BiocCheck::BiocCheck note: The recommended function length is 50 lines or less. There are 8 functions greater than 50 lines.
The longest 5 functions are:
driftCorr() (R/clustFuncBare3.r): 105 lines
normalizeBatches() (R/normalizeBatches.R):  85 lines
driftCalc() (R/clustFuncBare3.r):  80 lines
align() (R/batchAlign5.R):  79 lines
refCorr() (R/batchNorm1.r):  67 lines

BiocCheck::BiocCheck note: Avoid 1:...; use seq_len() or seq_along()

BiocCheck::BiocCheck note: Avoid 'cat' and 'print' outside of 'show' methods

BiocCheck::BiocCheck note: Avoid using '=' for assignment and use '<-' instead

Make a filePath argument to all plotting functionality, as Bioconductor does not allow writing to working directory. It could also be useful to return a list of plots if the filePath argument is not supplied.

Implement BiocParallel parallelization for within-batch drift correction

Remove cat() and print() calls, replace with message, warning or stop where applicable

Document the data in batchCorr.R

Write more code comments for maintainability

Comment/modify the vignette, Carl?

Formatting (at end):
R CMD check note:
examples lines wider than 100 characters (80 for Bioconductor!)
NOTE: Consider shorter lines; 251 lines (15%) are > 80 characters long.
BiocCheck::BiocCheck note: 
Consider shorter lines; 250 lines (15%) are > 80 characters long.
BiocCheck::BiocCheck note:
Consider 4 spaces instead of tabs; 276 lines (17%) contain tabs.

After Vilhelm, in the vignette, removed the mentions of the low-level functions as they are no longer exported, it might be good if Carl read through it and made sure it's good!




