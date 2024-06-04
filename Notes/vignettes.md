# Vignette leftovers
xcms -> the objects? Demonstrate in the Bioconductor context!
Include the package meta-information from the .docx tutorial?
The abundance plots of features in specific clusters is nice, put in vignette.
Plots in the original presentation vignette as a sort of test that functionality doesn't change. Maybe we aren't the best at writing unit tests, and if something is accidentally changed without any failure in the unit tests we can see that.
"no true samples are shipped with the sample data"?

Why dis? Reference samples are not considered sufficiently representative (limits of detection, linearity quantifiers etc.) of the sample population or "otherwise inadequate for normalization"

Where _Feature Intensity Ratio_ is the ratio of average long-term reference QC sample intensity for a feature in batches _i_ and _j_ and _Average Feature Intensity Ratio_ is the ratio of average all-sample intensity for a feature in in batches _i_ and _j_.

Ensuring orthogonal batch presence among feature alignment candidates
- two or more alignment candidates cannot be present in the same batch; candidates which are at both ends of the m/z and/or RT window relative to a windows in an another batch such that the ones in the same batch don't fit in the same m/z RT window

Decreased operator bias in clustering since it has few modifiable parameters

Don't use "reproducible features" so as not to confuse it with reproducibility at large?