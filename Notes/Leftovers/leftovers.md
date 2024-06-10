# Leftovers
Aha, the modelnames is an identifier for the parametrization of the cluster shapes consisting of three letters. The first letter is for volume, second for shape and third for orientation. Each position is coded as "E" for equal or "V" for variable, for example VVE. G is simply which numbers of clusters to try.

coefficient of variation (CV) = relative standard deviation (RSD) except for that relative standard deviation is always absolute

oldworkflow.R doesn't work anymore, bcause of data(PTfill), for example. Meaning that the low-level normalizeBatches functions which are not called from normalizeBatches used there are also old, should we remove them and just rely on normalizeBatches?

it doesn't really make sense to pass reportname as the reports from the different subfunctions would have the same name :))
