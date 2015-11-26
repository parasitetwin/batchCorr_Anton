# batchCorr
Within and between batch correction of LC-MS metabolomics data

This is a repository containing functions within three areas of batch correction. These algorithms were originally developed 
to increase quality and information content in data from LC-MS metabolomics. However, the algorithms should be applicable to 
other data structures/origins, where within and between batch irregularities occur.

The three areas indicated are:
area | abbreviation | description
:--- | :----------- | :----------
Batch alignment | BA | Functions to align features that are originally systematically misaligned across batches
Drift correction | DC | Functions to perform within batch intensity drift correction
Batch normalisation | BN | Funtions to perform between batch normalisation

Batch alignment is achieved based on three concepts:
-Aggregation of feature presence/missingness on batch level.
-Identifying features with missingness within "the box", i.e. sufficiently similar in retention time and m/z.
-Ensuring orthogonal batch presence among feature alignment candidates.

Drift corrections is achieved based on:

Batch normalisation is achieved based on:

These development and inner workings of these algorithms are reported in:

Brunius C, Shi L and Landberg R. Within and between batch correction of LC-MS metabolomics data. Submitted manuscript.

