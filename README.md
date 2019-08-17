# Tunable Approach for Median Polish of Ratio -- Biological Abundance Batch Effect Removal
## (TAMPOR), by Eric Dammer
For public use with citation of origin, and original author
####
## Purpose
####
#### - Robustly removes batch artifacts and batchwise variance from abundance data.
#### - Works well even if batch replicates for normalization have unusual variance compared to biological samples not used for normalization (e.g., due to differential peptide digestion, different tissue region, or different genetic background of control samples.)
####
## General Use Cases
####
#### 1. Batch replicates (e.g., global internal standard, GIS, in TMT proteomics) if present are baseline samples considered as ratio denominator for removal of batch artifacts.
#### 2. Consider biologically similar samples as baseline for the denominator.    [set a GIS column in traits input]
#### 3. Third, baseline samples present can be used, and all other samples can also be used. This assumes batches are generally trait-balanced and randomized, but GIS and non-GIS samples can be grossly different.    [useAllNonGIS=TRUE]
#### 4. Last, if batches consist of trait-randomized and generally trait-balanced biological samples, all samples can be considered for baseline/denominator.    [noGIS=TRUE]
####
####
####  *inputs*: 
1. (normalized or RAW) abundance
2. traits, each with sample names
####
####  *outputs*:  (as list, with the following elements)
####
### FOR DOWNSTREAM ANALYSIS
| list element | Description|
| ------------:|:-----------|
| cleanRelAbun | abundance within-batch across batch batch-corrected|
|     cleanDat | log2(cleanRelAbun/row's central tendency of abundance) |
|              | (rows with 50%+ missing batches/values removed from the above.)|
|       traits | case-sample traits, of used samples, sorted and annotated |
####
### FOR COMPARISON/DIAGNOSTICS, VISUALIZATION
|        list element      | Description|
|     --------------------:|:-----------|
|     cleanDat.oneIter     | log2 of naive ratio, applying algorithm without iteration|
|     cleanRelAbun.oneIter | naive ratio, converted back to abundance|
|     converged            | TRUE/FALSE, was convergence reached in the median polishing?|
|     iterations           | Number of median polish iterations performed|
|     convergencePlot      | Convergence plot (Frobenius norm 2-iteration difference)|
|     logConvergencePlot   | log10(Frobenius norm difference) convergence plot|
|     varPlot.input        | mean-SD plot of input abundance of same dimensions as output|
|     varplot.oneIter      | mean-SD plot of cleanRelAbun.oneIter|
|     varPlot.cleanRelAbun | mean-SD plot of corrected abundance|
|     MDSplot.input        | limma MDS plot of input abundance with same dimensions as output|
|     MDSplot.oneIter      | MDS plot of cleanRelAbun.oneIter|
|     MDSplot.cleanRelAbun | MDS plot of corrected abundance|
|                          |                                |
     *.noGIS versions of variance and MDS plots are also possible list elements.
####
#### Depends on these 5 packages and their dependencies:
 limma, vsn, doParallel, ggplot2, ggpubr
####
```R
`options(stringsAsFactors = FALSE)

# TAMPOR parameters (full list of 11 possible)
#####################################################################################
outputSuffix = "TAMPOR"      #Filename Suffix for some outputs, and name of automatically created subfolder for output
noGIS=FALSE        #If TRUE, the next setting and GISchannels settings are ignored, and all samples in each batch are taken as GIS
                   #(assumes batch randomization with trait balance)
useAllNonGIS=FALSE #if TRUE, all randomized samples in the batch (except any that are GIS) are used for step 1b
                   #(TAMPOR equation1, second term) rowwise multiplier calculation
                   #if FALSE, only GIS channels are used
batchPrefixInSampleNames=FALSE #If TRUE, sample names should have format "batch.channel.(...)"
                               #If FALSE, it is recommended you have a GIS column in traits provided,
                               #with GIS samples specified "GIS" or TRUE and other samples NA
GISchannels=c("126","131")     #TMT batch channel(s) that are GIS controls in every batch where they appear, or full unique
                               #sample name(s) from traits rownames and abundance column names. At least one should be present
                               #in every batch; additive with specification of GIS samples in provided traits GIS column.
iterations=250     #How many times (max) will abundance matrix (dat) be subjected to 2-way table median polish?
                   #If convergence is reached sooner, not all iterations will be run.
samplesToIgnore=c("NONE")      #These samples have all values initially set to NA, & are removed from matrices used for visualization.
                               #Useful for comparison of output without known outliers to prior output. If none, specify any string
                               #not equal to a column name of the inputAbundanceCSV; E.g. "NONE"
meanOrMedian="median"  #must be a valid R function, e.g. 'mean' or 'median'; median is recommended unless robustness to outliers
                       #is not required.
removeGISafter=FALSE   #Should samples designated as GIS be removed before visualization of variance and MDS?
minimumBatchSize=5     #batches with fewer samples will not be used or kept in data for batch variance removal [default=5]
parallelThreads=8      #doParallel library will run this many threads locally to split steps 1a and 1b batchwise calculations,
                       #speeding processing time [default=2]
#####################################################################################

# Use cases 1 & 2
TAMPORlist.GIS <- TAMPOR(dat, traits, noGIS=FALSE, useAllNonGIS=FALSE, GISchannels=c("126","131"), parallelThreads=parallelThreads,
                  outputSuffix="GISonly")
# Use case 3
TAMPORlist.GIShybrid <- TAMPOR(dat.real, traits.real, noGIS=FALSE, useAllNonGIS=TRUE,GISchannels=c("126","131"),
                    batchPrefixInSampleNames=TRUE, parallelThreads=parallelThreads, outputSuffix="GIShybrid",path=outputfigs)
# Use case 4
TAMPORlist.noGIS <- TAMPOR(dat, traits, noGIS=TRUE, batchPrefixInSampleNames=TRUE, parallelThreads=parallelThreads,
                    outputSuffix="noGIS")
`
