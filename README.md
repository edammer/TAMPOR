# Tunable Approach for Median Polish of Ratio -- Biological Abundance Batch Effect Removal
## (TAMPOR), by Eric Dammer
For public use with citation of origin, and original author
####
## Purpose
####
#### - Robustly removes batch artifacts and batchwise variance from abundance data.
#### - Works well even if batch replicates for normalization have unusual variance compared to biological samples not used for normalization (e.g., due to differential peptide digestion, different tissue region, or different genetic background of control samples.)
#### - See <A href="https://github.com/edammer/TAMPOR/blob/master/walkthrough.md">walkthrough using sample 3 batches x 10 protein data</a>, an example of TAMPOR use case 3.
####
## General Use Cases
####
#### 1. Batch replicates (e.g., global internal standard, GIS, in TMT proteomics) if present are baseline samples considered as ratio denominator for removal of batch artifacts.
#### 2. Consider biologically similar samples as baseline for the denominator.    [set a GIS column in traits input]
#### 3. Third, baseline samples present can be used, and all other samples can also be used. This assumes batches are generally trait-balanced and randomized, but GIS and non-GIS samples can be grossly different.    [useAllNonGIS=TRUE]
#### 4. Last, if batches consist of trait-randomized and generally trait-balanced biological samples, all samples can be considered for baseline/denominator.    [noGIS=TRUE]
####
####
####  **_inputs_**: 
| data frame or matrix              | Description                                       |
|:----------------------------------|:--------------------------------------------------|
| 1. (normalized or RAW) abundance  | columns named by samples, rows by measured species|
| 2. traits / phenotypes            | columns named by traits, rows named by sample     |

 *Sample names as column and row names in abundance and traits, respectively, must match exactly.*
####
####  **_outputs_**:  (as list, with the following elements)
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
|     cleanDat.oneIter     | log2 of naive ratio, applying algorithm without iteration      |
|     cleanRelAbun.oneIter | naive ratio, converted back to abundance                       |
|     converged            | TRUE/FALSE, was convergence reached in the median polishing?   |
|     iterations           | Number of median polish iterations performed                   |
|     **convergencePlots** | **Convergence plots of Frobenius norm 2-iteration differences**|
|     **meanSDplots**      | **vsn package mean-SD plots of variance for:**                 |
|                          |   a. input abundance of same dimensions as output              |
|                          |   b. cleanRelAbun.oneIter                                      |
|                          |   c. corrected abundance (cleanRelAbun)                        |
|     **MDSplots**         | **limma multi-dimensional scaling plots for:**                 |
|                          |   a. input abundance with same dimensions as output            |
|                          |   b. cleanRelAbun.oneIter                                      |
|                          |   c. corrected abundance (cleanRelAbun)                        |

**The 3 last list elements above render the graphics of the PDF output pages.**

     If removeGISafter option not enabled, variance and MDS plots each display with and without GIS.

##### Writes PDF visualizing variance and MDS for (a) original, (b) first iteration, and (c) final output abundances to file:
     /path/TAMPOR-Improvement.Vis(#iterations)-#ROWSx#COLUMNS_outputSuffix.pdf
     
####
#### Depends on these 5 packages and their dependencies:
 limma, vsn, doParallel, ggplot2, ggpubr
####
```R
options(stringsAsFactors = FALSE)

# TAMPOR parameters (full list)
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

path=getwd()       # Defaults to working directory path. This folder will contain output PDF of visualizations.
#####################################################################################

source("TAMPOR.R")

# Use cases 1 & 2
TAMPORlist.GIS <- TAMPOR(dat, traits, noGIS=FALSE, useAllNonGIS=FALSE, GISchannels=c("126","131"),
                  parallelThreads=parallelThreads, outputSuffix="GISonly")
# Use case 3
TAMPORlist.GIShybrid <- TAMPOR(dat, traits, noGIS=FALSE, useAllNonGIS=TRUE,GISchannels=c("126","131"),
                    batchPrefixInSampleNames=TRUE, parallelThreads=parallelThreads, outputSuffix="GIShybrid")
# Use case 4
TAMPORlist.noGIS <- TAMPOR(dat, traits, noGIS=TRUE, batchPrefixInSampleNames=TRUE, parallelThreads=parallelThreads,
                    outputSuffix="noGIS")
```
####
####  **_Sample Output:_**
### Multidimensional Scaling
![alt text](https://github.com/edammer/TAMPOR/blob/master/MDS50batchTMT.jpg "MDS improvement of 50 batch TMT proteomics data")
####
### Mean-SD Variance Plots
![alt text](https://github.com/edammer/TAMPOR/blob/master/meanSD50batchTMT.jpg "Variance removal visualized for same 50 batch TMT proteomics data")
####
### Convergence Tracking
![alt text](https://github.com/edammer/TAMPOR/blob/master/convergence50batchTMT.jpg "Convergence tracked for same 50 batch TMT brain proteome.")
####
