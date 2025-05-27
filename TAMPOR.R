TAMPOR <- function (dat, traits, noGIS = FALSE, useAllNonGIS = FALSE, batchPrefixInSampleNames = FALSE, GISchannels = "GIS", iterations = 250, skipMDS = FALSE, sampleMedianRows = "ALL", fractionNAmax = 0.500,
    samplesToIgnore = FALSE, meanOrMedian = "median", removeGISafter = FALSE, minimumBatchSize = 5, parallelThreads=2, outputSuffix = "TAMPOR", path=getwd()) {

cleanDat <- dat
outputfigs <- outputtabs <- path

require(doParallel, quietly=TRUE)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)

require(ggplot2, quietly=TRUE)
require(ggpubr, quietly=TRUE)
require(vsn, quietly=TRUE)
require(limma, quietly=TRUE)


if(Sys.getenv("RSTUDIO")=="1") cat("\nDetected TAMPOR run within RStudio.\nBe sure that the pane with your plots tab covers at least 1/2 the width of your screen and 3/4 of its height\nto avoid plot capture failure and failure to retain TAMPOR output.\n")

# Standardize "Batch" column in traits
if ("batch" %in% colnames(traits) | "Batch" %in% colnames(traits) | "BATCH" %in% colnames(traits)) {
  if ("batch" %in% colnames(traits)) colnames(traits)[which(colnames(traits)=="batch")] <- "Batch"
  if ("BATCH" %in% colnames(traits)) colnames(traits)[which(colnames(traits)=="BATCH")] <- "Batch"
} else {
  cat("WARNING: no 'Batch' column found in traits file. Failing over to batchPrefixInSampleNames=TRUE.\n")
  batchPrefixInSampleNames=TRUE
}

# Check that all sample names for abundance are found in traits
if (!identical(colnames(cleanDat),na.omit(rownames(traits)[na.omit(match(colnames(cleanDat),rownames(traits)))]))) {
  missingSampleTraits<-sort(setdiff(colnames(cleanDat), na.omit(rownames(traits)[na.omit(match(colnames(cleanDat),rownames(traits)))])))
  cat(paste0("WARNING: Missing traits for the following sample(s) in provided abundance data:\n",paste(missingSampleTraits,collapse=", "),"\n"))
  cat("         These samples have been removed.\n")
  cleanDat<-cleanDat[,-match(missingSampleTraits,colnames(cleanDat))]
}

# Check for and remove extra traits not in abundance data
if (length(intersect(rownames(traits),colnames(cleanDat))) < length(rownames(traits))) {
  extraSampleTraits<-sort(setdiff(rownames(traits), na.omit(colnames(cleanDat)[na.omit(match(rownames(traits),colnames(cleanDat)))])))
  cat(paste0("WARNING: Sample(s) in traits not found in provided abundance data:\n",paste(extraSampleTraits,collapse=", "),"\n"))
  cat("         These samples have been removed from traits.\n")
  traits<-traits[-match(extraSampleTraits,rownames(traits)),]
}


# Annotate which samples go with a batch, and identifiers for each sample, which can be checked against GISchannels list
if (batchPrefixInSampleNames) {
  if(!length(which(grepl("\\.",rownames(traits))))==nrow(traits)) {
    stop("ERROR: sample names in both input files expected to have at least one '.' separating batch.channel! (Not found in traits sample/row names.)\n")
  } else {
    sampleIndex <- as.data.frame(do.call(rbind,strsplit(rownames(traits),"[.]")))[,1:2]
    colnames(sampleIndex) <- c("batch", "channel")
    traits$Batch <- sampleIndex$batch
  }
}
  


# Check for possible GIS samples in other situations
if(!length(which(grepl("\\.",rownames(traits))))==nrow(traits)) {
  channelCandidates <- rownames(traits)
} else {
  channelCandidates <- as.data.frame(do.call(rbind,strsplit(rownames(traits),"[.]")))[,2]
}

if ( (!"GIS" %in% colnames(traits) & length(na.omit(match(GISchannels,channelCandidates)))==0 & length(na.omit(match(GISchannels,rownames(traits))))==0) 
   | ("GIS" %in% colnames(traits) & (length(which(names(table(traits$GIS))=="GIS"))+length(which(as.logical(names(table(traits$GIS))))))==0) ) {
  if (!noGIS) cat("WARNING: no GIS samples specified or not found if specified. Failing over to noGIS=TRUE!\n")
  noGIS=TRUE
  traits$GIS=rep("GIS",nrow(traits))
  GISchannels="GIS"
}
if (length(which(as.logical(names(table(traits$GIS)))))>=1) traits$GIS[which(as.logical(traits$GIS))]<-"GIS"
if (!exists("sampleIndex")) {
  if (!"GIS" %in% colnames(traits)) {
    sampleIndex <- data.frame(batch=traits$Batch,channel=rownames(traits))
    cat("NOTE: Sample names from traits will be checked for exact matching to GISchannels specified.\n")
  } else {
    sampleIndex <- data.frame(batch=traits$Batch,channel=traits$GIS)
  }
}


# Match abundance (cleanDat) columns to batch-ordered traits
sampleIndex<-sampleIndex[order(traits$Batch),]
rownames(sampleIndex)<-NULL   #renumber after ordering.
traits<-traits[order(traits$Batch),]
cleanDat<-cleanDat[,na.omit(match(rownames(traits),colnames(cleanDat)))]
cleanDat.original <- cleanDat 

batchIndex <- as.character(unique(sampleIndex$batch))
if (!"Batch" %in% traits) traits$Batch<-sampleIndex$batch


# Ignore samples in samplesToIgnore vector (set all their values to NA)
if (length(na.omit(match(samplesToIgnore, colnames(cleanDat)))) == length(samplesToIgnore)) {
  for (i in 1:length(samplesToIgnore)) {
    cleanDat[, samplesToIgnore[i]] <- as.vector(rep(NA, nrow(cleanDat)))
  }
} else {
  cat("NOTE: one or more samplesToIgnore do not match sample names (colnames) in input abundance data.\n      Not ignoring any samples.\n")
}


# Set parameters for runs with no GIS
if (noGIS) {
  if (useAllNonGIS) cat("WARNING: noGIS and useAllNonGIS options were both enabled. Keeping with noGIS enabled, useAllNonGIS disabled.\n")
  useAllNonGIS=FALSE
  GISchannels=rownames(traits) #all channels are specified in all batches for all denominators in equation 1.
  sampleIndex$channel<-rep("GIS",nrow(sampleIndex))

  traits$GIS<-sampleIndex$channel
}


# Check sample names, channels in all batches to confirm at least one of GISchannels exists in each batch
minimumBatchSize<-as.integer(minimumBatchSize)
batchwiseSampleCounts=sapply(as.character(batchIndex),function(batch) length(which(sampleIndex$batch==batch)))
batchwiseGIScandidateCounts=sapply(as.character(batchIndex),function(batch) length(which(any(c(GISchannels, "GIS") %in% sampleIndex$channel[sampleIndex$batch==batch]))) + length(which(any(c(GISchannels,"GIS") %in% rownames(traits)[sampleIndex$batch==batch]))) + length(which(traits$GIS[sampleIndex$batch==batch]=="GIS")))
names(batchwiseGIScandidateCounts) <- names(batchwiseSampleCounts) <- as.character(batchIndex)
if (!any( batchwiseGIScandidateCounts == 0) & !any( batchwiseSampleCounts < minimumBatchSize)) {
  if (!noGIS) cat("NOTE: Successfully checked that each batch has at least one sample designated as GIS control replicate.\n")
} else {
  badBatches<-as.character(batchIndex)[which(batchwiseGIScandidateCounts==0)]
  badBatches<-sort(unique(c(badBatches,as.character(batchIndex)[which(batchwiseSampleCounts < minimumBatchSize)])))
  cat(paste0("WARNING: Removing ",length(badBatches)," batches with <",minimumBatchSize," samples and/or no designated control (GIS) replicates found:\n",paste(badBatches,collapse=", "),"\n"))
  for (batch in badBatches) {
    traits<-traits[-which(sampleIndex$batch==batch),]
    sampleIndex<-sampleIndex[-which(sampleIndex$batch==batch),]
    cleanDat<-cleanDat[,match(rownames(traits),colnames(cleanDat))]
    batchIndex<-as.character(unique(sampleIndex$batch))
  }
}


## Specify Normalization 'channels' as batch-specific indices, for mean or median initial denominator
   # *Now allows setting GISindices even if your GIS is on different channels in different batches, essentially if channel order is not the same in each batch.
   # This enables using TAMPOR for batched LFQ or LFQ of different cohorts.
if ("GIS" %in% colnames(traits)) GISchannels<-unique(c(GISchannels,rownames(traits)[which(traits$GIS=="GIS")])) #fixes bug where GIS channels specified in traits$GIS not recognized after splitting channel names when batchPrefixInSampleNames = TRUE
GISindices<-list()
i.prev=0
offset=0
iter=0
unique_batches <- unique(sampleIndex$batch)
for ( i in as.character(batchIndex) ) {
  GISindices[[i]]<-vector()
  iter=iter+1
#  for (denomChannel in GISchannels) {
#    GISindices[[i]] <- c(GISindices[[i]], which(sampleIndex$batch == unique(sampleIndex$batch)[iter] & (sampleIndex$channel == denomChannel | rownames(traits) == denomChannel | sampleIndex$channel == "GIS") ))
#  }
# Above code can hang for large batched sets in the 10s of 1000s; Alternate:  #***
  # Precompute unique batches and current batch
  current_batch <- unique_batches[iter]
  
  # Precompute row indices for the current batch
  batch_indices <- which(sampleIndex$batch == current_batch)
  
  # Precompute channel matches
  channel_matches <- sort(unique( 
                                  c( na.omit(match(GISchannels,sampleIndex$channel[batch_indices])),
                                     na.omit(which(rownames(traits)[batch_indices] %in% GISchannels)),
                                     na.omit(which(sampleIndex$channel[batch_indices] == "GIS"))
                                   )
                                )  )
  
  # Indices within context of full data set for this batch meeting criteria
  GISindices[[i]] <- batch_indices[channel_matches]

  if (!iter==1) { offset=offset+length(which(sampleIndex$batch==i.prev)) }
  GISindices[[i]] <- unique(GISindices[[i]]) - offset
  i.prev=i
}
#names(GISindices)<-batchIndex  #in case batchIndexes are numeric


## Finalize traits$GIS column and throw error message if one or more batch has no GIS
if (!sum(table(unlist(lapply(GISindices,length))))==length(batchIndex)) stop("ERROR: ONE OR MORE BATCH IS MISSING DENOMINATOR SAMPLE DESIGNATIONS.")

globalExperimentGISsamples<-as.vector(unlist(sapply(names(GISindices),function(x) rownames(traits)[which(traits$Batch==x)[ GISindices[[x]] ]])))
traits$GIS<-rep(NA,nrow(traits))
traits$GIS[match(globalExperimentGISsamples,rownames(traits))]<-"GIS"

if (!identical(which(traits$GIS=="GIS"),match(globalExperimentGISsamples,rownames(traits)))) stop("ERROR: GIS samples now specified in traits do not match those being used.")


# Check for 0 and negative values, and tell the user to address this in their data themselves, or choose to replace with NA
badValues=as.vector(na.omit(cleanDat[cleanDat<=0]))
badRows=rownames(cleanDat[which(apply(cleanDat,1,function(x) min(c(x,1),na.rm=TRUE))<=0),])  #c(x,1) suppresses warnings for data with rows all NA
if (length(badValues)>0) {
  it=0
  cat("ABUNDANCES FORMAT WARNING/CHOICE:  you have ",length(badValues)," values <=0 in your input Abundances across ",length(badRows)," rows.\n\nNOTE: First bad rows (up to 10):\n")
  while(it<=min(10,length(badRows))) {
    cat(badRows[it+1])
    cat(paste(cleanDat[badRows,][(it*ncol(cleanDat)+1):((it+1)*ncol(cleanDat))],collapse=", "),"\n")
    it=it+1
  } 
  fix <- readline(prompt=paste0("What do you want to do? [ENTER]=Exit function to fix; OR type any value to set these to NA: "))
  if(!fix=="") stop("\n")
  #replace with <=0 values with NA; otherwise +/- Inf values propagate later!
  cleanDat<-apply(cleanDat,2,function(x) {
    x[x<=0] <- NA
    x
  })

  cat("__________________________________\n")
  it=0
  cat("First fixed rows (up to 10):\n")
  while(it<=min(10,length(badRows))) {
    cat(badRows[it+1])
    cat(paste(cleanDat[badRows,][(it*ncol(cleanDat)+1):((it+1)*ncol(cleanDat))],collapse=", "),"\n")
    it=it+1
  } 
  # Show user rows that had 0 or negative values fixed (to NA).
  cat(paste0(length(badValues)," bad values (<=0) in matrix replaced with NA.\nFirst fixed rows (up to 100) with at least one bad value are listed here:\n"))
  if(!is.null(badRows)) head(badRows,100)
}


if(!sampleMedianRows[1]=='ALL') {
  if(is.numeric(sampleMedianRows)) {
    sampleMedianRows.idx=sampleMedianRows
    if(min(sampleMedianRows.idx)<1) stop("\n ERROR: step 2 sampleMedianRows include a value of 0 or less. All values of sampleMedianRows vector must specify rows that exist in input abundance matrix.\n")
    if(max(sampleMedianRows.idx)>nrow(cleanDat)) stop("\n ERROR: step 2 sampleMedianRows specify a value greater than then number of rows remaining. All values of sampleMedianRows vector must specify rows that exist in input abundance matrix.\n")
  } else {
    sampleMedianRows.idx=as.vector(na.omit(match(sampleMedianRows,rownames(cleanDat))))
  }
  sampleMedianRows.found<-rownames(cleanDat)[sampleMedianRows.idx]
  
  cat(paste0("\n sampleMedianRows requested: ",length(sampleMedianRows),"\n sampleMedianRows found: ",length(sampleMedianRows.found),"\n - These rows will be used for determining log2(ratio/central tendency) median in each sample for step 2 subtraction, sample-wise.\n\n"))
}


cleanDat.original <- cleanDat 


#####################################################################################
iterationTrackingDF <- data.frame(Iteration = 1:iterations, FrobeniusNorm = NA, FrobenPrev = NA, FrobenDiff = NA, FrobenOverFirstFroben = NA)
cat(paste0("Starting # of Rows in data: ", nrow(cleanDat), ".\n"))
#####################################################################################
for (repeats in 1:iterations) {
  if (repeats==1 & exists("cleanDatNorm2")) rm(cleanDatNorm2)
  timer.start <- Sys.time()
  ###################################################################################################################################
  # STEP 1a. Ratio data and prepare to row-normalize
  ratioedBatches <- batchGISavgs <- list()
  ratioCleanDatUnnorm <- data.frame(row.names = rownames(cleanDat))
  withinBatchGISgeomeans <- withinBatchRowGeomeans <- data.frame(row.names = rownames(cleanDat))

#  comb <- function(x, ...) lapply(x, function(i) do.call(list,lapply(x, function(y) x[[y]])))
#  step1a <- foreach(batch=batchIndex, .combine='comb', .multicombine=TRUE, .init=list(list(), list(), list())) %dopar% {
  step1a <- foreach(batch=as.character(batchIndex)) %dopar% {
    tempForAvg <- matrix()
    tempForAvg <- as.data.frame(as.matrix(cleanDat[, which(sampleIndex$batch == batch)][, GISindices[[batch]] ], nrow = nrow(cleanDat), ncol = dim(cleanDat[, which(sampleIndex$batch == batch)][, GISindices[[batch]] ])[2]))
    batchGISavgs <- apply(tempForAvg, 1, function(x) eval(parse(text = paste0(meanOrMedian, "(x,na.rm=TRUE)")))) # ADDED na.rm v04 ##MEAN/MEDIAN FUNCTION CHOICE***
    ratioedBatches <- cleanDat[, which(sampleIndex$batch == batch)] / batchGISavgs

    ## Below unnormed ratio data are only assembled for graphing purposes, for comparison to step 1b and final step 2 output
    ## If batches are randomized channels distributing cases and controls evenly across all batches, useAllNonGIS==TRUE
    if (useAllNonGIS) {
      df3 <- as.data.frame(as.matrix(apply(ratioedBatches[, -GISindices[[batch]] ], 1, function(x) eval(parse(text = paste0("2^", meanOrMedian, "(log2(na.omit(x)))")))), ncol = dim(ratioedBatches)[2], nrow = dim(ratioedBatches)[1])) ## MEAN/MEDIAN FUNCTION CHOICE***
                                                                                                                   #as.matrix(), NOT matrix()
    } else {
      ## If we cannot rely on the robust assumption of batch-to-batch biological equivalence (with randomized sample order across all avalable channels in all batches), then use robust mean of GIS samples only
      df3  <- as.data.frame(as.matrix(apply(ratioedBatches, 1, function(x) eval(parse(text = paste0("2^", meanOrMedian, "(log2(na.omit(x[GISindices[[batch]] ])))")))), ncol = dim(ratioedBatches)[2], nrow = dim(ratioedBatches)[1]))
                                                                                                                   #as.matrix(), NOT matrix()
    } ## MEAN/MEDIAN FUNCTION CHOICE***
  return(list(batchGISavgs, ratioedBatches, df3))
  }

  # re-combine list elements from three outputs, over all batches
  batchGISavgs            = do.call(list,lapply(step1a,function(x){x[[1]]}))
  ratioedBatches          = do.call(list,lapply(step1a,function(x){x[[2]]}))
  withinBatchRowGeomeans  = as.data.frame(do.call(cbind,lapply(step1a,function(x){x[[3]]}))) #was also set to "withinBatchGISgeomeans" (not used below)
  names(batchGISavgs) <- names(ratioedBatches) <- colnames(withinBatchRowGeomeans)  <- batchIndex

  ratioCleanDatUnnorm <- do.call(cbind,ratioedBatches)


  ###################################################################################################################################
  ## Step 1b. Complete row-normalization. This step normalizes rows within batch by batchCorrFactors; ratioCleanDatUnnorm does not go through this step
  meanBatchGeomeans <- apply(withinBatchRowGeomeans, 1, function(x) eval(parse(text = paste0(meanOrMedian, "(na.omit(x))")))) ## MEAN/MEDIAN FUNCTION CHOICE***

  ## Rowwise (RW) relative abundances from GIS (or representative samples),
  # if we want to take the whole protein row (across batches) back to abundance after normalization step2 is complete**
  RW.relAbunFactors <- RW.GISavgs <- data.frame(row.names = rownames(cleanDat))
  RW.GISavgs <- cbind(apply(data.frame(column = as.character(batchIndex)), 1, function(x) batchGISavgs[[x]]))
  colnames(RW.GISavgs) <- as.character(batchIndex)
  RW.relAbunFactors <- apply(RW.GISavgs, 1, function(x) eval(parse(text = paste0("2^", meanOrMedian, "(log2(na.omit(x)))")))) #** relative abundance multipliers for recovery of relative abundance (all rows from input cleanDat) ##MEAN/MEDIAN FUNCTION CHOICE***
  rownames(RW.GISavgs) <- names(RW.relAbunFactors) <- rownames(cleanDat)

  ## Calculate Step 1b multipliers to complete RW normalization
  batchCorrFactors <- meanBatchGeomeans / withinBatchRowGeomeans
  colnames(batchCorrFactors) <- colnames(withinBatchRowGeomeans) <- as.character(batchIndex)

  # rowwise correct by multiplers (batchCorrFactors)
  ratioCleanDatNorm <- foreach(batch=as.character(batchIndex), .combine='cbind', .multicombine=TRUE) %dopar% as.data.frame(ratioedBatches[[batch]] * batchCorrFactors[, batch])

  ###################################################################################################################################

  ## log2 transform ratios output from step 1a, and 1b
  cleanDat.log2.ratioUnnorm <- log2(ratioCleanDatUnnorm)
  cleanDatNormNoColScaling <- log2(ratioCleanDatNorm)


  ## For cleanDat.log2.ratioUnnorm:
  ## Enforce <50% missingness (1 less than half of columns (or round down half if odd number of columns))
  LThalfSamples <- length(colnames(cleanDat.log2.ratioUnnorm)) * fractionNAmax
  LThalfSamples <- LThalfSamples - if ((length(colnames(cleanDat.log2.ratioUnnorm)) %% 2) == 1) {
    0.5
  } else {
    1.0
  }

  removedRownames1 <- rownames(cleanDat.log2.ratioUnnorm[which(rowSums(as.matrix(is.na(cleanDat.log2.ratioUnnorm))) > LThalfSamples), ]) # list rows to be removed
  removedRownames1

  # remove rows with >=50% missing values (only if there are some rows to be removed)
  if (length(na.omit(match(removedRownames1, rownames(cleanDat.log2.ratioUnnorm)))) == length(removedRownames1) & length(removedRownames1) > 0) {
    cleanDat.log2.ratioUnnorm <- cleanDat.log2.ratioUnnorm[-match(removedRownames1, rownames(cleanDat.log2.ratioUnnorm)), ]
  } else {
    cat("")
  } # "no rows removed.\n"); }
  dim(cleanDat.log2.ratioUnnorm)


  ## For cleanDatNormNoColScaling and companion relative abundance factors:
  ## Enforce <50% missingness (1 less than half of cleanDatNormNoColScaling columns (or round down half if odd number of columns))
  LThalfSamples <- length(colnames(cleanDatNormNoColScaling)) * fractionNAmax
  LThalfSamples <- LThalfSamples - if ((length(colnames(cleanDatNormNoColScaling)) %% 2) == 1) {
    0.5
  } else {
    1.0
  }

  removedRownames <- rownames(cleanDatNormNoColScaling[which(rowSums(as.matrix(is.na(cleanDatNormNoColScaling))) > LThalfSamples), ]) # list rows to be removed
#  removedRownames

  # remove rows with >=50% missing values (only if there are some rows to be removed)
  if (length(as.vector(na.omit(match(removedRownames, rownames(cleanDatNormNoColScaling))))) == length(removedRownames) & length(removedRownames) > 0) {
    cleanDatNormNoColScaling <- cleanDatNormNoColScaling[-match(removedRownames, rownames(cleanDatNormNoColScaling)), ]
    RW.relAbunFactors.HiMissRmvd <- RW.relAbunFactors[-match(removedRownames, names(RW.relAbunFactors))]
    #nrow(cleanDatNormNoColScaling)-length(removedRownames1) #, on iteration 1 is the final row count for the matrix; for this CSF data it should be 2875.
    cat(paste0("[iter_", repeats, "] Removed ", length(removedRownames), " high missingness rows. ", nrow(cleanDatNormNoColScaling), " rows remaining."))
  } else {
    cat(paste0("[iter_", repeats, "] No rows removed with >=50% missing values."))
    RW.relAbunFactors.HiMissRmvd <- RW.relAbunFactors
  }
  dim(cleanDatNormNoColScaling)


  if(repeats==1) {
    ratioCleanDatUnnorm.iter1 <- 2^cleanDat.log2.ratioUnnorm
    relAbundanceUnnorm.iter1  <- as.matrix(RW.relAbunFactors.HiMissRmvd*ratioCleanDatUnnorm.iter1)
    colnames(ratioCleanDatUnnorm.iter1)<-colnames(relAbundanceUnnorm.iter1)<- colnames(ratioCleanDatNorm)
  }


  prevIterCleanDatNorm2 <- data.frame(matrix(0, nrow = nrow(cleanDatNormNoColScaling), ncol = ncol(cleanDatNormNoColScaling)))
  if (exists("cleanDatNorm2")) prevIterCleanDatNorm2 <- cleanDatNorm2


  ###################################################################################################################################
  ## Step 2, Enforce equal loading assumption on output of step 1b (all well-quantified proteins equally considered/weighted)
  colMeans(ratioCleanDatUnnorm, na.rm = TRUE)
  # should be zero, but they are usually not.

  ## Set all column means to ~0 (10^-16 or less), essentially 0=log2(ratio/GIS) for all column means or an average ratio/GIS=1
#  cleanDatNorm <- scale(cleanDatNormNoColScaling, scale = FALSE)
#  colMeans(cleanDatNorm, na.rm = TRUE)

  if(sampleMedianRows[1]=="ALL") {
    ## alternative columnwise normalization operation using mean or median [equivalent to scale() function above, if colAvg=mean(x,na.rm=TRUE) ]
    cleanDatNorm2 <- apply(cleanDatNormNoColScaling, 2, function(x) {
      colAvg <- eval(parse(text = paste0(meanOrMedian, "(x,na.rm=TRUE)"))) ## MEAN/MEDIAN FUNCTION CHOICE***
      outputCol <- x - colAvg # rep(colAvg,length(x));
      outputCol
    })
  } else {
    sampleMedianRows.idx=as.vector(na.omit(match(sampleMedianRows.found,rownames(cleanDatNormNoColScaling))))
    if(!length(sampleMedianRows.idx)>5) stop("\n ERROR: step 2 requires at least 5 sampleMedianRows specified that survived removal due to 50%+ missingness rowwise. Specify different row (indices) for parameter sampleMedianRows, or leave as 'ALL'.\n\n")
    cat(paste0(" sampleMedianRows: ",length(sampleMedianRows.idx),"/",length(sampleMedianRows)," (remaining/requested)  "))

    ## alternative columnwise normalization operation using mean or median [equivalent to scale() function above, if colAvg=mean(x,na.rm=TRUE) ]
    cleanDatNorm2 <- apply(cleanDatNormNoColScaling, 2, function(x) {
      colAvg <- eval(parse(text = paste0(meanOrMedian, "(x[sampleMedianRows.idx],na.rm=TRUE)"))) ## MEAN/MEDIAN FUNCTION CHOICE***
      outputCol <- x - colAvg # rep(colAvg,length(x));
      outputCol
    })
  }

  ## show columnwise normalization/scaling methods are equivalent (with rounding to a few decimals, at least)
#  colMeans(cleanDatNorm, na.rm = TRUE) == colMeans(cleanDatNorm2, na.rm = TRUE) # TRUE if meanOrMedian="mean"

  ## GIS-derived relative RW abundances can be applied back to the Step2 RW & CW-normalized data's rows.
  ## values of cleanDatNorm2 are log2(measurement/batch-stabilized GIS abundance).
  if (paste(names(RW.relAbunFactors.HiMissRmvd), collapse = ",") == paste(rownames(cleanDatNorm2), collapse = ",")) {
    relAbundanceNorm2 <- RW.relAbunFactors.HiMissRmvd * 2^cleanDatNorm2
  } else {
    stop(paste0("\n[iter_", repeats, "] ERROR: step 2 data with removed high missingness rows does not match relative abundance factors after trying to remove same rows.\n"))
  }


  ###################################################################################################################################
  ## Prepare for convergence check and next iteration

  cleanDat <- relAbundanceNorm2
  DFforFrobCurrent <- apply(cleanDatNorm2, 2, as.numeric)
  DFforFrobPrev <- apply(prevIterCleanDatNorm2, 2, as.numeric)
  removeColumnsCurrent=which(apply(DFforFrobCurrent, 2, function(x) sum(is.na(x))) == nrow(DFforFrobCurrent))
  removeColumnsPrev=which(apply(DFforFrobPrev, 2, function(x) sum(is.na(x))) == nrow(DFforFrobPrev))

  if(!all(rowSums(is.na(DFforFrobCurrent)) >=1)) {  # norm function will work if a row or more have no missing values.
    if (length(removeColumnsCurrent)>0 ) {
      frobeniusNormCurrent <- norm(na.omit(DFforFrobCurrent[,-removeColumnsCurrent]), type = "F")
    } else {
      frobeniusNormCurrent <- norm(matrix((na.omit(DFforFrobCurrent))), type = "F")
    }
    if (length(removeColumnsPrev)>0 ) {
      frobeniusNormPrev <- norm(na.omit(DFforFrobPrev[,-removeColumnsPrev]), type = "F")
    } else {
      frobeniusNormPrev <- norm(na.omit(DFforFrobPrev), type = "F")
    }
  } else {  # Calculate a norm ignoring missing values without throwing out rows.
    if (length(removeColumnsCurrent)>0 ) {
      frobeniusNormCurrent <- altFrobenius_norm(DFforFrobCurrent[,-removeColumnsCurrent])
    } else {
      frobeniusNormCurrent <- altFrobenius_norm(matrix(DFforFrobCurrent))
    }
    if (length(removeColumnsPrev)>0 ) {
      frobeniusNormPrev <- altFrobenius_norm(DFforFrobPrev[,-removeColumnsPrev])
    } else {
      frobeniusNormPrev <- altFrobenius_norm(DFforFrobPrev)
    }
  }

  initialFrobeniusNorm <- if (repeats == 1) {
      frobeniusNormCurrent
  } else {
    initialFrobeniusNorm
  }
  iterationTrackingDF[repeats, ] <- c(repeats, frobeniusNormCurrent, frobeniusNormPrev, frobeniusNormPrev - frobeniusNormCurrent, frobeniusNormCurrent / initialFrobeniusNorm)
  timer <- difftime(Sys.time(), timer.start, units="secs")
  time.sec<-round(as.numeric(timer),2)
  cat(paste0("     ",time.sec," sec     iteration convergence tracking (Frobenius Norm Difference):  ", signif(iterationTrackingDF[repeats, 4], 3), "\n"))
  if (repeats>1 & abs(iterationTrackingDF[repeats,4])<0.00000001) { cat("...Reached convergence criterion (Frobenius Norm Difference)<1e-8!\n"); break; }
} # closes "for (repeats in 1:iterations)"
iterations.intended=iterations
iterations=repeats
converged=as.logical(abs(iterationTrackingDF[repeats,4])<0.00000001)

## Output plots showing convergence approach over the iterations run.
if (iterations>1) {
  par(mfrow=c(1,2),oma=c(0,0,2,0))
  plot(iterationTrackingDF[2:nrow(iterationTrackingDF),1],abs(iterationTrackingDF[2:nrow(iterationTrackingDF),4]),xlab="Iteration",ylab="Frobenius Norm Diff from Previous", main="Linear Scale")
  logAbsFrobenDiff<-log10(abs(iterationTrackingDF[2:nrow(iterationTrackingDF),4]))
  logAbsFrobenDiff[logAbsFrobenDiff< -10] <- -10    #censor log10(0) and all norm differences of less than 1e-10 in next plot
  plot(iterationTrackingDF[2:nrow(iterationTrackingDF),1],logAbsFrobenDiff,xlab="Iteration",ylab="log10(|Frobenius Norm Diff from Previous|)", main="Log10 Scale")
  mtext(paste0("Convergence Tracking, Iterations 2-",iterations),line=0.25,outer=TRUE,cex=1.5)
  convergencePlots<-recordPlot()
  PDFpage3=TRUE
} else {
#  convergencePlot=FALSE
#  logConvergencePlot=FALSE
  convergencePlots=FALSE
  PDFpage3=FALSE
}

rm(step1a)

# END NORMALIZATION
###################################################################################################################################

## Below unnormed ratio (from iteration 1) data are only assembled for graphing purposes, for comparison to step 1b and final step 2 output
ratioCleanDatUnnorm <- ratioCleanDatUnnorm.iter1 #as.matrix(cbind(ratioCleanDatUnnorm, ratioedBatches[[batch]]))
relAbundanceUnnorm  <- relAbundanceUnnorm.iter1  #as.matrix(RW.relAbunFactors.HiMissRmvd*ratioCleanDatUnnorm)


#saveRDS(cleanDatNorm2,paste0(outputtabs,"cleanDat.TAMPOR_",iterations,"iter.RDS"))
#saveRDS(relAbundanceNorm2,paste0(outputtabs,"relAbundanceNorm2.alternateTAMPOR.output_",iterations,"iter.RDS"))
#saveRDS(ratioCleanDatUnnorm,paste0(outputtabs,"ratioCleanDatUnnorm.naiveRatioOverGIS_",iterations,"iter.RDS"))
#saveRDS(relAbundanceUnnorm,paste0(outputtabs,"relAbundanceUnnorm.naiveRelAbun.AfterRatioOverGIS_",iterations,"iter.RDS"))

###################################################################################################################################
##Prepare sample traits

#library(WGCNA) #only used for labels2colors() function; fn reproduced instead of requiring install of all WGCNA dependencies.

## Insure that traits column order are correct.
traits<-traits[match(colnames(cleanDatNorm2),rownames(traits)),]

traits$BatchColor<-labels2colors(traits$Batch)

numericMeta<-traitsWithGIS<-traits

###################################################################################################################################
## Sample Removal and finalization of cleanDat as log2(ratio).

#cleanDat<-cleanDatNorm2<-readRDS(paste0(outputtabs,"cleanDat.TAMPOR_",iterations,"iter.RDS"))
#relAbundanceNorm2<-readRDS(paste0(outputtabs,"relAbundanceNorm2.alternateTAMPOR.output_",iterations,"iter.RDS"))
#ratioCleanDatUnnorm<-readRDS(paste0(outputtabs,"ratioCleanDatUnnorm.naiveRatioOverGIS_",iterations,"iter.RDS"))
#relAbundanceUnnorm<-readRDS(paste0(outputtabs,"relAbundanceUnnorm.naiveRelAbun.AfterRatioOverGIS_",iterations,"iter.RDS"))

  ## One can remove GIS here, depending on whether these are (not) biologically meaningful for comparisons downstream
  if (length(which(traits$GIS == "GIS"))>0) if(removeGISafter) numericMeta <- traits <- traits[-which(traits$GIS == "GIS"), ] #***
  cleanDat <- cleanDatNorm2[, match(rownames(traits), colnames(cleanDatNorm2))] #NOTE: overwrites cleanDat (intermediate structure during cleanup)

  ratioCleanDatUnnorm <- ratioCleanDatUnnorm[, match(rownames(traits), colnames(ratioCleanDatUnnorm))]
  relAbundanceNorm2 <- relAbundanceNorm2[, match(rownames(traits), colnames(relAbundanceNorm2))]


  # remove any ignored columns (all NA)
  removeColumns.cleanDat<-samplesToIgnore[as.vector(na.omit(match(colnames(cleanDat),samplesToIgnore)))] #which(apply(cleanDat, 2, function(x) sum(is.na(x))) == nrow(cleanDat))
  if (length(removeColumns.cleanDat)>0 ) {
    cleanDat<-cleanDat[,-match(removeColumns.cleanDat,colnames(cleanDat))]
    numericMeta<-traits<-traits[match(colnames(cleanDat),rownames(traits)),]

    ratioCleanDatUnnorm <- ratioCleanDatUnnorm[, -match(removeColumns.cleanDat,colnames(ratioCleanDatUnnorm))]
    relAbundanceNorm2 <- relAbundanceNorm2[, -match(removeColumns.cleanDat,colnames(relAbundanceNorm2))]
  }
  cleanDat.orig<-as.data.frame(cleanDat.original)[, match(rownames(traits), colnames(cleanDat.original))]

  # Only keep rows in original cleanDat (cleanDat.orig) matching final normed data
#  dim(cleanDat.orig)
  cleanDat.orig<-cleanDat.orig[match(rownames(cleanDat),rownames(cleanDat.orig)),]
#  dim(cleanDat.orig)


###################################################################################################################################
## Check the data.

## Check the data integrity and generalize improvements

file1 <- paste0(outputfigs,"/TAMPOR-Improvement.Vis(", iterations, "iterations)-", nrow(cleanDat), "x", ncol(cleanDat), "_", outputSuffix, ".pdf")
cat(paste0("Writing PDF output to ",file1,"\n"))

# Has variance been stabilized?
force.inf2NA <- function(x) ifelse(is.finite(x), x, NA)
plot1 <- vsn::meanSdPlot(force.inf2NA(log2(as.matrix(relAbundanceNorm2))),xlab="rank(mean):  log2(relAbundanceNorm2)")
#plot2 <- vsn::meanSdPlot(force.inf2NA(cleanDat),xlab="rank(mean):  Final TAMPOR log2(ratio)") # Ratio data has u-shape.
#plot3 <- vsn::meanSdPlot(force.inf2NA(log2(relAbundanceUnnorm)),xlab="rank(mean):  log2(relabundanceUnnorm)") # Abundance. 
plot4 <- vsn::meanSdPlot(force.inf2NA(log2(as.matrix(ratioCleanDatUnnorm*RW.relAbunFactors.HiMissRmvd))),xlab="rank(mean):  Naive log2(abundance/GIS *rel. Abun.)") # ratio converted to rel abun does not have u-shape.
plot5 <- vsn::meanSdPlot(force.inf2NA(log2(as.matrix(cleanDat.orig))),xlab="rank(mean): log2(orig. abun)") #RAW

ymax=1.25*max(c(plot1$sd,plot4$sd,plot5$sd))  #(scale all the meanSD plots so they are comparable)
horizMin=min(c(plot1$sd,plot4$sd,plot5$sd))   #set to y=minimum of lowest-reaching trendline form mean-sd plots being printed

#meanSDplots2<-ggpubr::ggarrange(plotlist=list(plot5$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(paste0("Variance (mean-SD) Plots for Abundance\nMatrices of Dimensions ",nrow(relAbundanceNorm2)," x ",ncol(relAbundanceNorm2))), plot4$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n "), plot1$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n ") ),nrow=2,ncol=3,common.legend=TRUE) #frame fill order is columwise, but is rowwise for marrangeGrob. 


# Compare same plots without GIS, if we did not remove GIS after TAMPOR.
if(!length(which(traits$GIS=="GIS"))==nrow(traits) & !removeGISafter) {
  plot1.noGIS <- vsn::meanSdPlot(log2(as.matrix(relAbundanceNorm2))[,-which(traits$GIS=="GIS")],xlab="rank(mean):  log2(relAbundanceNorm2)")
  plot4.noGIS <- vsn::meanSdPlot(log2(as.matrix(ratioCleanDatUnnorm*RW.relAbunFactors.HiMissRmvd))[,-which(traits$GIS=="GIS")],xlab="rank(mean):  Naive log2(abundance/GIS *rel. Abun.)") # ratio converted to rel abun does not have u-shape.
  plot5.noGIS <- vsn::meanSdPlot(log2(as.matrix(cleanDat.orig))[,-which(traits$GIS=="GIS")],xlab="rank(mean): log2(orig. abun)") #RAW

  ymax.noGIS=1.25*max(c(plot1.noGIS$sd,plot4.noGIS$sd,plot5.noGIS$sd))
  horizMin.noGIS=min(c(plot1.noGIS$sd,plot4.noGIS$sd,plot5.noGIS$sd))

  # all 6 meanSDplots on one page:
  meanSDplots<-ggpubr::ggarrange(plotlist=list(plot5$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(paste0("Variance (mean-SD) Plots for Abundance\nMatrices of Dimensions ",nrow(relAbundanceNorm2)," x ",ncol(relAbundanceNorm2))), plot4$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n "), plot1$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n "),
                                               plot5.noGIS$gg + scale_y_continuous(limits=c(0,ymax.noGIS)) + geom_hline(yintercept=horizMin.noGIS, linetype="dashed", color = "yellow", size=1.2) + ggtitle(paste0("Variance (mean-SD) Plots for Abundance\nMatrices* of Dimensions ",nrow(relAbundanceNorm2)," x ",ncol(relAbundanceNorm2[,-which(traits$GIS=="GIS")]))), plot4.noGIS$gg + scale_y_continuous(limits=c(0,ymax.noGIS)) + geom_hline(yintercept=horizMin.noGIS, linetype="dashed", color = "yellow", size=1.2) + ggtitle("*No GIS control samples included\n "), plot1.noGIS$gg + scale_y_continuous(limits=c(0,ymax.noGIS)) + geom_hline(yintercept=horizMin.noGIS, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n ") ),nrow=2,ncol=3,common.legend=TRUE) #frame fill order is columwise, but is rowwise for marrangeGrob. 
} else {
  #just 3 meanSDplots:
  meanSDplots<-ggpubr::ggarrange(plotlist=list(plot5$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(paste0("Variance (mean-SD) Plots for Abundance\nMatrices of Dimensions ",nrow(relAbundanceNorm2)," x ",ncol(relAbundanceNorm2))), plot4$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n "), plot1$gg + scale_y_continuous(limits=c(0,ymax)) + geom_hline(yintercept=horizMin, linetype="dashed", color = "yellow", size=1.2) + ggtitle(" \n ") ),nrow=2,ncol=3,common.legend=TRUE)
}

# Output page 1 for capture
par(mfrow=c(2,3))
print(meanSDplots)
meanSDplots.rec<-recordPlot()


# Generate MDS plot pre/post normalization.
if(!skipMDS) {
  cat("Generating MDS plots in PDF output...\n")
} else {
  cat("Skipping MDS plotting [skipMDS=TRUE]...\n")
}

MDSplotter.fallback<-function(data,colors,title) {
	skip=FALSE
	x <- as.matrix(data)
	nsamples <- ncol(x)
	if(nsamples < 3) {
	  cat(paste0("x Only ",nsamples," columns of data: need at least 3. Skipping MDS plot.\n"))
	  skip=TRUE
	}
	# Remove rows with missing or Inf values
	bad <- rowSums(is.finite(x)) < nsamples
	if(any(bad)) x <- x[!bad,,drop=FALSE]
	nprobes <- nrow(x)
	if(nsamples < 2) {
	  cat("x Number of samples is less than number of dimensions (2). Skipping MDS plot.\n")
	  skip=TRUE
	}
	#cat(paste0("nprobes: ",nprobes,"\nndim: ",ndim,"\n\n"))
	if(nprobes < 2) {
	  cat("x Number of rows of data with no missing values across all samples is less than 2. Skipping MDS plot.\n")
	  skip=TRUE
	}
	if(!skip) { 
	  limma::plotMDS(data,col=colors,main=title)
	} else {
	  frame()
	  mtext(paste0(title,"\n\n\n\ninput not suitable for MDS\n\n-See console messages."), adj=0.5, line=-4, cex=0.85, font=2);
	}
}

if(!skipMDS) {
  
  par(mfrow=c(2.3,3))
  
  # The starkest comparisons -- from original norm abundance (no ratio), to Naive abun/GIS (converted back to rel abundance), to polished output
  MDS1.plot <- MDSplotter.fallback(data=log2(cleanDat.orig),colors=traits$BatchColor,title="ORIGINAL log2(abundance)")
  MDS2.plot <- MDSplotter.fallback(data=log2(ratioCleanDatUnnorm*RW.relAbunFactors.HiMissRmvd),colors=traits$BatchColor,title="Naive log2(abundance/GIS *Rel Abun.)")
  MDS3.plot <- MDSplotter.fallback(data=log2(relAbundanceNorm2),colors=traits$BatchColor,title=paste0("TAMPOR log2(abundance) [",iterations,"iterations]"))
  
  # Print same plots without GIS, if we did not remove GIS after TAMPOR.
  if(!length(which(traits$GIS=="GIS"))==nrow(traits) & !removeGISafter) {
    MDS1.noGIS.plot <- MDSplotter.fallback(data=log2(cleanDat.orig)[,-which(traits$GIS=="GIS")],colors=traits$BatchColor[-which(traits$GIS=="GIS")],title="ORIGINAL log2(abundance)")
    mtext("GIS Removed")
    MDS2.noGIS.plot <- MDSplotter.fallback(data=log2(ratioCleanDatUnnorm*RW.relAbunFactors.HiMissRmvd)[,-which(traits$GIS=="GIS")],colors=traits$BatchColor[-which(traits$GIS=="GIS")],title="Naive log2(abundance/GIS *Rel Abun.)")
    mtext("GIS Removed")
    MDS3.noGIS.plot <- MDSplotter.fallback(data=log2(relAbundanceNorm2)[,-which(traits$GIS=="GIS")],colors=traits$BatchColor[-which(traits$GIS=="GIS")],title=paste0("TAMPOR log2(abundance) [",iterations,"iterations]"))
    mtext("GIS Removed")
  }
  
  MDSplots.rec<-recordPlot()
}


#Output PDF with same pages, and convergencePlots on page 3
pdf(file1, width = 18, height = 12)

#par(mfrow=c(2,3))
print(meanSDplots.rec) #page 1

#par(mfrow=c(2.3,3))
if (!skipMDS) print(MDSplots.rec) #page 2

# Page 3 (convergence tracking)
if (PDFpage3) print(convergencePlots)

dev.off() #closes file1


# Output list of all data and visualizations

#Auto notes output may be a good idea-- X (rows/proteins/transcripts/metabolites) remain, Y iterations (characteristic of convergence). 
#mean-SD plots have improved mean and minimum SD, going from no norm, to naive ratio*rel abun, to TAMPOR rel abundance...

if(!length(which(traits$GIS=="GIS"))==nrow(traits) & !removeGISafter) { #single return statement now ok if only outputting recorded PDF pages...
  if(!skipMDS) {
    return(list(cleanDat=cleanDatNorm2,cleanRelAbun=relAbundanceNorm2,traits=traits,cleanDat.oneIter=ratioCleanDatUnnorm,cleanRelAbun.oneIter=relAbundanceUnnorm,
                convergencePlots=convergencePlots,meanSDplots=meanSDplots.rec,MDSplots=MDSplots.rec,
  #              convergencePlot=convergencePlot,logConvergencePlot=logConvergencePlot,varPlot.input=plot5,varPlot.oneIter=plot4,varPlot.cleanRelAbun=plot1,MDSplot.input=MDS1.plot,MDSplot.oneIter=MDS2.plot,MDSplot.cleanRelAbun=MDS3.plot,
  #              varPlot.input.noGIS=plot5.noGIS,varPlot.oneIter.noGIS=plot4.noGIS,varPlot.cleanRelAbun.noGIS=plot1.noGIS,MDSplot.input.noGIS=MDS1.noGIS.plot,MDSplot.oneIter.noGIS=MDS2.noGIS.plot,MDSplot.cleanRelAbun.noGIS=MDS3.noGIS.plot,
                iterations=iterations,converged=converged))
               } else {  # skipMDS=TRUE
    return(list(cleanDat=cleanDatNorm2,cleanRelAbun=relAbundanceNorm2,traits=traits,cleanDat.oneIter=ratioCleanDatUnnorm,cleanRelAbun.oneIter=relAbundanceUnnorm,
                convergencePlots=convergencePlots,meanSDplots=meanSDplots.rec, # MDSplots=MDSplots.rec,
                iterations=iterations,converged=converged))
               }
} else {
  if(!skipMDS) {
    return(list(cleanDat=cleanDatNorm2,cleanRelAbun=relAbundanceNorm2,traits=traits,cleanDat.oneIter=ratioCleanDatUnnorm,cleanRelAbun.oneIter=relAbundanceUnnorm,
                convergencePlots=convergencePlots,meanSDplots=meanSDplots.rec,MDSplots=MDSplots.rec,
  #              convergencePlot=convergencePlot,logConvergencePlot=logConvergencePlot,varPlot.input=plot5,varPlot.oneIter=plot4,varPlot.cleanRelAbun=plot1,MDSplot.input=MDS1.plot,MDSplot.oneIter=MDS2.plot,MDSplot.cleanRelAbun=MDS3.plot,
                iterations=iterations,converged=converged))
               } else {  # skipMDS=TRUE
    return(list(cleanDat=cleanDatNorm2,cleanRelAbun=relAbundanceNorm2,traits=traits,cleanDat.oneIter=ratioCleanDatUnnorm,cleanRelAbun.oneIter=relAbundanceUnnorm,
                convergencePlots=convergencePlots,meanSDplots=meanSDplots.rec, # MDSplots=MDSplots.rec,
                iterations=iterations,converged=converged))
               }
}

} #close TAMPOR function


altFrobenius_norm <- function(mat) {
  sqrt(sum(mat^2, na.rm = TRUE))
}



labels2colors <- function (labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey", commonColorCode = TRUE)  #small function from WGCNA
{
    if (is.null(colorSeq)) colorSeq = 
            c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue",
            "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen",
            "skyblue3", "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue", "plum2", "thistle2", "thistle1", "salmon4",
            "palevioletred3", "navajowhite2", "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4", "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2", "yellow4", "skyblue1", "plum", "orangered3",
            "mediumpurple2", "lightsteelblue", "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", "brown2", "blue2", "darkviolet", "plum3", "thistle3", "thistle", "salmon2", "palevioletred2", "navajowhite1", "magenta4",
            "lightpink3", "lavenderblush2", "honeydew", "darkseagreen3", "coral", "antiquewhite2", "coral3", "mediumpurple4", "skyblue4", "yellow3", "sienna4", "pink4", "orangered1", "mediumpurple1", "lightslateblue",
            "lightblue4", "indianred3", "firebrick3", "darkolivegreen2", "blueviolet", "blue4", "deeppink", "plum4", "thistle4", "tan4", "salmon1", "palevioletred1", "navajowhite", "magenta3", "lightpink2", "lavenderblush1",
            "green4", "darkseagreen2", "chocolate4", "antiquewhite1", "coral4", "mistyrose", "slateblue", "yellow2", "sienna2", "pink3", "orangered", "mediumpurple", "lightskyblue4", "lightblue3", "indianred2", "firebrick2",
            "darkolivegreen1", "blue3", "brown1", "deeppink1", "powderblue", "tomato", "tan3", "royalblue3", "palevioletred", "moccasin", "magenta2", "lightpink1", "lavenderblush", "green3", "darkseagreen1", "chocolate3",
            "aliceblue", "cornflowerblue", "navajowhite3", "slateblue1", "whitesmoke", "sienna1", "pink2", "orange4", "mediumorchid4", "lightskyblue3", "lightblue2", "indianred1", "firebrick", "darkgoldenrod4", "blue1",
            "brown3", "deeppink2", "purple2", "tomato2", "tan2", "royalblue2", "paleturquoise4", "mistyrose4", "magenta1", "lightpink", "lavender", "green2", "darkseagreen", "chocolate2", "antiquewhite", "cornsilk",
            "navajowhite4", "slateblue2", "wheat3", "sienna", "pink1", "orange3", "mediumorchid3", "lightskyblue2", "lightblue1", "indianred", "dodgerblue4", "darkgoldenrod3", "blanchedalmond", "burlywood", "deepskyblue", "red1",
            "tomato4", "tan1", "rosybrown4", "paleturquoise3", "mistyrose3", "linen", "lightgoldenrodyellow", "khaki4", "green1", "darksalmon", "chocolate1", "antiquewhite3", "cornsilk2", "oldlace", "slateblue3", "wheat1",
            "seashell4", "peru", "orange2", "mediumorchid2", "lightskyblue1", "lightblue", "hotpink4", "dodgerblue3", "darkgoldenrod1", "bisque3", "burlywood1", "deepskyblue4", "red4", "turquoise2", "steelblue4", "rosybrown3",
            "paleturquoise1", "mistyrose2", "limegreen", "lightgoldenrod4", "khaki3", "goldenrod4", "darkorchid4", "chocolate", "aquamarine", "cyan1", "orange1", "slateblue4", "violetred4", "seashell3", "peachpuff4",
            "olivedrab4", "mediumorchid1", "lightskyblue", "lemonchiffon4", "hotpink3", "dodgerblue1", "darkgoldenrod", "bisque2", "burlywood2", "dodgerblue2", "rosybrown2", "turquoise4", "steelblue3", "rosybrown1",
            "palegreen4", "mistyrose1", "lightyellow4", "lightgoldenrod3", "khaki2", "goldenrod3", "darkorchid3", "chartreuse4", "aquamarine1", "cyan4", "orangered2", "snow", "violetred2", "seashell2", "peachpuff3",
            "olivedrab3", "mediumblue", "lightseagreen", "lemonchiffon3", "hotpink2", "dodgerblue", "darkblue", "bisque1", "burlywood3", "firebrick1", "royalblue1", "violetred1", "steelblue1", "rosybrown", "palegreen3",
            "mintcream", "lightyellow3", "lightgoldenrod2", "khaki1", "goldenrod2", "darkorchid2", "chartreuse3", "aquamarine2", "darkcyan", "orchid", "snow2", "violetred", "seashell1", "peachpuff2", "olivedrab2",
            "mediumaquamarine", "lightsalmon4", "lemonchiffon2", "hotpink1", "deepskyblue3", "cyan3", "bisque", "burlywood4", "forestgreen", "royalblue4", "violetred3", "springgreen3", "red3", "palegreen1", "mediumvioletred",
            "lightyellow2", "lightgoldenrod1", "khaki", "goldenrod1", "darkorchid1", "chartreuse2", "aquamarine3", "darkgoldenrod2", "orchid1", "snow4", "turquoise3", "seashell", "peachpuff1", "olivedrab1", "maroon4",
            "lightsalmon3", "lemonchiffon1", "hotpink", "deepskyblue2", "cyan2", "beige", "cadetblue", "gainsboro", "salmon3", "wheat", "springgreen2", "red2", "palegreen", "mediumturquoise", "lightyellow1", "lightgoldenrod",
            "ivory4", "goldenrod", "darkorchid", "chartreuse1", "aquamarine4", "darkkhaki", "orchid3", "springgreen1", "turquoise1", "seagreen4", "peachpuff", "olivedrab", "maroon3", "lightsalmon2", "lemonchiffon", "honeydew4",
            "deepskyblue1", "cornsilk4", "azure4", "cadetblue1", "ghostwhite", "sandybrown", "wheat2", "springgreen", "purple4", "palegoldenrod", "mediumspringgreen", "lightsteelblue4", "lightcyan4", "ivory3", "gold3",
            "darkorange4", "chartreuse", "azure", "darkolivegreen3", "palegreen2", "springgreen4", "tomato3", "seagreen3", "papayawhip", "navyblue", "maroon2", "lightsalmon1", "lawngreen", "honeydew3", "deeppink4", "cornsilk3",
            "azure3", "cadetblue2", "gold", "seagreen", "wheat4", "snow3", "purple3", "orchid4", "mediumslateblue", "lightsteelblue3", "lightcyan3", "ivory2", "gold2", "darkorange3", "cadetblue4", "azure1", "darkorange1",
            "paleturquoise2", "steelblue2", "tomato1", "seagreen2", "palevioletred4", "navy", "maroon1", "lightsalmon", "lavenderblush4", "honeydew2", "deeppink3", "cornsilk1", "azure2", "cadetblue3", "gold4", "seagreen1",
            "yellow1", "snow1", "purple1", "orchid2", "mediumseagreen", "lightsteelblue2", "lightcyan2", "ivory1", "gold1") #WGCNA ordered standardColors()

    if (is.numeric(labels)) {
        if (zeroIsGrey) minLabel = 0
        else minLabel = 1
        if (any(labels < 0, na.rm = TRUE)) minLabel = min(c(labels), na.rm = TRUE)
        nLabels = labels
    }
    else {
        if (commonColorCode) {
            factors = factor(c(as.matrix(as.data.frame(labels))))
            nLabels = as.numeric(factors)
            dim(nLabels) = dim(labels)
        }
        else {
            labels = as.matrix(as.data.frame(labels))
            factors = list()
            for (c in 1:ncol(labels)) factors[[c]] = factor(labels[, c])
            nLabels = sapply(factors, as.numeric)
        }
    }
    if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
        nRepeats = as.integer((max(labels) - 1)/length(colorSeq)) + 1
        warning(paste0("Number of labels exceeds number of available colors. Some colors will be repeated ", nRepeats, " times."))
        extColorSeq = colorSeq
        for (rep in 1:nRepeats) extColorSeq = c(extColorSeq, paste(colorSeq, ".", rep, sep = ""))
    }
    else {
        nRepeats = 1
        extColorSeq = colorSeq
    }
    colors = rep("grey", length(nLabels))
    fin = !is.na(nLabels)
    colors[!fin] = naColor
    finLabels = nLabels[fin]
    colors[fin][finLabels != 0] = extColorSeq[finLabels[finLabels != 0]]
    if (!is.null(dim(labels))) dim(colors) = dim(labels)
    colors
}
