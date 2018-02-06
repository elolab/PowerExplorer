# calculate power for each sample size
# INPUT
# simData: the simulated data under null and alternative hypothesis
# alpha: desired false positive rate control
# saveResultData: whether to save the estimated power of each entry in *.RData file
# OUTPUT
# powerEst: a entry-named vector of power values (rounded as xx.xx)
# cutoffStatistics will be saved into as cutoffStats_[time].RData
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 19th, Dec, 2017
calPwr <- function(simData,
                           alpha=0.05,
                           dataType=c("RNASeq", "Proteomics"),
                           saveResultData=FALSE,
                           showOverallPower=TRUE) {
  # determine the dataType
  if(identical(dataType, "RNASeq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  } else {
    stop("dataType is not specified correctly.")
  }
  # extract statistics from simluated data
  statMatrix0 <- 
    do.call(cbind, lapply(simData, function(x) tail.matrix(x, 2)[1,]))
  statMatrix1 <- 
    do.call(cbind, lapply(simData, function(x) tail.matrix(x, 2)[2,]))
  
  nullData <- simData$nullData
  alterData <- simData$alterData
  entryNames <- rownames(statMatrix0) # IDs/identifiers for genes

  # Determine the cut-off statistics for each sample size and mean count
  # Based on the alpha value (default: 0.05)
  # control false positive rate as alpha
  cutoffStats <- apply(statMatrix0, 1,
                       function(x) quantile(x,
                                            1-alpha,
                                            na.rm=TRUE,
                                            type=9))

  # Estmate the power based on the rate of TP/NP
  powerEst <-
    vapply(entryNames, function(x)
                       mean(statMatrix1[x, ] > cutoffStats[x],
                       na.rm=TRUE), double(1))

  # save data
  names(powerEst) <- names(cutoffStats)
  if(saveResultData) {
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.power <- RDataName(ifelse(dataTypeSelect,
                                       "[RNASeq] Estimated Power",
                                       "[Proteomics] Estimated Power"))
    save(powerEst, file=paste0(getwd(), "/savedRData/", filename.power))
    message(paste0(">> Estimated power saved in directory savedRData."))
    message(paste0(">> Size: ",
                   round(file.size(paste0(getwd(),
                                          "/savedRData/",
                                          filename.power))/2^10, 2), " KB"))
  }
  if(showOverallPower) cat(paste0("\nOVERALL ESTIMATED POWER: ",
                                  round(mean(powerEst), 4), "\n\n"))
  powerEst
}
