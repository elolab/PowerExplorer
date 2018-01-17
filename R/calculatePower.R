# calculatePower
# INPUT
# simData: the simulated data under null and alternative hypothesis
# alpha: desired false positive rate control
# saveRData: whether to save the estimated power of each entry in *.RData file
# OUTPUT
# powerEst: a entry-named vector of power values (rounded as xx.xx)
# cutoffStatistics will be saved into as cutoffStats_[time].RData
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 19th, Dec, 2017
calculatePower <- function(simData,
                           alpha=0.05,
                           dataType=c("RNA-Seq", "Proteomics"),
                           saveRData=FALSE,
                           showOverallPower=TRUE) {
  # determine the dataType
  if(identical(dataType, "RNA-Seq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  } else {
    stop("dataType is not specified correctly.")
  }

  nullData <- simData$nullData
  alterData <- simData$alterData
  entryNames <- names(nullData$parameter0) # IDs/identifiers for genes

  # Determine the cut-off statistics for each sample size and mean count
  # Based on the alpha value (default: 0.05)

  if(showOverallPower) utils::timestamp()
  if(showOverallPower)
  cat(paste("Estimating cut-off statistics with false positive rate:", alpha))

  # control false positive rate as alpha
  cutoffStats <- apply(nullData$testStatistics, 1,
                       function(x) quantile(x,
                                            1-alpha,
                                            na.rm=TRUE,
                                            type=9))

  # Estmate the power based on the rate of TP/NP
  powerEst <-
    sapply(entryNames, function(x)
                       mean(alterData$testStatistics[x, ] > cutoffStats[x],
                       na.rm=TRUE))

  # save data
  names(powerEst) <- names(cutoffStats)
  if(saveRData == TRUE) {
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.power <- RDataName(ifelse(dataTypeSelect,
                                       "[RNA-Seq] Estimated Power",
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
