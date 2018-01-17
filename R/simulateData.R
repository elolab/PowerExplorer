# Simulation
# INPUT
# para0.list - a list of two mean vectors of genes in two case conditions
# para1.list - a list of two theta vectors of genes in two case conditions
# simNumRep - the number of replicates of each group for simulation
# ST - simulation times
# showProcess - TRUE: show the info of each simulation,
# FALSE: only the process percentage bar
# saveRData - whether to save the simulated data into *.RData file
# OUTPUT
# simulatedData - a list with two list objects: nullData and alterData both
# contains count data, Wald statistics, mean and dispersion of each gene
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 29th, Dec, 2017
#' @import utils
simulateData <- function(para0.list,
                         para1.list,
                         dataType=c("RNA-Seq", "Proteomics"),
                         simNumRep,
                         ST=100,
                         showProcess=FALSE,
                         saveRData=FALSE) {

  # determine the dataType
  if(identical(dataType, "RNA-Seq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  } else{
    stop("dataType is not specified correctly.")
  }

  # extract parameters: para0 - mean, para1 - dispersion/sd
  para00 <- para0.list[[1]]
  para01 <- para0.list[[2]]
  para10 <- para1.list[[1]]
  para11 <- para1.list[[2]]
  if(length(para00)!=length(para10))
    stop("Error: Unequal length between parameter vectors\n
         Please check the input.")

  numSim <- ST
  numEntry <- length(para00)
  entryList <- names(para00)
  if(is.null(entryList))
    stop("Error: unnamed mean vectors.")

  # initiate data object
  # count/abundance data
  # D1: each gene
  # D2: each sample
  # D3: each simulation
  dataMatrix <- array(data=0,
                    dim=c(numEntry, 2 * simNumRep, numSim),
                    dimnames=list(entryList,
                                    paste("Sample ID", 1:(2 * simNumRep)),
                                    paste("Sim", 1:numSim)))
  # Wald statistics for NB distribution / t Statistics for normal distribution
  # row: each gene / protein
  # column: each simulation
  testStatistics <- array(0, c(numEntry, numSim),
                           dimnames=list(entryList, paste("Sim", 1:numSim)))

  # packup into lists
  nullData <- list(dataMatrix=dataMatrix,
                   testStatistics=testStatistics,
                   parameter0=para00,
                   parameter1=para10,
                   T0=ST)

  alterData <- list(dataMatrix=dataMatrix,
                    testStatistics=testStatistics,
                    parameter0=para01,
                    parameter1=para11,
                    T1=ST)
  rm(dataMatrix, testStatistics)

  count <- 0
  totalIter <- numSim * numEntry
  for(i in 1:numEntry) {
    # Show Progress Details - show each distribution under simulation
    if(dataTypeSelect) {
      if(showProcess == TRUE)
        message(sprintf(">> Distributions: NB(%s, %s)\tand\t NB(%s, %s)",
                        round(para00[i],2),
                        round(para10[i],2),
                        round(para01[i],2),
                        round(para11[i],2)))
    }else{
        if(showProcess == TRUE)
          message(sprintf(">> Distributions: N(%s, %s)\tand\t N(%s, %s)",
                          round(para00[i],2),
                          round(para10[i],2),
                          round(para01[i],2),
                          round(para11[i],2)))
      }

    for(j in 1:numSim) {
      # simulate each distribution T0 times
      # null hypohesis: two groups should follow the same distribution
      simData0 <- if(dataTypeSelect) { simCounts(simNumRep,
                                                 para00[i],
                                                 para00[i],
                                                 para10[i],
                                                 para10[i] )
                     } else { simAbundance( simNumRep,
                                           para00[i],
                                           para00[i],
                                           para10[i],
                                           para10[i]) }
      simData1 <- if(dataTypeSelect) { simCounts(simNumRep,
                                                 para00[i],
                                                 para01[i],
                                                 para10[i],
                                                 para11[i] )
                     } else { simAbundance( simNumRep,
                                           para00[i],
                                           para01[i],
                                           para10[i],
                                           para11[i]) }
      # store simulated data in each simulation
      nullData$dataMatrix[i,,j] <- if(dataTypeSelect) simData0$counts
                                    else simData0$abundance
      alterData$dataMatrix[i,,j] <- if(dataTypeSelect) simData1$counts
                                     else simData1$abundance

      # try to calculate the test statistic,
      tmp0 <- if(dataTypeSelect) try(WaldTest(simData0))
               else try(tTestGLM(simData0))
      tmp1 <- if(dataTypeSelect) try(WaldTest(simData1))
               else try(tTestGLM(simData1))

      # if error occurs, test statistics remain as 0
      if (!inherits(tmp0,"try-error")) {
        nullData$testStatistics[i, j] <- tmp0
      }
      if (!inherits(tmp1,"try-error")) {
        alterData$testStatistics[i, j] <- tmp1
      }
      # progress bar
      if(showProcess == FALSE) {
        count <- count + 1
        progress(now=count, max=totalIter, word="of All Simulations Completed")
      }
    }
  }

  simulatedData <- list(nullData=nullData, alterData=alterData)
  attr(simulatedData, "dataType") <- dataType
  attr(simulateData, "repNum") <- simNumRep
  # Save data for later studies
  timestamp()
  message(">> Simulations Complete!")
  if(saveRData == TRUE) {
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    # simulated data name indicates the basic information
    filename.simulatedData <- RDataName(paste("simulated_Data","numRep",
                                              simNumRep, "numSim",
                                              numSim, sep="_"))
    save(simulatedData, file=paste0(getwd(), "/savedRData/",
                                    filename.simulatedData))
    message(">> Data saved in directory savedRData.")
    message(paste0(">> Size: ", round(file.size(
    paste0(getwd(), "/savedRData/", filename.simulatedData))/2^20, 2), " MB"))
  }
  return(simulatedData)
}
