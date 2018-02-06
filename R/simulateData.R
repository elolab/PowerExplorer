# Simulation
# INPUT
# para0.list - a list of two mean vectors of genes in two case conditions
# para1.list - a list of two theta vectors of genes in two case conditions
# simNumRep - the number of replicates of each group for simulation
# ST - simulation times
# saveResultData - whether to save the simulated data into *.RData file
# OUTPUT
# simulatedData - a list with two list objects: nullData and alterData both
# contains count data, Wald statistics, mean and dispersion of each gene
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 29th, Dec, 2017
#' @import utils
#' @import plyr
simulateData <- function(paraMatrix,
                         dataType=c("RNASeq", "Proteomics"),
                         simNumRep,
                         ST=100,
                         showSimProcess = FALSE,
                         saveResultData=FALSE) {

  # determine dataType
  dataTypeSelect <- switch (dataType, 
                            RNASeq=TRUE, 
                            Proteomics=FALSE, 
                            stop("Incorrect dataType."))
  
  numEntry <- nrow(paraMatrix)
  entryList <- row.names(paraMatrix)
  if(is.null(entryList))
    stop("Error: gene/protein names not found from input paramter matrix.")
  
  # nodes <- detectCores()
  # cl <- makeCluster(nodes)
  # registerDoParallel(cl)
  simulatedData <- llply(seq_len(ST), 
                         .progress=progress_win(title = "Simulation in progress..."), 
                         function(x){
      tempMatrix <- apply(paraMatrix, 1, function(x) {
        # extract parameters: para0 - mean, para1 - dispersion/sd
        para0_0 <- x[1]
        para0_1 <- x[3] 
        para1_0 <- x[2]
        para1_1 <- x[4]
        # Show Progress Details - show each distribution under simulation
        
        distName <- ifelse(dataTypeSelect,"NB","N")
        if(showSimProcess == TRUE)
          message(sprintf(">> Distributions: %s(%s, %s)\tand\t %s(%s, %s)",
                          distName,
                          round(para0_0,2),
                          round(para1_0,2),
                          distName,
                          round(para0_1,2),
                          round(para1_1,2)))
        
      # simulate data for each gene/protein
      # null hypohesis: two groups should follow the same distribution
      simData0 <- 
        if(dataTypeSelect) { 
          simCounts(simNumRep, para0_0, para0_0, para1_0, para1_0)
        } else { 
          simAbundance(simNumRep, para0_0, para0_0, para1_0, para1_0) }
      simData1 <- 
        if(dataTypeSelect) { 
          simCounts(simNumRep, para0_0, para0_1, para1_0, para1_1)
        } else { 
          simAbundance(simNumRep, para0_0, para0_1, para1_0, para1_1) }
        
      statistics0 <- 
        if(dataTypeSelect) try(WaldTest(simData0))
      else try(tTestGLM(simData0))
      
      statistics1 <- 
        if(dataTypeSelect) try(WaldTest(simData1))
      else try(tTestGLM(simData1))
        
      # if error occurs, test statistics remain as 0
      if (inherits(statistics0,"try-error")) {
        statistics0 <- 0
      }
      if (inherits(statistics0,"try-error")) {
        statistics1 <- 0
      }
      # try to calculate the test statistic
      if(dataTypeSelect) {
        data0 <- simData0$counts
        data1 <- simData1$counts
      }
      else {
        data0 <- simData0$abundance
        data1 <- simData1$abundance
      }
      
      # collection all the simulated counts and statistics
      res_data <- c(data0, data1, statistics0, statistics1)
      return(res_data)
      })
      return(tempMatrix)
  })
  # stopCluster(cl)
  # record the attributes
  attr(simulatedData, "dataType") <- dataType
  attr(simulatedData, "repNum") <- simNumRep
  # Save data for later studies
  if(saveResultData) {
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    # simulated data name indicates the basic information
    filename.simulatedData <- RDataName(paste("simulated_Data","numRep",
                                              simNumRep, "numSim",
                                              ST, sep="_"))
    save(simulatedData, file=paste0(getwd(), "/savedRData/",
                                    filename.simulatedData))
    message(">> Data saved in directory savedRData.")
    message(paste0(">> Size: ", round(file.size(
    paste0(getwd(), "/savedRData/", filename.simulatedData))/2^20, 2), " MB"))
  }
  return(simulatedData)
}
