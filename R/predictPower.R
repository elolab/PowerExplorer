#' Estimate Power Under Increasing Sample Sizes
#' @description Simlilar to \code{\link{estimatePower}},
#' power estimations are performed under multiple increasing sample sizes
#' @param inputObject a numeric raw Proteomics abundance data matrix,
#' in which rows correspond to proteins and columns correspond to samples.
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; if \code{TRUE}, input data will be used
#' directly without further transformation
#' @param dataType "RNASeq" or "Proteomics" indictes the data type of
#' the input data matrix.
#' @param rangeSimNumRep a vector of sample sizes under which power
#' will be estimated
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each protein (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param showProcess logical; if \code{TRUE}, show the detailed information of
#' each simulation, used for debugging only.
#' @param saveResultData logical; if \code{TRUE}, save the simulated data
#' into RData with name pattern
#' "simulated_Data_numRep_X_numSim_XXX_XXXXX.RData"
#' @return a list of power predictions for each sample size, grouped in
#' comparisons between each two groups
#' @seealso \code{\link{estimatePower}} estimate power based on
#' actual data
#' @import DESeq2
#' @import utils
#' @import SummarizedExperiment
#' @export
#' @examples
#' # Example 1: a random generated Proteomics dataset (10 DE, 100 non-DE)
#' data(exampleProteomicsData)
#' dataMatrix <- exampleProteomicsData$dataMatrix
#' groupVec <- exampleProteomicsData$groupVec
#'
#' # Run estimation
#' # Note: Simulation times(ST) is specified as 5 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' predictedPower <- predictPower(dataMatrix, groupVec,
#'                                          isLogTransformed=FALSE,
#'                                          dataType="Proteomics",
#'                                          minLFC=0,
#'                                          rangeSimNumRep=c(5, 10, 15),
#'                                          alpha=0.05, ST=5, seed=123,
#'                                          showSimProcess=FALSE,
#'                                          saveResultData=FALSE)

#' # Example 2: a random generated RNASeq dataset (10 DE, 100 non-DE)
#' data(exampleRNASeqData)
#' dataMatrix <- exampleRNASeqData$dataMatrix
#' groupVec <- exampleRNASeqData$groupVec
#'
#' # Run estimation
#' # Note: Simulation times(ST) is specified as 5 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' predictedPower <- predictPower(dataMatrix, groupVec,
#'                                          isLogTransformed=FALSE,
#'                                          dataType="RNASeq",
#'                                          minLFC=0,
#'                                          rangeSimNumRep=c(5, 10, 15),
#'                                          alpha=0.05, ST=5, seed=123,
#'                                          showSimProcess=FALSE,
#'                                          saveResultData=FALSE)
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 10th, Jan, 2018

predictPower <- function(inputObject, groupVec,
                                   isLogTransformed=FALSE,
                                   dataType=c("RNASeq", "Proteomics"),
                                   minLFC = 0.5,
                                   rangeSimNumRep=NA,
                                   alpha=0.05, ST=100, seed=123,
                                   showSimProcess=FALSE, 
                                   saveResultData=FALSE) {
  dataMatrix <- switch (class(inputObject),
                        RangedSummarizedExperiment = assays(inputObject)$counts,
                        SummarizedExperiment = assays(inputObject)$counts,
                        ExpressionSet = exprs(inputObject),
                        inputObject
  )
  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numGroup <- length(unique(groupVec))
  timestamp()
  cat("Num. of groups:\t", numGroup,
      "\nNum. of replicates:\t", paste(rangeSimNumRep, collapse=", "),
      "\nNum. of simulations:\t", ST,
      "\nMin. Log Fold Change:\t", minLFC,
      "\nFalse Postive Rate:\t", alpha,
      "\nTransformed:\t", isLogTransformed,
      "\n\n")
  # Simulation procedure
  paraMatrices <- extParaMatrix(dataMatrix = dataMatrix, 
                                groupVec = groupVec, 
                                isLogTransformed = isLogTransformed, 
                                dataType = dataType, 
                                minLFC = minLFC, 
                                seed = seed)
  
  predictedPower <- lapply(rangeSimNumRep, function(eachRepNum) {
    cat(sprintf('\n##--Simulation with %s replicates per group--##\n', eachRepNum))
    
    eachRepPower <- lapply(paraMatrices, function(eachParaMatrix) {
      cat(sprintf("\n[repNum:%s] Simulation in process, it may take a few minutes...\n", 
                  eachRepNum))
      comp_idx <- attributes(eachParaMatrix)$Comparison
      # start simulation and power estimation
      cat(sprintf("\n[repNum:%s] Power Estimation between groups %s:\n",
                  eachRepNum, comp_idx))
      simData <- simulateData(eachParaMatrix,
                              dataType=dataType,
                              simNumRep=eachRepNum,
                              ST=ST,
                              showSimProcess=showSimProcess,
                              saveResultData=saveResultData)
      
      powerest <- calPwr(simData, alpha=0.05,
                         dataType=dataType,
                         saveResultData=FALSE,
                         showOverallPower=TRUE)
      temp <- apply(dataMatrix, 1, function(x) NA)
      temp[names(powerest)] <- powerest
      return(temp)
    })
    eachRepPower <- do.call(cbind, eachRepPower)
    return(eachRepPower)
  })
  names(predictedPower) <- paste0("repNum: ", rangeSimNumRep)
  
  switch (class(inputObject),
          SummarizedExperiment = {
            SE <- inputObject},
          {
            SE <- SummarizedExperiment(
              assays=list(counts=dataMatrix),
              rowData=DataFrame(row.names=rownames(dataMatrix)))
            colData(SE) <- 
              DataFrame(sampleName=factor(colnames(dataMatrix)),
                        group=factor(groupVec),
                        row.names = colnames(dataMatrix))
          }
  )
  
  resObject <- new("PEObject", SE, groupVec=groupVec,
                   LFCRes=DataFrame(attributes(paraMatrices)$LFCRes),
                   minLFC=minLFC,
                   alpha=alpha,
                   ST=ST,
                   dataType=dataType,
                   resultType="predictedPower",
                   simRepNumber=rangeSimNumRep,
                   predPwr=predictedPower)
  if(saveResultData){
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.power <- RDataName(ifelse(dataTypeSelect,
                                       "[RNASeq] Predicted Power",
                                       "[Proteomics] Predicted Power"))
    savedRDataDir <- paste0(getwd(), "/savedRData/")
    save(resObject, file=paste0(savedRDataDir, filename.power))
    message(">> Data saved in savedRData directory.")
    message(paste0(">> Size: ", round(file.size(
      paste0(savedRDataDir, filename.power))/2^10, 2), " KB")) 
  }
  return(resObject)
}
