#' @title Estimate Power of the Actual Data
#' @description Estimate power of comparison between each two groups based on
#' the data simulated from estimated normal distributions of entrys
#' in the entire dataset
#' @param inputObject a numeric raw data matrix or SummarizedExperiment object
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; if \code{TRUE},
#' input data will be used directly without further transformation
#' @param dataType "RNASeq" or "Proteomics"
#' indictes the data type of the input data matrix.
#' @param minLFC the threshold for log2 fold change,
#' entrys with lower LFC are not included in the power calculation,
#' set to 0 if no threshold is needed.
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each entry (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param showSimProcess logical; if \code{TRUE}, show the detailed 
#' information of each simulation, used for debugging only.
#' @param saveResultData logical; if \code{TRUE},
#' save the simulated data into RData with name pattern
#' "simulated_Data_numRep_X_numSim_XXX_XXXXX.RData"
#' @return a list of power estimates grouped in comparisons
#' between each two groups
#' @seealso \code{\link{predictPower}}
#' predict power with incresing sample sizes
#' @export
#' @import utils
#' @import methods
#' @import SummarizedExperiment
#' @importFrom Biobase exprs
#' @importFrom S4Vectors DataFrame
#' @examples
#' # Example 1: a random generated Proteomics dataset (10 DE, 100 non-DE)
#' # Note: Simulation times(ST) is specified as 10 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' data(exampleProteomicsData)
#' dataMatrix <- exampleProteomicsData$dataMatrix
#' groupVec <- exampleProteomicsData$groupVec
#'
#' # Run estimation without LFC filtration
#' resObject <- estimatePower(dataMatrix, groupVec,
#'                            dataType="Proteomics",
#'                            isLogTransformed=FALSE,
#'                            minLFC=0, alpha=0.05, 
#'                            ST=10, seed=123)
#'
#' # Run estimation with minLFC=1 (Fold Change=2)
#' resObject2 <- estimatePower(dataMatrix, groupVec,
#'                            dataType="Proteomics",
#'                            isLogTransformed=FALSE,
#'                            minLFC=1, alpha=0.05, 
#'                            ST=10, seed=123)
#'
#'# Example 2: a random generated RNASeq dataset (10 DE, 100 non-DE)
#' # Note: Simulation times(ST) is specified as 10 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' data(exampleRNASeqData)
#' dataMatrix <- exampleRNASeqData$dataMatrix
#' groupVec <- exampleRNASeqData$groupVec
#'
#' # Run estimation without LFC filtration
#' resObject <- estimatePower(dataMatrix, groupVec,
#'                              dataType="RNASeq",
#'                              isLogTransformed=FALSE,
#'                              minLFC=0, alpha=0.05, 
#'                              ST=10, seed=123)
#'
#' # Run estimation with minLFC=1 (Fold Change=2)
#' resObject2 <- estimatePower(dataMatrix, groupVec,
#'                             dataType="RNASeq",
#'                             isLogTransformed=FALSE,
#'                             minLFC=1, alpha=0.05, 
#'                             ST=10, seed=123)
#'
#'
# Author: Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modifed: 9th, Feb, 2018

estimatePower <- function(inputObject, groupVec,
                                 isLogTransformed=FALSE,
                                 dataType=c("RNASeq","Proteomics"),
                                 minLFC=0.5, alpha=0.05, ST=100, seed=123,
                                 showSimProcess=FALSE,
                                 saveResultData=FALSE) {
  if(missing(dataType) & length(dataType)!=1) 
    stop("Please tell the data type of the input.")
  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))

  dataMatrix <- switch (class(inputObject),
                        RangedSummarizedExperiment = assays(inputObject)$counts,
                        SummarizedExperiment = assays(inputObject)$counts,
                        ExpressionSet = Biobase::exprs(inputObject),
                        PEObject = {
                          if(dataType != inputObject@dataType)
                            stop("Incorrect data type.")
                          assays(inputObject)$counts},
                        as.matrix(inputObject)
  )
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup) #[!] unequal repNum
  timestamp()
  cat("Num. of groups:\t", numGroup,
      "\nNum. of replicates:\t", numRep,
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
  estimatedPower <- lapply(paraMatrices, function(eachParaMatrix) {
    cat("\nSimulation in process, it may take a few minutes...\n")
    comp_idx <- attributes(eachParaMatrix)$Comparison
    # start simulation and power estimation
    cat(sprintf("\nPower Estimation between groups %s:\n",
                comp_idx))
    simData <- simulateData(eachParaMatrix,
                            dataType=dataType,
                            simNumRep=numRep,
                            minLFC=minLFC,
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
  estimatedPower <- do.call(cbind, estimatedPower)
  
switch (class(inputObject),
    SummarizedExperiment = {
      SE <- inputObject},
    PEObject = {
      SE <- inputObject
    },
    {
      SE <- SummarizedExperiment(
        assays=list(counts=dataMatrix),
        rowData=S4Vectors::DataFrame(row.names=rownames(dataMatrix)))
      colData(SE) <- 
        S4Vectors::DataFrame(sampleName=factor(colnames(dataMatrix)),
                 group=factor(groupVec),
                 row.names = colnames(dataMatrix))
    }
  )
  
  resObject <- new("PEObject", SE, groupVec=groupVec,
                  LFCRes=S4Vectors::DataFrame(attributes(paraMatrices)$LFCRes),
                  minLFC=minLFC,
                  alpha=alpha,
                  ST=ST,
                  dataType=dataType,
                  simRepNumber=numRep,
                  estPwr=S4Vectors::DataFrame(estimatedPower))
  
  if(saveResultData){
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.power <- RDataName(
      ifelse(dataTypeSelect,
             "[RNASeq] Estimated Power",
             "[Proteomics] Estimated Power"))
    savedRDataDir <- paste0(getwd(), "/savedRData/")
    save(resObject, file=paste0(savedRDataDir, filename.power))
    message(">> Estimated power saved in savedRData directory.")
    message(paste0(">> Size: ", round(file.size(
      paste0(savedRDataDir, filename.power))/2^10, 2), " KB"))
  }
  plotEstPwr(resObject)
  return(resObject)
}
