#' @title Estimate Power of the Actual Data
#' @description Estimate power of comparison between each two groups based on
#' the data simulated from estimated normal distributions of entrys
#' in the entire dataset
#' @param inputObject a numeric raw data matrix or SummarizedExperiment object
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; set to \code{TRUE},
#' if the input data is log transformed.
#' @param dataType "RNASeq" or "Proteomics"
#' indictes the data type of the input data matrix.
#' @param minLFC the threshold for log2 fold change,
#' entrys with lower LFC are not included in the power calculation,
#' set to 0 if no threshold is needed.
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each entry (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param enableROTS logical; if \code{TRUE}, Reproducibility-Optimized
#' Test Statistic (ROTS) will be used as the statistical model.
#' used as the statistical model.
#' @param paraROTS a \code{list} object containing addtional parameters 
#' passed to ROTS (if enabled), see \link{ROTS}.
#' @param showProcess logical; if \code{TRUE}, show the detailed
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
#' # Note: Simulation times(ST) is specified as 10 for shorter example runtime,
#' # ST > 50 is recommended
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
#'
#'
# Author: Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modifed: 13rd, March, 2018

estimatePower <- function(inputObject, groupVec,
                          isLogTransformed=FALSE,
                          dataType=c("RNASeq","Proteomics"),
                          minLFC=0.5, alpha=0.05, ST=100, seed=123,
                          enableROTS=FALSE,
                          paraROTS=list(B=1000,
                                        K=NULL,
                                        paired=FALSE,
                                        a1=NULL,
                                        a2=NULL,
                                        progress=FALSE),
                          showProcess=FALSE,
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
  ### numRep <- floor(length(groupVec)/numGroup) #[!] unequal repNum
  numRep <- as.numeric(table(groupVec))
  numGroup <- length(numRep)
  timestamp()
  cat("Sample groups:\t\t\t", paste(unique(groupVec), collapse=", "),
      "\nNum. of replicates:\t\t", paste(numRep, collapse = ", "),
      "\nNum. of simulations:\t\t", ST,
      "\nMin. Log Fold Change:\t\t", minLFC,
      "\nFalse Postive Rate:\t\t", alpha,
      "\nLog-transformed:\t\t", isLogTransformed,
      "\nROTS enabled:\t\t\t", enableROTS,
      "\n\n")

  # remove entries with too many zero counts
  numEntries.old <- nrow(dataMatrix)
  is.na(dataMatrix) <- !dataMatrix
  # keep rows with at least two reads in each group
  exZero <-
    apply(dataMatrix, 1, function(x){
      t_groupV <- table(factor(groupVec))
      c_zero <- table(factor(groupVec)[is.na(x)])
      c_read <- t_groupV - c_zero
      return(sum(c_read <2)!=0)
    })

  dataMatrix <- dataMatrix[!exZero, ]
  cat(sprintf("%s of %s entries are filtered due to excessive zero counts\n",
              numEntries.old - nrow(dataMatrix), numEntries.old))

  # Simulation procedure
  paraMatrices <- extParaMatrix(dataMatrix = dataMatrix,
                                groupVec = groupVec,
                                enableROTS=enableROTS,
                                isLogTransformed = isLogTransformed,
                                dataType = dataType,
                                minLFC = minLFC,
                                paraROTS=paraROTS,
                                seed = seed)
  
  estimatedPower <- lapply(paraMatrices, function(eachParaMatrix) {
    cat("\nSimulation in process, it may take a few minutes...\n")
    comp_idx <- attributes(eachParaMatrix)$Comparison
    repNumVec <- attributes(eachParaMatrix)$numRep
    # start simulation and power estimation
    cat(sprintf("\nPower Estimation between groups %s:\n",
                comp_idx))
    set.seed(seed)
    simData <- simulateData(eachParaMatrix,
                            dataType=dataType,
                            isLogTransformed=isLogTransformed,
                            enableROTS=enableROTS,
                            simNumRep=repNumVec,
                            minLFC=minLFC,
                            ST=ST,
                            showProcess=showProcess,
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
                  LFCRes=S4Vectors::DataFrame(attributes(paraMatrices)$LFCRes,
                                              check.names = FALSE),
                  minLFC=minLFC,
                  alpha=alpha,
                  ST=ST,
                  dataType=dataType,
                  simRepNumber=numRep,
                  estPwr=S4Vectors::DataFrame(estimatedPower,
                                              check.names = FALSE))

  if(saveResultData){
    if(!("savedRData" %in% list.files())) dir.create("savedRData")
    filename.result <-
      sprintf(
      ifelse(dataTypeSelect,
             "RNASeq_Results_%s.RData",
             "Proteomics_Results_%s.RData"),
      format(Sys.time(), "%H%M%S")
      )
    savedRDataDir <- paste0(getwd(), "/savedRData/")
    save(resObject, file=paste0(savedRDataDir, filename.result))
    message(">> Results saved in savedRData directory.")
    message(paste0(">> Size: ", round(file.size(
      paste0(savedRDataDir, filename.result))/2^10, 2), " KB"))
  }
  plotEstPwr(resObject)
  return(resObject)
}
