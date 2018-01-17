#' Estimate Power Under Increasing Sample Sizes
#' @description Simlilar to \code{\link{estimateCurrentPower}},
#' power estimations are performed under multiple increasing sample sizes
#' @param inputDataMatrix a numeric raw Proteomics abundance data matrix,
#' in which rows correspond to proteins and columns correspond to samples.
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; if \code{TRUE}, input data will be used
#' directly without further transformation
#' @param dataType "RNA-Seq" or "Proteomics" indictes the data type of
#' the input data matrix.
#' @param rangeSimNumRep a vector of sample sizes under which power
#' will be estimated
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each protein (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param showProcess logical; if \code{TRUE}, show the detailed information of
#' each simulation, used for debugging only.
#' @param saveSimulatedData logical; if \code{TRUE}, save the simulated data
#' into RData with name pattern
#' "simulated_Data_numRep_X_numSim_XXX_XXXXX.RData"
#' @return a list of power predictions for each sample size, grouped in
#' comparisons between each two groups
#' @seealso \code{\link{estimateCurrentPower}} estimate power based on
#' actual data
#' @import DESeq2
#' @import utils
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
#' predictedPower <- predictSampleSizePower(dataMatrix, groupVec,
#'                                          isLogTransformed=FALSE,
#'                                          dataType="Proteomics",
#'                                          rangeSimNumRep=c(5, 10, 15),
#'                                          alpha=0.05, ST=5, seed=123,
#'                                          showProcess=FALSE,
#'                                          saveSimulatedData=TRUE)

#' # Example 2: a random generated RNA-Seq dataset (10 DE, 100 non-DE)
#' data(exampleRNASeqData)
#' dataMatrix <- exampleRNASeqData$dataMatrix
#' groupVec <- exampleRNASeqData$groupVec
#'
#' # Run estimation
#' # Note: Simulation times(ST) is specified as 5 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' predictedPower <- predictSampleSizePower(dataMatrix, groupVec,
#'                                          isLogTransformed=FALSE,
#'                                          dataType="RNA-Seq",
#'                                          rangeSimNumRep=c(5, 10, 15),
#'                                          alpha=0.05, ST=5, seed=123,
#'                                          showProcess=FALSE,
#'                                          saveSimulatedData=TRUE)
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 10th, Jan, 2018

predictSampleSizePower <- function(inputDataMatrix, groupVec,
                                   isLogTransformed=FALSE,
                                   dataType=c("RNA-Seq", "Proteomics"),
                                   rangeSimNumRep=NA,
                                   alpha=0.05, ST=100, seed=123,
                                   showProcess=FALSE, saveSimulatedData=TRUE) {
  # determine the dataType
  if(identical(dataType, "RNA-Seq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  } else{
    stop("dataType is not specified correctly.")
  }

  timestamp()
  set.seed(seed)
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) #unfactorize the group vector
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup)
  # when no simulation replicate number is provided, abord function
  if(identical(rangeSimNumRep, NA))
    stop("the range of replicate numbers is not specified")

  # remove genes with too many zero counts(>=numRep)
  numEntries.old <- nrow(inputDataMatrix)
  is.na(inputDataMatrix) <- !inputDataMatrix
  inputDataMatrix <- na.omit(inputDataMatrix)
  cat(sprintf("%s of %s entries are filtered due to excessive zero counts\n",
              numEntries.old - nrow(inputDataMatrix), numEntries.old))


  if(dataTypeSelect) {
    # construct a colData for DESeq2
    colData <- data.frame(group=groupVec, row.names=colnames(inputDataMatrix))
    res <- estimationByDESeq2(inputDataMatrix, colData)
    dispersionVec <- res$dispersionVec
    LFCList <- res$LFCList
    normalized.countData <- res$normCounts
  }
  else{

    # calcualte log fold changes between each two groups
    LFCList <- calculateLFC.prot(inputDataMatrix,
                                 groupVec,
                                 isLogTransformed=isLogTransformed)
    normalized.countData <- inputDataMatrix
  }

  # group names
  groups <- unique(groupVec)
  compIndex <- apply(combn(groups, 2), 2, function(x)
                                          paste(x, collapse=".vs."))
  repIndex <- sapply(rangeSimNumRep, function(x) paste("repNum:", x))

  predictedPower <- array(data=NA,
                       dim=c(nrow(normalized.countData),
                             length(repIndex),
                             length(compIndex)),
                       dimnames=list(rownames(normalized.countData),
                                     repIndex,
                                     compIndex))

  # loop more than once when rangeSimNumRep is a vector of values
  cat("Number of groups:\t", numGroup,
      "\nNumber of replicates:\t", paste(rangeSimNumRep, collapse=", "),
      "\nNumber of simulations:\t", ST,
      "\nFalse Postive Rate:\t", alpha,
      "\nTransformed:\t", isLogTransformed,
      "\n\n")

  for(i in 1:length(rangeSimNumRep)) {
    simNumRep <- rangeSimNumRep[i]
    cat(sprintf('\n(%s / %s) Simulation with %s replicates per group:\n', i,
                length(rangeSimNumRep), simNumRep))
    for(g1 in 1:(numGroup-1)) {
      for(g2 in (g1+1):numGroup) {
        idx00 <- numRep * (g1 - 1) + 1
        idx01 <- g1*numRep
        idx10 <- numRep * (g2 - 1) + 1
        idx11 <- g2*numRep

        # extract data from each group
        data.case.g1 <- normalized.countData[, idx00:idx01]
        data.case.g2 <- normalized.countData[, idx10:idx11]

        if(dataTypeSelect) { #RNA-Seq
          # calculate the means of all genes
          mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE)
          mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE)
          mean.list <- list(mean.g1, mean.g2)
          theta.list <- list(dispersionVec, dispersionVec)
        }
        else{ #Proteomics
          # log-transform the data if not isLogTransformed
          if(isLogTransformed == FALSE) {
            data.case.g1 <- log2(data.case.g1 +1)
            data.case.g2 <- log2(data.case.g2 +1)
          }

          meanSD.g1 <- apply(data.case.g1, 1, normDistMLE)
          meanSD.g2 <- apply(data.case.g2, 1, normDistMLE)

          mean.list <- list(mean.g1=meanSD.g1[1,], mean.g2=meanSD.g2[1,])
          sd.list <- list(sd.g1=meanSD.g1[2,], sd.g2=meanSD.g2[2,])
        }
        idx <- paste0(groups[g1], ".vs.", groups[g2]) # list name
        cat(sprintf("\nPower Estimation between group %s and group %s:\n",
                    groups[g1], groups[g2]))

        # start simulation and power estimation
        simData <- simulateData(para0.list=mean.list,
                                para1.list=if(dataTypeSelect) get("theta.list")
                                          else get("sd.list"),
                                dataType=dataType,
                                simNumRep=simNumRep,
                                ST=ST,
                                showProcess=showProcess,
                                saveRData=saveSimulatedData)

        powerest <- calculatePower(simData, alpha=0.05, dataType=dataType,
                                  saveRData=FALSE, showOverallPower=FALSE)
        predictedPower[, repIndex[i], idx] <- powerest
        # free memory
        rm(simData, powerest)
        gc()
      }
    }
  }
  attr(predictedPower, 'LFCList') <- LFCList
  attr(predictedPower, 'simRepNumber') <- rangeSimNumRep
  attr(predictedPower, 'alpha') <- alpha
  attr(predictedPower, 'ST') <- ST
  attr(predictedPower, "dataType") <- dataType
  if(!("savedRData" %in% list.files())) dir.create("savedRData")
  filename.power <- RDataName(ifelse(dataTypeSelect,
                                     "[RNA-Seq] Predicted Power",
                                     "[Proteomics] Predicted Power"))
  savedRDataDir <- paste0(getwd(), "/savedRData/")
  save(predictedPower, file=paste0(savedRDataDir, filename.power))
  message(">> Data saved in savedRData directory.")
  message(paste0(">> Size: ", round(file.size(
    paste0(savedRDataDir, filename.power))/2^10, 2), " KB"))
  return(predictedPower)
}
