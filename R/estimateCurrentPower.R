#' @title Estimate Power of the Actual Data
#' @description Estimate power of comparison between each two groups based on
#' the data simulated from estimated normal distributions of entrys
#' in the entire dataset
#' @param inputDataMatrix a numeric raw Proteomics abundance data matrix,
#' in which rows correspond to entrys and columns to samples.
#' @param groupVec a vector indicating the grouping of samples
#' @param isLogTransformed logical; if \code{TRUE},
#' input data will be used directly without further transformation
#' @param dataType "RNA-Seq" or "Proteomics"
#' indictes the data type of the input data matrix.
#' @param minLFC the threshold for log2 fold change,
#' entrys with lower LFC are not included in the power calculation,
#' set to 0 if no threshold is needed.
#' @param alpha controlled false positive rate.
#' @param ST the number of simulations of abundance data generation and
#' repeated times of statistical test for each entry (>=100 recommended).
#' @param seed an integer seed for the random number generator.
#' @param showProcess logical; if \code{TRUE}, show the detailed information of
#' each simulation, used for debugging only.
#' @param saveSimulatedData logical; if \code{TRUE},
#' save the simulated data into RData with name pattern
#' "simulated_Data_numRep_X_numSim_XXX_XXXXX.RData"
#' @return a list of power estimates grouped in comparisons
#' between each two groups
#' @seealso \code{\link{predictSampleSizePower}}
#' predict power with incresing sample sizes
#' @export
#' @import stats
#' @import utils
#' @import methods
#' @examples
#' # Example 1: a random generated Proteomics dataset (10 DE, 100 non-DE)
#' # Note: Simulation times(ST) is specified as 10 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' data(exampleProteomicsData)
#' dataMatrix <- exampleProteomicsData$dataMatrix
#' groupVec <- exampleProteomicsData$groupVec
#'
#' # Run estimation without LFC filtration
#' currentPower <- estimateCurrentPower(dataMatrix, groupVec,
#'                                      dataType="Proteomics",
#'                                      isLogTransformed=FALSE,
#'                                      minLFC=0, alpha=0.05, ST=10, seed=123,
#'                                      showProcess=FALSE,
#'                                      saveSimulatedData=TRUE)
#'
#' # Run estimation with minLFC=1 (Fold Change=2)
#' currentPower2 <- estimateCurrentPower(dataMatrix, groupVec,
#'                                       dataType="Proteomics",
#'                                       isLogTransformed=FALSE,
#'                                       minLFC=1, alpha=0.05, ST=10, seed=123,
#'                                       showProcess=FALSE,
#'                                       saveSimulatedData=TRUE)
#'
#'# Example 2: a random generated RNA-Seq dataset (10 DE, 100 non-DE)
#' # Note: Simulation times(ST) is specified as 10 for shorter example runtime
#' #       For better performence, ST > 50 is recommended
#' data(exampleRNASeqData)
#' dataMatrix <- exampleRNASeqData$dataMatrix
#' groupVec <- exampleRNASeqData$groupVec
#'
#' # Run estimation without LFC filtration
#' currentPower <- estimateCurrentPower(dataMatrix, groupVec,
#'                                      dataType="RNA-Seq",
#'                                      isLogTransformed=FALSE,
#'                                      minLFC=0, alpha=0.05, ST=10, seed=123,
#'                                      showProcess=FALSE,
#'                                      saveSimulatedData=TRUE)
#'
#' # Run estimation with minLFC=1 (Fold Change=2)
#' currentPower2 <- estimateCurrentPower(dataMatrix, groupVec,
#'                                      dataType="RNA-Seq",
#'                                      isLogTransformed=FALSE,
#'                                      minLFC=1, alpha=0.05, ST=10, seed=123,
#'                                      showProcess=FALSE,
#'                                      saveSimulatedData=TRUE)
#'
#'
# Author: Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modifed: 9th, Jan, 2018

estimateCurrentPower <- function(inputDataMatrix, groupVec,
                                 isLogTransformed=FALSE,
                                 dataType=c("RNA-Seq","Proteomics"),
                                 minLFC=0.5, alpha=0.05, ST=100, seed=123,
                                 showProcess=FALSE,
                                 saveSimulatedData=TRUE) {
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
  if(length(groupVec) != ncol(inputDataMatrix))
    stop("groupVec and sample number vary in length.")
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup) #[!] unequal repNum

  # remove entries with too many zero counts(>=numRep)
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

  cat("Number of groups:\t", numGroup,
      "\nNumber of replicates:\t", numRep,
      "\nNumber of simulations:\t", ST,
      "\nMin. Log Fold Change:\t", minLFC,
      "\nFalse Postive Rate:\t", alpha,
      "\nTransformed:\t", isLogTransformed,
      "\n\n")

  currentPower <- list()
  for(g1 in 1:(numGroup-1)) {
    for(g2 in (g1+1):numGroup) {
      idx00 <- numRep * (g1 - 1) + 1
      idx01 <- g1*numRep
      idx10 <- numRep * (g2 - 1) + 1
      idx11 <- g2*numRep
      groups <- unique(groupVec)
      idx <- paste0(groups[g1], ".vs.", groups[g2]) # list name

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
      cat(sprintf("\nPower Estimation between group %s and group %s:\n",
                  groups[g1], groups[g2]))
      cat("\nQuantiles of absolute Log2 Fold Change:\n")
      show(quantile(round(abs(LFCList[[idx]]), 2),
                    probs=seq(0,1,0.1),
                    type=9))

      # LFC filter
      if(!is.na(minLFC) && minLFC > 0) {
        # filter the entrys with log fold change below minLFC
        abs.lfc <- abs(LFCList[[idx]])
        filter <- names(abs.lfc[abs.lfc >= minLFC])
        # apply filter
        mean.list <- lapply(mean.list, function(x) x[filter])
        if(dataTypeSelect) {
          theta.list <- lapply(theta.list, function(x) x[filter])
          cat(sprintf("\n%s of %s genes are over minLFC threshold %s:\n",
                      length(filter), length(abs.lfc), minLFC))
          }
        else{
          sd.list <- lapply(sd.list, function(x) x[filter])
          cat(sprintf("\n%s of %s proteins are over minLFC threshold %s:\n",
                      length(filter), length(abs.lfc), minLFC))
          }
      }
      # start simulation and power estimation
      simData <- simulateData(para0.list=mean.list,
                              para1.list=if(dataTypeSelect) get("theta.list")
                                          else get("sd.list"),
                              dataType=dataType,
                              simNumRep=numRep,
                              ST=ST,
                              showProcess=showProcess,
                              saveRData=saveSimulatedData)

      powerest <- calculatePower(simData, alpha=0.05,
                                 dataType=dataType,
                                 saveRData=FALSE,
                                 showOverallPower=!is.na(minLFC))
      currentPower[[idx]] <- powerest
      # free memory
      rm(simData, powerest)
    }
  }
  attr(currentPower, 'LFCList') <- LFCList
  attr(currentPower, 'minLFC') <- minLFC
  attr(currentPower, 'alpha') <- alpha
  attr(currentPower, 'ST') <- ST
  attr(currentPower, "dataType") <- dataType
  if(!("savedRData" %in% list.files())) dir.create("savedRData")
  filename.power <- RDataName(
    ifelse(dataTypeSelect,
           "[RNA-Seq] Estimated Power",
           "[Proteomics] Estimated Power"))
  savedRDataDir <- paste0(getwd(), "/savedRData/")
  save(currentPower, file=paste0(savedRDataDir, filename.power))
  message(">> Estimated power saved in savedRData directory.")
  message(paste0(">> Size: ", round(file.size(
    paste0(savedRDataDir, filename.power))/2^10, 2), " KB"))
  plotEstimatedPower(currentPower)
  return(currentPower)
}
