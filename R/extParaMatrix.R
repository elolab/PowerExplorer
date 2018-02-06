# extract distribution parameters from input data
extParaMatrix <- function(dataMatrix, groupVec,
                          isLogTransformed=FALSE,
                          dataType=c("RNASeq","Proteomics"),
                          minLFC=0.5,
                          seed=123){
  # determine dataType
  dataTypeSelect <- switch (dataType,
                            RNASeq=TRUE,
                            Proteomics=FALSE,
                            stop("Incorrect dataType."))
  timestamp()
  set.seed(seed)
  if(length(groupVec) != ncol(dataMatrix))
    stop("groupVec and sample number vary in length.")
  # determine input group and replicate numbers
  groupVec <- as.character(groupVec) # unfactorize the group vector
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup) #[!] unequal repNum

  # remove entries with too many zero counts
  numEntries.old <- nrow(dataMatrix)
  is.na(dataMatrix) <- !dataMatrix
  dataMatrix <- na.omit(dataMatrix)
  cat(sprintf("%s of %s entries are filtered due to excessive zero counts\n",
              numEntries.old - nrow(dataMatrix), numEntries.old))

  if(dataTypeSelect) {
    # construct a colData for DESeq2
    colData <- data.frame(group=groupVec, row.names=colnames(dataMatrix))
    res <- estimationByDESeq2(dataMatrix, colData)
    dispersion <- res$dispersion
    LFCRes <- res$LFCRes
    normalized.countData <- res$normCounts
  }
  else{

    # calcualte log fold changes between each two groups
    LFCRes <- calculateLFC.prot(dataMatrix,
                                 groupVec,
                                 isLogTransformed=isLogTransformed)
    normalized.countData <- dataMatrix
  }

  groups <- unique(groupVec)
  comp_index <- combn(seq_len(numGroup), 2)
  comp_index <- split(comp_index, col(comp_index))
  names(comp_index) <- lapply(comp_index, function(x) paste0(groups[x[1]], ".vs.",groups[x[2]]))
  paraMatrices <- lapply(comp_index, function(x){
    g1 <- x[1]; g2 <- x[2]
    idx00 <- numRep * (g1 - 1) + 1
    idx01 <- g1*numRep
    idx10 <- numRep * (g2 - 1) + 1
    idx11 <- g2*numRep
    comp_idx <- paste0(groups[g1], ".vs.",groups[g2])
    # extract data from each group
    data.case.g1 <- normalized.countData[, idx00:idx01]
    data.case.g2 <- normalized.countData[, idx10:idx11]

    if(dataTypeSelect) { #RNASeq
      # calculate the means of all genes
      mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE)
      mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE)
      paraMatrix <- cbind(mean.g1, dispersion,
                          mean.g2, dispersion)
    }
    else{ #Proteomics
      # log-transform the data if not isLogTransformed
      if(isLogTransformed == FALSE) {
        data.case.g1 <- log2(data.case.g1 +1)
        data.case.g2 <- log2(data.case.g2 +1)
      }

      meanSD.g1 <- apply(data.case.g1, 1, normDistMLE)
      meanSD.g2 <- apply(data.case.g2, 1, normDistMLE)
      paraMatrix <- t(rbind(meanSD.g1, meanSD.g2))
    }

    # LFC filter
    if(!is.na(minLFC) && minLFC > 0) {
      # filter the entrys with log fold change below minLFC
      abs.lfc <- abs(LFCRes[,comp_idx])
      filter <- names(abs.lfc[abs.lfc >= minLFC])
      # apply filter
      paraMatrix <- paraMatrix[filter, ]
      if(dataTypeSelect) {
        cat(sprintf("\n[%s] %s of %s genes are over minLFC threshold %s:\n",
                    comp_idx, length(filter), length(abs.lfc), minLFC))
      }
      else{
        cat(sprintf("\n[%s] %s of %s proteins are over minLFC threshold %s:\n",
                    comp_idx, length(filter), length(abs.lfc), minLFC))
      }
    }
    cat(sprintf("\n[%s] Log2 Fold Change Quantiles:\n", comp_idx))
    show(round(quantile(abs(LFCRes[,comp_idx]),
                  probs=seq(0,1,0.1),
                  type=9),2))
    attr(paraMatrix, "Comparison") <- paste0(groups[x[1]], ".vs.",groups[x[2]])
    return(paraMatrix)
})
  attr(paraMatrices, "LFCRes") <- LFCRes
  return(paraMatrices)
}
