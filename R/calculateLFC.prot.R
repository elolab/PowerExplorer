# Estimate Log2FoldChange between each two groups (Proteomics)
# OUTPUT
# a data frame containing LFC of each comparison
# Author: Xu Qiao (xu.qiao@utu.fi)
# Created: 22nd, Sep, 2017
# Last Modified: 5th, Jan, 2017
calculateLFC.prot <- function(abundanceData,
                              groupVec,
                              isLogTransformed=FALSE) {
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup)

  # remove proteins with too many zero counts(>=numRep)
  is.na(abundanceData) <- !abundanceData
  abundanceData <- na.omit(abundanceData)
  groupnames <- unique(groupVec)
  proteinNames <- rownames(abundanceData)
  output <- list()
  for(g1 in seq_len(numGroup-1)) {
    for(g2 in seq(g1+1, numGroup)) {
      idx00 <- numRep * (g1 - 1) + 1
      idx01 <- g1*numRep
      idx10 <- numRep * (g2 - 1) + 1
      idx11 <- g2*numRep

      # extract data from each group
      data.case.g1 <- abundanceData[, idx00:idx01]
      data.case.g2 <- abundanceData[, idx10:idx11]
      # log-transform the data if not isLogTransformed
      if(isLogTransformed == FALSE) {
        data.case.g1 <- log2(data.case.g1 +1)
        data.case.g2 <- log2(data.case.g2 +1)
      }

      mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE) # means of group 1
      mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE) # means of group 2
      # names
      groupOne <- groupnames[g1]
      groupTwo <- groupnames[g2]
      idx <- paste0(groupOne, ".vs.", groupTwo)
      lfc <- mean.g2 - mean.g1
      output[[idx]] <- lfc[proteinNames]
    }
  }
  output
}
