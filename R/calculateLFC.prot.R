# Estimate Log2FoldChange between each two groups (Proteomics)
# OUTPUT
# a data frame containing LFC of each comparison
# Author: Xu Qiao (xu.qiao@utu.fi)
# Created: 22nd, Sep, 2017
# Last Modified: 5th, Jan, 2017
#' @importFrom stats na.omit
calculateLFC.prot <- function(abundanceData,
                              groupVec) {
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup)

  # remove proteins with too many zero counts(>=numRep)
  is.na(abundanceData) <- !abundanceData
  abundanceData <- na.omit(abundanceData)
  groupnames <- unique(groupVec)
  proteinNames <- rownames(abundanceData)
  
  comp_index <- combn(seq_len(numGroup), 2)
  colnames(comp_index) <- apply(comp_index, 2, 
                                function(x) 
                                  paste0(groupnames[x[1]], 
                                         ".vs.",
                                         groupnames[x[2]]))
  res <- apply(comp_index, 2, function(x){
    g1 <- x[1]; g2 <- x[2]
    idx00 <- numRep * (g1 - 1) + 1
    idx01 <- g1*numRep
    idx10 <- numRep * (g2 - 1) + 1
    idx11 <- g2*numRep
    
    # extract data from each group
    data.case.g1 <- abundanceData[, idx00:idx01]
    data.case.g2 <- abundanceData[, idx10:idx11]
    mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE) # means of group 1
    mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE) # means of group 2
    # names
    # groupOne <- groupnames[g1]
    # groupTwo <- groupnames[g2]
    # idx <- paste0(groupOne, ".vs.", groupTwo)
    lfc <- mean.g2 - mean.g1
    return(lfc)
  })
  return(res)
}
