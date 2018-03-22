# Estimate Log2FoldChange between each two groups
# OUTPUT
# a data frame containing LFC of each comparison
# Author: Xu Qiao (xu.qiao@utu.fi)
# Created: 22nd, Sep, 2017
# Last Modified: 5th, Jan, 2017
#' @importFrom stats na.omit
calLFC <- function(dataMatrix, groupVec) {
  numGroup <- length(unique(groupVec))
  numRep <- floor(length(groupVec)/numGroup)
  groupnames <- unique(groupVec)
  proteinNames <- rownames(dataMatrix)
  idx <- split(seq_along(groupVec),
               rep(unique(groupVec), table(groupVec)))
  comp_index <- combn(seq_len(numGroup), 2)
  colnames(comp_index) <- apply(comp_index, 2,
                                function(x)
                                  paste0(groupnames[x[1]],
                                         ".vs.",
                                         groupnames[x[2]]))
  res <- apply(comp_index, 2, function(x){
    g1 <- x[1]; g2 <- x[2]

    # extract data from each group
    data.case.g1 <- dataMatrix[, idx[[g1]]]
    data.case.g2 <- dataMatrix[, idx[[g2]]]
    mean.g1 <- rowMeans(data.case.g1, na.rm = TRUE) # means of group 1
    mean.g2 <- rowMeans(data.case.g2, na.rm = TRUE) # means of group 2
    lfc <- mean.g2 - mean.g1
    return(lfc)
  })
  return(res)
}

