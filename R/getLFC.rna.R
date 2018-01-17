# retrieve LFC from generated DESeq object
# only can be called when dds has been calculated by DESeq2
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 19th, Dec, 2017

getLFC.rna <- function(dds, colData) {
  groupnames <- levels(colData$group)
  numGroup <- length(groupnames)
  numComp <- numGroup * (numGroup - 1) / 2

  output <- list()
  for(g1 in 1:(numGroup-1)) {
    for(g2 in (g1+1):numGroup) {
      groupOne <- groupnames[g1]
      groupTwo <- groupnames[g2]
      res <- results(dds,contrast=c("group", groupTwo, groupOne))
      # lfc <- res$log2FoldChange
      # names(lfc) <- rownames(res)
      lfc <- res$log2FoldChange
      names(lfc) <- rownames(res)
      idx <- paste0(groupOne, ".vs.", groupTwo)
      output[[idx]] <- lfc
    }
  }
  output
}
