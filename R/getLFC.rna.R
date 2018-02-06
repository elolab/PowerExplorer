# retrieve LFC from generated DESeq object
# only can be called when dds has been calculated by DESeq2
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 19th, Dec, 2017

getLFC.rna <- function(dds) {
  coldata <- colData(dds) 
  groupnames <- unique(as.character(coldata[, "group"]))
  numGroup <- length(groupnames)
  comp_index <- combn(seq_len(numGroup), 2)
  colnames(comp_index) <- apply(comp_index, 2, function(x) 
    paste0(groupnames[x[1]], ".vs.",groupnames[x[2]]))
  res <- apply(comp_index, 2, function(x){
      groupOne <- groupnames[x[1]]
      groupTwo <- groupnames[x[2]]
      ext <- results(dds,contrast=c("group", groupTwo, groupOne))
      lfc <- ext$log2FoldChange
      names(lfc) <- rownames(ext)
      return(lfc)
    })
  return(res)
}

