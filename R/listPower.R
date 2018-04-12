#' List Estimated Power
#' @description show a top-list of power in numerical order,
#' or list the power selected genes/proteins.
#' @param inputObject the input inputObject.
#' @param decreasing logical; TRUE, decreasing order; FALSE, increasing order.
#' @param top the number of genes/proteins in the top list
#' @param selected default as NA; specify as a list of geneID or protein ID
#' to show power of a list of interested records.
#' @return a top list of power / power of a list of interested genes or proteins
#' @export
#'
#' @examples
#' data(exampleObject)
#' # show 10 top genes with high power (decreasing order)
#' listEstPwr(exampleObject, decreasing = TRUE, top = 10)
#' # show a list of interested genes
#' listEstPwr(exampleObject, 
#'             selected = c("ENSMUSG00000000303", 
#'                          "ENSMUSG00000087272", 
#'                          "ENSMUSG00000089921"))
#'
listEstPwr <- function(inputObject, decreasing = TRUE,  top = 20,  selected = NA){
  all.power <- data.frame(estPwr(inputObject), check.names = FALSE)
  result <- all.power[do.call(order, c(as.list(all.power),
                                       decreasing = decreasing)),,drop=FALSE]
  if(identical(selected, NA)) return(result[seq_len(top),, drop=FALSE])
  else{
    if(FALSE %in% (selected %in% rownames(inputObject)))
      stop(sprintf("ID %s not found in results.", 
                   paste(selected[!selected %in% rownames(inputObject)],
                         collapse = ", ")))
    return(all.power[selected,, drop=FALSE ]) 
  }
  return(result)
}


#' List Predicted Power
#' @description show a top-list of predicted power in numerical order,
#' or list the power selected genes/proteins.
#' @param inputObject the input inputObject.
#' @param decreasing logical; TRUE, decreasing order;
#' FALSE, increasing order.
#' @param top the number of genes/proteins in the top list
#' @param selected default as NA;
#' specify as a list of geneID or protein ID
#' to show power of a list of interested records.
#' @return a top list of power / power of a list of interested
#' genes or proteins for each sample size
#' @export
#'
#' @examples
#' data(exampleObject)
#' # show 10 top genes with high power (decreasing order)
#' listPredPwr(exampleObject, decreasing = TRUE, top = 10)
#'
#' # show a list of interested genes
#' listPredPwr(exampleObject, 
#'             selected = c("ENSMUSG00000000303", 
#'                          "ENSMUSG00000087272", 
#'                          "ENSMUSG00000089921"))
#'
listPredPwr <- function(inputObject, 
                        decreasing = TRUE, 
                        top = 20, 
                        selected = NA){
  result <- lapply(predPwr(inputObject), function(x) {
    temp <- data.frame(x, check.names = FALSE)
    temp <- temp[do.call(order, c(as.list(temp), 
                                  decreasing = TRUE)),,drop=FALSE]

    if(identical(selected, NA)) return(temp[seq_len(top),, drop=FALSE])
    else{
      if(FALSE %in% (selected %in% rownames(inputObject)))
        stop(sprintf("ID %s not found in results.", 
             paste(selected[!selected %in% rownames(inputObject)],
             collapse = ", ")))
      return(temp[selected,, drop=FALSE ]) 
    }
  })
  return(result)
}

