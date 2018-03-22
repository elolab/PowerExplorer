#' List Estimated Power
#' @description show a top-list of power in numerical order,
#' or list the power selected genes/proteins.
#' @param PEObject the input PEObject.
#' @param decreasing logical; TRUE, decreasing order; FALSE, increasing order.
#' @param top the number of genes/proteins in the top list
#' @param selected default as NA; specify as a list of geneID or protein ID
#' to show power of a list of interested records.
#' @return a top list of power / power of a list of interested genes or proteins
#' @export
#'
#' @examples
#' data(examplePEObject)
#' # show 10 top proteins with high power (decreasing order)
#' listEstPwr(examplePEObject, decreasing = TRUE, top = 10)
#' # show a list of interested proteins
#' listEstPwr(examplePEObject, selected = c("Protein_1", "Protein_11"))
#'
listEstPwr <- function(PEObject, decreasing = TRUE,  top = 20,  selected = NA){
  all.power <- data.frame(PEObject@estPwr, check.names = FALSE)
  result <- all.power[do.call(order, c(as.list(all.power),
                                       decreasing = decreasing)),,drop=FALSE]
  if(identical(selected, NA)) return(result[seq_len(top),, drop=FALSE])
  return(result[selected,, drop=FALSE])
}


#' List Predicted Power
#' @description show a top-list of predicted power in numerical order,
#' or list the power selected genes/proteins.
#' @param PEObject the input PEObject.
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
#' data(examplePEObject)
#' # show 10 top proteins with high power (decreasing order)
#' listPredPwr(examplePEObject, decreasing = TRUE, top = 10)
#'
#' # show a list of interested proteins
#' listPredPwr(examplePEObject, selected = c("Protein_1", "Protein_11"))
#'
listPredPwr <- function(PEObject, 
                        decreasing = TRUE, 
                        top = 20, 
                        selected = NA){
  result <- lapply(PEObject@predPwr, function(x) {
    temp <- data.frame(x, check.names = FALSE)
    temp <- temp[do.call(order, c(as.list(temp), 
                                  decreasing = TRUE)),,drop=FALSE]

    if(identical(selected, NA)) return(temp[seq_len(top),, drop=FALSE])
    return(temp[selected,, drop=FALSE ])
  })
  return(result)
}

