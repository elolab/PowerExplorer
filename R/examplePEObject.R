#' @title An Example Predicted Power Object For Plot Tests
#'
#' @description This is an example PEObject resulted from
#' an example run on a publically available spiked protein 
#' dataset. The input dataset was subsetted, only concertration
#' 2fmol, 10fmol and 50fmol are kept, and all the 47 spike UPS1 
#' yeast proteins are kept, while only 200 background human proteins
#' are randomly picked to lower the overall dataset size.
#' Details: https://www.btk.fi/research/computational-biomedicine/
#' downloads/datasets/ups1-yeast-data/
#' 
#'
#' @docType data
#'
#' @usage data(examplePEObject)
#'
#' @keywords datasets
#'
#' @examples
#'data(examplePEObject)
#'
#'show(examplePEObject)
#'plotEstPwr(examplePEObject)
"examplePEObject"
