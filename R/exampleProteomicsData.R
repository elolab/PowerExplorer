#' @title Randomly Generated Proteomics Dataset
#'
#' @description This is a randomly generated proteomics dataset with
#' 130 protein entries (rows) and 15 samples (columns) in
#' 3 sample groups A, B and C, the log2 fold change (LFC)
#' between group B and A is specified as 1,
#' between C and B is also 1, thus the LFC is 2 between C and A.
#'
#' @docType data
#'
#' @usage data(exampleProteomicsData)
#'
#' @format An list contains \code{"dataMatrix"} and \code{"groupVec"}
#'
#' @keywords datasets
#'
#' @examples
#' data(exampleProteomicsData)
#' head(exampleProteomicsData$dataMatrix)
#'
"exampleProteomicsData"
