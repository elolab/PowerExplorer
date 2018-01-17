#' @title An Example Predicted Power Object For Plot Tests
#'
#' @description This is an example predicted power object calculated from
#' a spiked-in proteomics dataset, the predictions are performed under
#' replicate number 10, 15, 20, 25, and 30.
#'
#' @docType data
#'
#' @usage data(examplePredictedPower)
#'
#' @keywords datasets
#'
#' @examples
#'data(examplePredictedPower)
#'plotPredictedPower(examplePredictedPower,
#'                   plotType="lineplot",
#'                   minLFC=0, maxLFC=4,
#'                   LFCscale=0.5, savePlot=FALSE)
#'
"examplePredictedPower"
