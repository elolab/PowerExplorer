## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----showDataExample, message=FALSE, warning=FALSE-----------------------
library(PowerExplorer)
data("exampleRNASeqData")
head(exampleRNASeqData$dataMatrix[,1:5])

## ------------------------------------------------------------------------
show(exampleProteomicsData$groupVec)

## ------------------------------------------------------------------------
colnames(exampleProteomicsData$dataMatrix)

## ----estimation run, fig.keep='none', message=FALSE, warning=FALSE, results='hide'----
library(PowerExplorer)
data("exampleProteomicsData")
estimatedPower <- estimateCurrentPower(inputDataMatrix = exampleProteomicsData$dataMatrix, 
                                       groupVec = exampleProteomicsData$groupVec, 
                                       isLogTransformed = FALSE, 
                                       dataType = "Proteomics", 
                                       minLFC = 0.5, 
                                       alpha = 0.05, 
                                       ST = 1, 
                                       seed = 345, 
                                       showProcess = FALSE, 
                                       saveSimulatedData = FALSE)

## ----plot estimated power, echo = FALSE----------------------------------
plotCurrentPower(powerList = estimatedPower, savePlot = FALSE)

## ----prediction run, fig.keep='none', message=FALSE, warning=FALSE, results='hide'----
data("exampleProteomicsData")
predictedPower <- predictSampleSizePower(inputDataMatrix = exampleProteomicsData$dataMatrix,
                                         groupVec = exampleProteomicsData$groupVec,
                                         isLogTransformed = FALSE,
                                         dataType = "Proteomics",
                                         rangeSimNumRep = c(5, 10, 15, 20),
                                         alpha = 0.05,
                                         ST = 1,
                                         seed = 345,
                                         showProcess = FALSE,
                                         saveSimulatedData = FALSE)

## ----LinePlot, fig.height=3.5--------------------------------------------
plotPredictedPower(predictedPower = predictedPower,
                   plotType = "lineplot",
                   LFCscale = 0.5,
                   savePlot = FALSE)

## ----Heatmap, fig.height=3.5---------------------------------------------
plotPredictedPower(predictedPower = predictedPower,
                   plotType = "heatmap",
                   LFCscale = 0.5,
                   savePlot = FALSE)

