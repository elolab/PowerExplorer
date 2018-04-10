#' Plot A Summary of Estimated Power
#' @description Produce a plot to summary the power estimated by function
#' \code{\link{estimatePower}}, the plot function is called in
#' \code{\link{estimatePower}}, but using it manually is possible
#' @param inputObject a result container object 
#' \code{\link{PowerExplorerStorage}} returned from 
#' \code{\link{estimatePower}}.
#' @return plot(s) of the summarised information on the estimated power
#' @importFrom  data.table melt
#' @import ggplot2
#' @import gridExtra
#' @export
#' @examples
#' data(exampleObject)
#' plotEstPwr(exampleObject)
#'
# Author: Xu Qiao
# Created: 30th, Oct, 2017
# Last Modifed: 18th, Feb, 2018

plotEstPwr <- function(inputObject) {
  estimatedPower <- estPwr(inputObject)
  if(is.null(nrow(estimatedPower)))
    stop("Estimated power results not found in input inputObject!")
  #determine dataType
  dataType <- parameters(inputObject)[["dataType"]]
  if(identical(dataType, "RNASeq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  }

  boxplot.data <- melt(as.matrix(estimatedPower))
  boxplot.data <- na.omit(boxplot.data)
  colnames(boxplot.data) <- c("entryID", "comparison", "power")
  
  numSum <- apply(estimatedPower, 2,
                  function(x) c(sum(x >= 0.8, na.rm = TRUE),
                                sum(x < 0.8 & x >= 0.4, na.rm = TRUE),
                                sum(x < 0.4, na.rm = TRUE)))
  rownames(numSum) <- c("high (0.8, 1)",
                        "low (0.4, 0.8)",
                        "mere (0, 0.4)")
  barplot.data <- melt(numSum)

  colnames(barplot.data) <- c("power.level", "comparison", "num")

  SummaryTable <- data.frame(
    Comparison=colnames(estimatedPower),
    Gnum=apply(numSum, 2, sum),
    A.P=apply(estimatedPower, 2, function(x)
                                    round(mean(x, na.rm=TRUE), 2)),
    H.P=apply(numSum, 2, function(x)
      sprintf("%s (%s%%)", x[1],
              round(x[1]*100/sum(x)),2)),
    L.P=apply(numSum, 2, function(x)
      sprintf("%s (%s%%)", x[2],
              round(x[2]*100/sum(x)),2)),
    M.P=apply(numSum, 2, function(x)
      sprintf("%s (%s%%)", x[3],
              round(x[3]*100/sum(x)),2))
  )
  names(SummaryTable) <- c("Comp.",
                           ifelse(dataTypeSelect,
                                  "Gene Num.",
                                  "Protein Num."),
                           "Avg. Power",
                           "H (0.8, 1)",
                           "L (0.4, 0.8)",
                           "M (0, 0.4)")

  power.boxplot <-
    ggplot(data=boxplot.data,
           aes_string("comparison",
                      "power",
                      group="comparison",
                      colour="comparison")) +
    geom_boxplot() +
    ggtitle(label="Boxplot",
            subtitle=paste("minLFC threshold:",
                           parameters(inputObject)[["minLFC"]])) +
    xlab("Comparison pair") + ylab("Power") +
    theme_light(base_size = 9)+
    theme(axis.text.x=element_text(angle=60, hjust=1),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.5, "cm"))

  power.barplot <- ggplot(data=barplot.data,
                          aes_string("comparison",
                                     "num",
                                     fill="power.level")) +
    geom_bar(stat="identity") +
    ggtitle(label="Barplot", subtitle=paste("minLFC threshold:",
                                            parameters(
                                              inputObject
                                              )[["minLFC"]])) +
    ylab(ifelse(dataTypeSelect,
                "Num of genes above minLFC",
                "Num of proteins above minLFC")) +
    xlab("Comparison pair") +
    theme_light(base_size = 9)+
    theme(axis.text.x=element_text(angle=60, hjust=1),
          legend.text=element_text(size=8),
          legend.key.size = unit(0.5, "cm"))

  sumTable <- tableGrob(SummaryTable,
                        rows=NULL,
                        theme=ttheme_default(base_size=8))

  grid.arrange(power.barplot,
               power.boxplot,
               sumTable,
               top="Estimated Power Summary",
               newpage=TRUE,
               layout_matrix=matrix(c(1,2,1,2,3,3),
                                    nrow=3, byrow=TRUE))
}
