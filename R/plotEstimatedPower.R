#' Plot A Summary of Estimated Power
#' @description Produce a plot to summary the power estimated by function
#' \code{estimateCurrentPower}, the plot function is called in
#' \code{estimateCurrentPower}, but using it manually is possible
#' @param powerList the \code{power.list}
#' returned from \code{estimateCurrentPower}.
#' @return plot(s) of the summarised information on the estimated power
#' @import data.table
#' @import ggplot2
#' @import gridExtra
#'
#' @export
#' @examples
#' data(exampleEstimatedPower)
#' plotEstimatedPower(exampleEstimatedPower)
#'
# Author: Xu Qiao
# Created: 30th, Oct, 2017
# Last Modifed: 2nd, Jan, 2018
plotEstimatedPower <- function(powerList) {
  #determine dataType
  dataType <- attributes(powerList)$dataType
  if(identical(dataType, "RNA-Seq")) {
    dataTypeSelect <- TRUE
  } else if(identical(dataType, "Proteomics")) {
    dataTypeSelect <- FALSE
  }

  boxplot.data <- melt(powerList)
  colnames(boxplot.data) <- c("power", "comparison")
  boxplot.data$comparison <-  factor(boxplot.data$comparison,
                                     levels= unique(boxplot.data$comparison))

  barplot.data <- melt(lapply(powerList, function(x) c(sum(x >= 0.8),
                                                       sum(x < 0.8 & x >= 0.4),
                                                       sum(x < 0.4))))
  barplot.data$power.level <- rep(c("high (0.8, 1)",
                                    "low (0.4, 0.8)",
                                    "mere (0, 0.4)"),
                                  length(unique(barplot.data$L1)))
  colnames(barplot.data) <- c("num", "comparison", "power.level")
  barplot.data$comparison <-  factor(barplot.data$comparison,
                                     levels= unique(barplot.data$comparison))

  SummaryTable <- data.frame(
    Comparison=names(powerList),
    Gnum=sapply(powerList, length),
    A.P=sapply(powerList, function(x) round(mean(x, na.rm=TRUE), 2)),
    H.P=sapply(powerList, function(x)
                          sprintf("%s (%s%%)", sum(x >= 0.8),
                                  round((sum(x>=0.8)*100/length(x)),2))),
    L.P=sapply(powerList, function(x)
                          sprintf("%s (%s%%)", sum(x < 0.8 & x >= 0.4),
                                  round((sum(x < 0.8 & x >= 0.4)*
                                           100/length(x)), 2))),
    M.P=sapply(powerList, function(x)
                          sprintf("%s (%s%%)", sum(x < 0.4),
                                  round((sum(x < 0.4)*100 / length(x)), 2)))
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
                      color="comparison")) +
    geom_boxplot() +
    ggtitle(label="Boxplot",
            subtitle=paste("minLFC threshold:",
                           attributes(powerList)$minLFC)) +
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
                                            attributes(powerList)$minLFC)) +
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
               layout_matrix=matrix(c(1,2,3,3),
                                    nrow=2, byrow=TRUE))

  # if(savePlot) {
  #   if(!("savedPlots" %in% list.files())) dir.create("savedPlots")
  #   nam <- ifelse(dataTypeSelect,
  #                 "[RNA-Seq] Current Power Summary",
  #                 "[Proteomics] Current Power Summary")
  #   ggsave(filename=pngName(nam),
  #          plot=arranged.plots,
  #          path=paste0(getwd(), "/savedPlots"),
  #          width=6, height=4)
  #}
}
