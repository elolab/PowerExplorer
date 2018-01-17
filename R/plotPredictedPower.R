#' Observe Predicted Power Within Different LFC Scales
#' @description With a complete power list and LFC list returned from
#' \code{predictSampleSizePower}, power estimates can be observed dynamically
#' within specified LFC ranges.
#' @param predictedPower the \code{ predictedPower} returned from
#' \code{predictSampleSizePower}.
#' @param plotType choose graph type between lineplot and heatmap to
#' observe the predicted power from different aspects
#' @param minLFC default as 0, the left edge of the LFC range within
#' which genes will be included in the graph.
#' @param maxLFC default as NA (the maximun value in data will be used) the
#' right edge of the LFC range within which genes will be included in the graph.
#' @param LFCscale the size of each unit when segmenting predicted power by LFC
#' @return plot(s) of power density under multiple sample sizes
#' @export
#' @import ggplot2
#' @import data.table
#' @import gridExtra
#' @examples
#' # load an example onject containing
#' # predicted power under different sample sizes
#' data(examplePredictedPower)
#' # plot a heatmap
#' plotPredictedPower(examplePredictedPower,
#'                    plotType="heatmap",
#'                    LFCscale=1)
#' plotPredictedPower(examplePredictedPower,
#'                    plotType="lineplot",
#'                    LFCscale=1)
#' #It is possible to observe power trend in different scales and ranges of LFCs
#' plotPredictedPower(examplePredictedPower,
#'                    plotType="lineplot",
#'                    minLFC=0,
#'                    maxLFC=3,
#'                    LFCscale=0.5)
#' plotPredictedPower(examplePredictedPower,
#'                    plotType="heatmap",
#'                    minLFC=0,
#'                    maxLFC=3,
#'                    LFCscale=1)
# Author: Xu Qiao
# Created: 25th, Sep, 2017
# Last Modifed: 10th, Jan, 2018
plotPredictedPower <- function(predictedPower,
                               plotType=c("lineplot", "heatmap"),
                               minLFC=0, maxLFC=NA,
                               LFCscale=NA) {
  if (missing(plotType)) stop("plotType is missing, with no default")
  if (missing(LFCscale)) stop("LFCscale is missing, with no default")
  lfc <- attributes(predictedPower)$LFCList
  lfc <- as.data.frame.list(lfc)
  lfc$entryName <- row.names(lfc)
  molten.pp <- melt(predictedPower)
  colnames(molten.pp) <- c("entry", "repNum", "comp", "power")
  molten.lfc <- melt(lfc, id.vars=c("entryName"))
  colnames(molten.lfc) <- c("entry", "comp", "lfc")
  molten.lfc$lfc <- abs(as.numeric(molten.lfc$lfc))
  molten.pp.lfc <- merge(molten.pp, molten.lfc, by=c("entry", "comp"))
  if(is.na(maxLFC)) maxLFC <- round(max(molten.pp.lfc$lfc))
  if(minLFC >= maxLFC) stop("Error: minLFC should be smaller than maxLFC.")
  if(maxLFC > round(max(molten.pp.lfc$lfc))) {
    maxLFC <- round(max(molten.pp.lfc$lfc))
    message(paste("[Correction] The actual maxLFC in data is", maxLFC))}

  molten.pp.lfc$lfc.rg <- cut(molten.pp.lfc$lfc,
                              breaks=seq(minLFC, maxLFC+LFCscale, LFCscale),
                              right=FALSE)
  plot.data <- aggregate(molten.pp.lfc$power,
                         by=list(molten.pp.lfc$lfc.rg,
                                 molten.pp.lfc$comp,
                                 molten.pp.lfc$repNum), FUN=mean)
  plot.data <- within(plot.data, {Group.3 <- substr(Group.3, 9, 100)})
  colnames(plot.data) <- c("lfc.range", "comp", "repNum", "power")
  plot.data$lfc.range <- factor(plot.data$lfc.range,
                                levels = unique(plot.data$lfc.range))
  plot.data$repNum <- factor(plot.data$repNum,
                             levels = unique(plot.data$repNum))
  SummaryTable <- aggregate(x = plot.data$power, 
                            by = list(plot.data$comp, plot.data$repNum), 
                            FUN = mean)
  SummaryTable <- SummaryTable[order(SummaryTable$Group.1), ]
  repNums <- paste0("repNum:", unique(SummaryTable$Group.2))
  SummaryTable <- do.call(rbind, split(round(SummaryTable$x,2), SummaryTable$Group.1))
  sumTable <- tableGrob(SummaryTable,
                        rows=row.names(SummaryTable),
                        cols=repNums,
                        theme=ttheme_default(base_size=8))
  numPlot <- dim(predictedPower)[3]
  if(identical(plotType,"lineplot")) {
    p.trend <-
      ggplot(plot.data, aes_string("lfc.range",
                                   "power",
                                   colour="repNum",
                                   group="repNum")) +
      geom_point(alpha=0.5) +
      geom_line() +
      facet_wrap(~comp) +
      ggtitle(label="Average Predicted Power within LFC ranges",
          subtitle=
            sprintf(
              "segmented by every %s Log2FoldChange (minLFC: %s, maxLFC: %s)",
                    LFCscale, minLFC, maxLFC)) +
      xlab("LFC segment") + ylab("Average power") +
      guides(colour=guide_legend(title="repNum")) +
      theme_light(base_size = 9) +
      theme(axis.text.x=element_text(angle=60, hjust=1),
            legend.text=element_text(size=8),
            legend.position="right",
            legend.key.size = unit(0.5, "cm"))
    
    grid.arrange(p.trend,
                 sumTable,
                 newpage=TRUE,
                 layout_matrix=matrix(c(1,2),
                                      nrow=2, byrow=TRUE))
    # if(savePlot) {
    #   if(!("savedPlots" %in% list.files())) dir.create("savedPlots")
    #   ggsave(filename=pngName("Prediected Power Summary(lineplot)"),
    #          plot=p.trend, path=paste0(getwd(), "/savedPlots"),
    #          width=7, height=5)
    #}
  }
  else if(identical(plotType, "heatmap")) {
    p.heat <- ggplot(plot.data, aes_string("repNum", "lfc.range")) +
      geom_tile(aes_string(fill="power")) +
      scale_fill_gradient2(low="blue", mid="yellow", high="red",
                           midpoint=0.6,
                           guide_colorbar(title="power", title.vjust = 1)) +
      xlab("Replicate number") + ylab("LFC range") +
      facet_wrap(~comp) +
      theme_light(base_size = 9)+
      theme(axis.text.x=element_text(angle=60, hjust=1),
            legend.text=element_text(size=8),
            legend.position="right",
            legend.key.size = unit(0.5, "cm"))+
      ggtitle(label="Average Predicted Power within LFC ranges",
              subtitle=
                sprintf(
              "segmented by every %s Log2FoldChange (minLFC: %s, maxLFC: %s)",
                        LFCscale, minLFC, maxLFC))
    
    grid.arrange(p.heat,
                 sumTable,
                 newpage=TRUE,
                 layout_matrix=matrix(c(1,2),
                                      nrow=2, byrow=TRUE))
    # if(savePlot) {
    #   if(!("savedPlots" %in% list.files())) dir.create("savedPlots")
    #   ggsave(filename=pngName("Prediected Power Summary(heatmap)"),
    #          plot=p.heat,
    #          path=paste0(getwd(), "/savedPlots"),
    #          width=7, height=5)
    #}
  }
  else {stop("Incorrect plotType.")}
}
