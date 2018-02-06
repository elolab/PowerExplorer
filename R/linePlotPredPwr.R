#' Observe Predicted Power Within Different LFC Scales
#' @description With a complete power list and LFC list returned from
#' \code{predictPower}, power estimates can be observed dynamically
#' within specified LFC ranges.
#' @param PEObject the result object \code{PEObject} returned from
#' \code{predictPower}.
#' @param minLFC default as 0, the left edge of the LFC range within
#' which genes will be included in the graph.
#' @param maxLFC default as NA (the maximun value in data will be used) the
#' right edge of the LFC range within which genes will be included in the graph.
#' @param LFCscale the size of each unit when segmenting predicted power by LFC
#' @return plot(s) of power density under multiple sample sizes
#' @export
#' @import ggplot2
#' @import gridExtra
#' @examples
#' # load an example onject containing
#' # predicted power under different sample sizes
#' data(examplePredictedPower)
#' # plot a heatmap
#' plotPredPwr(examplePredictedPower,
#'                    LFCscale=1)
#' plotPredPwr(examplePredictedPower,
#'                    LFCscale=1)
#' #It is possible to observe power trend in different scales and ranges of LFCs
#' plotPredPwr(examplePredictedPower,
#'                    minLFC=0,
#'                    maxLFC=3,
#'                    LFCscale=0.5)
#' plotPredPwr(examplePredictedPower,
#'                    minLFC=0,
#'                    maxLFC=3,
#'                    LFCscale=1)
# Author: Xu Qiao
# Created: 25th, Sep, 2017
# Last Modifed: 10th, Jan, 2018
setGeneric(name="linePlotPredPwr",
           def=function(PEObject,
                        minLFC,
                        maxLFC,
                        LFCscale)
           {
             standardGeneric("linePlotPredPwr")
           }
)

setMethod("linePlotPredPwr", "PEObject", function(PEObject,
                                                  minLFC,
                                                  maxLFC,
                                                  LFCscale){
            plot.data <- extPlotData(PEObject, minLFC, maxLFC, LFCscale)
            SummaryTable <- aggregate(x = plot.data$power, 
                                      by = list(plot.data$comp, 
                                                plot.data$repNum), 
                                      FUN = function(x) mean(x, na.rm=TRUE))
            SummaryTable <- SummaryTable[order(SummaryTable$Group.1), ]
            repNums <- paste0("repNum:", unique(SummaryTable$Group.2))
            SummaryTable <- do.call(rbind, split(round(SummaryTable$x,2), 
                                                 SummaryTable$Group.1))
            sumTable <- tableGrob(SummaryTable,
                                  rows=row.names(SummaryTable),
                                  cols=repNums,
                                  theme=ttheme_default(base_size=8))
            minLFC <- attributes(plot.data)$info["minLFC"]
            maxLFC <- attributes(plot.data)$info["maxLFC"]
            LFCscale <- attributes(plot.data)$info["LFCscale"]  
            p.trend <-
                ggplot(plot.data, aes_string("repNum",
                                             "power",
                                             colour="lfc.range",
                                             group="lfc.range")) +
                geom_point(alpha=0.5, na.rm = TRUE) +
                geom_line(na.rm = TRUE) +
                facet_wrap(~comp) +
                ggtitle(label="Average Predicted Power within LFC ranges",
                        subtitle=
                          sprintf(
                            "segmented by every %s Log2FoldChange (minLFC: %s, maxLFC: %s)",
                            LFCscale, minLFC, maxLFC)) +
                xlab("Replicate Number") + ylab("Average power") +
                guides(colour=guide_legend(title="LFC range")) +
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
        })
