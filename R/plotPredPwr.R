#' Observe Predicted Power Within Different LFC Scales
#' @description With a complete power list and LFC list returned from
#' \code{\link{predictPower}}, power estimates can be observed dynamically
#' within specified LFC ranges.
#' @param inputObject a result container object 
#' \code{\link{PowerExplorerStorage}}
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
#' data(exampleObject)
#' plotPredPwr(exampleObject)
#' plotPredPwr(exampleObject)
#' #It is possible to observe power trend in different scales and ranges of LFCs
#' plotPredPwr(exampleObject, minLFC=0, maxLFC=2, LFCscale=0.5)
#' 
# Author: Xu Qiao
# Created: 25th, Sep, 2017
# Last Modifed: 21st, Feb, 2018

plotPredPwr <- function(inputObject, minLFC, maxLFC, LFCscale=1){
  res <- extPlotData(inputObject, minLFC, maxLFC, LFCscale)
  plot.data <- res[[1]]
  SumT <- res[[2]]
  # SumT2 <- res[[3]]
  repStr <- paste0("repNum:", levels(plot.data$repNum))
  sumTable1 <- tableGrob(SumT,
                        rows=NULL,
                        cols=c("Comp.", repStr),
                        theme=ttheme_default(base_size=8))
  
  # sumTable2 <- tableGrob(SumT2,
  #                       rows=NULL,
  #                       cols= colnames(SumT2),
  #                       theme=ttheme_default(base_size=8))
  
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
    #scale_color_brewer(palette="Spectral") +
    facet_wrap(~comp) +
    ggtitle(label="Average Predicted Power within LFC ranges",
            subtitle=
              sprintf(
    "segmented by every %s Log2FoldChange (minLFC: %s, maxLFC: %s)",
                LFCscale, minLFC, round(maxLFC, 2))) +
    xlab("Replicate Number") + ylab("Average power") +
    guides(colour=guide_legend(title="LFC range")) +
    theme_light(base_size = 9) +
    theme(axis.text.x=element_text(angle=60, hjust=1),
          legend.text=element_text(size=8),
          legend.position="right",
          legend.key.size = unit(0.5, "cm"))
  grid.arrange(p.trend,
               sumTable1,
               #sumTable2,
               newpage=TRUE,
               layout_matrix=matrix(c(1,1,2),
                                    nrow=3, byrow=TRUE))
  # return(list(SumT, SumT2))
}
