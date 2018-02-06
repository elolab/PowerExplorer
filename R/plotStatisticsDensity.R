#' Plot t-Statistics Density Between Null and Alternative Hypothesis
#' @description With saved simulated data, the t statistics under both null and
#' alternative hypotheses of a single protein can be observed by giving the
#' protein ID name.
#' @param simulatedData the simulated data object
#' @param entryID the gene name (string) or the row number (numeric)
#' in the countData oject
#' @param alpha controlled false positive rate
#' @return A plot presents the t-statistics density between null and
#' alternative hypothesis
#' @export
#' @import ggplot2
#' @import methods
#' @examples
#' #if the simulated data were saved
#' data("exampleSimulatedData")
#' plotStatisticsDensity(simulatedData=exampleSimulatedData,
#'                      entryID="A3KMH1",
#'                      alpha=0.05)
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 16th, Jan, 2018

plotStatisticsDensity <- function(simulatedData,
                                  entryID,
                                  alpha=0.05) {
  attrs <- attributes(simulatedData)
  if(attrs$dataType == "RNASeq") {
    dataTypeSelect <- TRUE
  }
  else if(attrs$dataType == "Proteomics") {
    dataTypeSelect <- FALSE
  }
  else {
    stop("Invalid simulatedData, no dataType found.")
  }
  entryList <- dimnames(simulatedData$nullData$testStatistics)[[1]]
  if(!entryID %in% entryList) {
    stop(sprintf("Incorrect Protein/Gene ID...\nTry something like %s",
                 paste(entryList[seq_len(3)], collapse=", ")))
  }
  x<-y<-NULL
  nullStatistics <- sqrt(simulatedData$nullData$testStatistics[entryID, ])
  alterStatistics <- sqrt(simulatedData$alterData$testStatistics[entryID, ])
  dens.null.xy <- with(density(nullStatistics), data.frame(x,y))
  dens.alter.xy <- with(density(alterStatistics), data.frame(x,y))
  dens.null.xy$hypothesis <- "null"
  dens.alter.xy$hypothesis <- "alternative"
  cutoff.statistics <- quantile(x=nullStatistics,
                                probs=1-alpha,
                                na.rm = TRUE,
                                type = 9)
  pwr <- mean(alterStatistics > cutoff.statistics, na.rm = TRUE)
  dens.all <- rbind(dens.null.xy, dens.alter.xy)
  p <- ggplot(data=dens.all, aes_string("x", "y")) +
    ggtitle(label=paste(if(dataTypeSelect)
                         {"Wald Statistics Density Plot of Gene"}
                        else
                         {"t Statistics Density Plot of Protein"}, entryID),
            subtitle=sprintf("Replicate number: %s\nSimulation depth: %s",
                             attrs$repNum, simulatedData$nullData$T0)) +
    geom_line(aes_string(colour="hypothesis")) +
    scale_colour_manual(values=c("blue", "red")) +
    geom_ribbon(data=dens.null.xy[dens.null.xy$x > cutoff.statistics, ],
                aes_string(ymax="y"),
                ymin=0, fill="red", colour=NA, alpha=0.4) +
    geom_ribbon(data=dens.alter.xy[dens.null.xy$x > cutoff.statistics, ],
                aes_string(ymax="y"),
                ymin=0, fill="blue", colour=NA, alpha=0.4)+
    xlab(if(dataTypeSelect) {"Wald Statistics"}
          else {"t Statistics"}) +
    ylab("Density") +
    geom_vline(xintercept=cutoff.statistics,
               colour="red", linetype="dashed") +
    geom_label(x=cutoff.statistics, y=max(dens.all$y)/2,
               label=paste("alpha: ", alpha), colour="red") +
    geom_label(x=max(dens.all$x)/2, y=max(dens.all$y)/3,
               label=paste("power: ", pwr), colour="orange")
  # if(savePlot) {
  #   if(!("savedPlots" %in% list.files())) dir.create("savedPlots")
  #   ggsave(filename=paste0("Statistics Plot of ", entryID, ".png"),
  #          plot=p, width=7, height=5, path=paste0(getwd(), "/savedPlots"))
  #   cat("Image saved to /savedPlots.")
  # }
  show(p)
}

