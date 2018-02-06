#' @import data.table
extPlotData <- function(PEObject, minLFC, maxLFC, LFCscale){
  predPwr <- PEObject@predPwr
  if(length(predPwr) == 0)
    stop("Predicted power results not found in PEObject!")
  if (missing(LFCscale)) stop("LFCscale is missing, with no default")
  lfc <- attributes(predPwr)$LFCRes
  molten.pp <- melt(predPwr)
  colnames(molten.pp) <- c("entry", "comp", "power", "repNum")
  molten.lfc <- melt(lfc)
  colnames(molten.lfc) <- c("entry", "comp", "lfc")
  molten.lfc$lfc <- abs(as.numeric(molten.lfc$lfc))
  molten.pp.lfc <- merge(molten.pp, molten.lfc, by=c("entry", "comp"))
  if (missing(minLFC)) minLFC <- round(min(molten.pp.lfc$lfc))
  if (missing(maxLFC)) maxLFC <- round(max(molten.pp.lfc$lfc))
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
                                 molten.pp.lfc$repNum), 
                         FUN=function(x) mean(x, na.rm=TRUE))
  colnames(plot.data) <- c("lfc.range", "comp", "repNum", "power")
  plot.data$repNum <- vapply(plot.data$repNum, 
                             function(x) as.numeric(substr(x, 9, 100)),
                             numeric(1))
  plot.data$lfc.range <- factor(plot.data$lfc.range,
                                levels = unique(plot.data$lfc.range))
  plot.data$repNum <- factor(plot.data$repNum,
                             levels = sort(unique(plot.data$repNum)))
  attr(plot.data, "info") <- c(minLFC=minLFC, 
                               maxLFC=maxLFC, 
                               LFCscale=LFCscale)
  return(plot.data)
}