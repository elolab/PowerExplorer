#' @importFrom data.table melt
#' @importFrom stats aggregate
extPlotData <- function(inputObject, minLFC, maxLFC, LFCscale){
  predPwr <- predPwr(inputObject)
  if(length(predPwr) == 0)
    stop("Predicted power results not found in inputObject!")
  if (missing(LFCscale)) stop("LFCscale is missing, with no default")
  lfc <- as.matrix(LFCRes(inputObject))
  molten.pp <- melt(predPwr)
  colnames(molten.pp) <- c("entry", "comp", "power", "repNum")
  molten.lfc <- melt(lfc)
  colnames(molten.lfc) <- c("entry", "comp", "lfc")
  molten.lfc$lfc <- abs(as.numeric(molten.lfc$lfc))
  molten.pp.lfc <- merge(molten.pp, molten.lfc, by=c("entry", "comp"))
  if (missing(minLFC)) minLFC <- parameters(inputObject)[["minLFC"]]
  if (missing(maxLFC)) maxLFC <- max(molten.pp.lfc$lfc)
  if(minLFC >= maxLFC) stop("\nError: minLFC should be smaller than maxLFC.\n")
    
  if(minLFC < parameters(inputObject)[["minLFC"]]) {
    minLFC <- parameters(inputObject)[["minLFC"]]
    message(paste("\n[Correction] The actual minLFC in data is", minLFC))}
  
  molten.pp.lfc <- subset(molten.pp.lfc, lfc>=minLFC & lfc<=maxLFC)
  
  molten.pp.lfc$lfc.rg <- cut(molten.pp.lfc$lfc,
                              breaks=seq(minLFC, round(maxLFC+0.5), LFCscale),
                              right=FALSE)

  plot.data <- aggregate(molten.pp.lfc$power,
                         by=list(molten.pp.lfc$lfc.rg,
                                 molten.pp.lfc$comp,
                                 molten.pp.lfc$repNum), 
                         FUN=function(x) mean(x, na.rm=TRUE))
  SummaryTable <- aggregate(molten.pp.lfc$power, 
                            by=list(molten.pp.lfc$repNum, molten.pp.lfc$comp), 
                            FUN = function(x) mean(x, na.rm=TRUE))
  
  # SummaryTable2 <- vapply(levels(molten.pp.lfc$comp), 
  #                         function(x) unlist(table(
  #                           subset(molten.pp.lfc, comp == x)$lfc.rg)), 
  #                         FUN.VALUE = numeric(length(
  #                           unique(molten.pp.lfc$lfc.rg))))
  # SummaryTable2 <- cbind(levels(molten.pp.lfc$lfc.rg), SummaryTable2)
  # colnames(SummaryTable2)[1] <- "LFC range"
  # 
  colnames(SummaryTable) <- c("repNum", "Comp", "AP")
  
  colnames(plot.data) <- c("lfc.range", "comp", "repNum", "power")
  
  repStr <- SummaryTable$repNum
  compStr <- levels(SummaryTable$Comp)
  repNum <- vapply(strsplit(repStr, ": "), 
                   function(x) as.numeric(x[2]), numeric(1))
  SummaryTable <- SummaryTable[order(repNum), ]
  repStr <- repStr[order(repNum)]
  SummaryTable <- do.call(rbind, split(round(SummaryTable$AP,2),
                                       SummaryTable$Comp))
  SummaryTable <- cbind(compStr, SummaryTable)
  
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
  return(list(pd = plot.data, st1 = SummaryTable))
}