#' PowerExplorer Object
#' @description An extended \code{SummarizedExpriment} object to 
#' contain input dataMatrix, grouping information, estimated power, 
#' predicted power, fold change estimates and other estimation
#' parameters.
#' @export
PEObject <- setClass("PEObject",
         contains="SummarizedExperiment",
         representation=representation(
           groupVec="character",
           estPwr="DataFrame",
           predPwr="list",
           LFCRes="DataFrame",
           simRepNumber="numeric",
           minLFC="numeric",
           alpha = "numeric",
           ST = "numeric",
           dataType = "character"
         )
)

#' show method for PEObject
#' @describeIn show method for PEObject objects
#' @param object a PEObject object as input
#' @return a summary of input PEObject object
#' @examples
#' data(examplePEObject)
#' show(examplePEObject)
setMethod("show", "PEObject", function(object) {
  cat("##--Parameters--##\n")
  cat("-dataType:", object@dataType, "\n")
  cat(sprintf("-original repNum: %s\n-Comparison groups: %s\n",
              length(object@groupVec)/length(unique(object@groupVec)),
              paste(unique(object@groupVec), collapse=", ")))
  cat(
    sprintf("-False positive rate: %s\n-LFC threshold: %s\n-Simulations: %s\n",
              object@alpha, object@minLFC, object@ST))
  cat("\n##--Log2 Fold Change Range--##\n")
  LFC <-
    apply(object@LFCRes, 2, function(x) round(range(abs(x), na.rm=TRUE), 2))
  rownames(LFC) <- c("minLFC", "maxLFC")
  print(LFC)
  cat("\n##--Average Estimated Power--##\n")
  if(nrow(object@estPwr) == 0){
    cat("NA\n")
  }else{
    avgEstPwr <- apply(object@estPwr, 2, function(x) mean(x, na.rm=TRUE))
    print(avgEstPwr)
  }

  cat("\n##--Average Predicted Power--##\n")
  if(length(object@predPwr) == 0){
    cat("NA\n")
  }else{
    avgPredPwr <- lapply(object@predPwr, function(x)
      apply(x, 2, function(y) mean(y, na.rm=TRUE)))
    print(avgPredPwr)
            }
          }
)
