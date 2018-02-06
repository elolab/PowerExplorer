setClass("PEObject",
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
           dataType = "character",
           resultType = "character"
         )
)



