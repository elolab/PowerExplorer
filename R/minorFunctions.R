# Minor Functions
# progress bar
# Author: Xu Qiao
# Created: 21st, Sep, 2017
# Last Modifed: 29th, Dec, 2017

progress <- function(now, max, word="Completed") {
  percent <- now / max * 100
  cat(sprintf(paste0('\r>> [%-50s] %d%% ', word,'...'),
              paste(rep('=', percent / 2), collapse=''),
              floor(percent)))
  if (now == max)
    cat('\n')
}

pngName <- function(prefix) {
  postfix <- paste0("_", format(Sys.time(), "%H%M%S"), sample(0:9, 1))
  nam <- paste0(prefix, postfix, ".png")
  nam
}

RDataName <- function(prefix) {
  postfix <- paste0("_", format(Sys.time(), "%H%M%S"))
  nam <- paste0(prefix, postfix, ".RData")
  nam
}



# save dispersion plot from DESeq2 object
# saveGeneDispPlot <- function(input_dds, ymin=1e-04) {
#   png(filename=pngName("DispEstiPlot"), width=700, height=700, res=100)
#   plotDispEsts(input_dds, ymin=ymin)
#   dev.off()
# }
