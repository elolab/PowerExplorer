RDataName <- function(prefix) {
  postfix <- paste0("_", format(Sys.time(), "%H%M%S"))
  nam <- paste0(prefix, postfix, ".RData")
  nam
}
