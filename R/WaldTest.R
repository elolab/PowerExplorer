# Wald statistical test
# fit count data to Negative Binomial Generalized Linear Model
# and obtain the Wald statistics (z value)

# nbdata - simulated read counts from simCounts
# formula - design formula for NB GLM model
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 28th, Dec, 2017
#' @import MASS
WaldTest <- function(nbdata, formula=counts ~ group) {
  # fit a NB generalized linear model
  tmp <- suppressWarnings(glm.nb(data=nbdata, formula=formula))
  if (!inherits(tmp,"try-error")) {
    GLMfit <- tmp
    statistics <- summary(GLMfit)$coefficients[2,"z value"]
    return(statistics^2)
  }
  else return(0)
}
