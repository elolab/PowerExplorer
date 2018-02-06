
# Author: Xu Qiao 
# Created: 5th, Dec, 2017
# Last Modifed: 19th, Dec, 2017
# formula reference:
# https://www.statlect.com/fundamentals-of-statistics/
#' @importFrom stats optim
normDistMLE <- function(x) {
  logLikeFunNorm <- function(paraVec) {
    # Log of likelihood of a normal distribution
    # paravec[1] - mean
    # paravec[2] - standard deviation
    # x - set of observations.
    x <- na.omit(x)
    n <- length(x)
    -(n/2)*log(2*pi)-(n/2)*log(paraVec[2]^2)-
      (1/(2*paraVec[2]^2))*sum((x-paraVec[1])^2)
  }
  # maximum likelihood estimation to find the most likely parameters
  MLE <- optim(c(0.1,0.1), # initial values for mu and sigma
              fn=logLikeFunNorm, # function to maximize
              method="L-BFGS-B", # this method lets set lower bounds
              lower=0.00001, # lower limit for parameters
              control=list(fnscale=-1), # maximize the function
              hessian=TRUE # calculate Hessian matricce because
              # we will need for confidence intervals
  )
  # return the mean and sd
  MLE$par
}
