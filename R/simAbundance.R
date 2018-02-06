# Generate abandance data from normal distribution with specified mean and sd
# n: sample size
# mu1: mean of group 1
# mu2: mean of group 2
# sd1: standard deviation of group 1
# sd2: standard deviation of group 2
# OUTPUT:
# a data frame with abundance and group
# Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modification: 3rd, Jan, 2018
#' @importFrom  stats rnorm
simAbundance <- function(n, mu1, mu2, sd1, sd2) {
  if(mu1 < 1) mu1 <- 1 # mean count lower than 1 is not allowed
  if(mu2 < 1) mu2 <- 1
  while(1) { #only non-negative values are allowed
    abd.group1 <- rnorm(n=n, mean=mu1, sd=sd1)
    abd.group2 <- rnorm(n=n, mean=mu2, sd=sd2)
    neg.g1 <- sum(abd.group1 < 0)
    neg.g2 <- sum(abd.group2 < 0)
    if(sum(neg.g1, neg.g2) == 0) break()
    }
  abds <- c(abd.group1, abd.group2)
  group.vec <- rep(c(1,2), each=n)
  abds.btw.2groups <- data.frame(abundance=abds, group=group.vec)
  return(abds.btw.2groups)
}
