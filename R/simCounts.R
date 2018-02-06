# simulate read counts with given mean and dispersion (RNASeq)
# between two selected groups
# n: number of replicates of each group
# mu1: mean of group 1
# theta1: dispersion of group 1
# mu2: mean of group 2
# theta2: dispersion of group 2
# simCounts(5, 300,400,0.3,0.2) to see an example simulated reads
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 20th, Dec, 2017
simCounts <- function(n, mu1, mu2, theta1, theta2) {
  if(mu1 < 1) mu1 <- 1 # mean count lower than 1 is not allowed
  if(mu2 < 1) mu2 <- 1
  while(1) {
    counts.group1 <- rnbinom(n=n, mu=mu1, size=1/theta1)
    counts.group2 <- rnbinom(n=n, mu=mu2, size=1/theta2)
    # avoid zero/negative counts from both groups
    if(var(counts.group1) > 0 &
       var(counts.group2) > 0 &
       (sum(counts.group1 == 0) + sum(counts.group2 == 0)) < n) break()
  }
  cts <- c(counts.group1, counts.group2)
  group.vec <- rep(c(1,2), each=n)
  counts.btw.2groups <- data.frame(counts=cts, group=group.vec)
  return(counts.btw.2groups)
}
