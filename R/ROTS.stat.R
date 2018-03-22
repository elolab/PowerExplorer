# ROTS - reproducibility-optimized statistical test
# simdata: simulated data
# opt.para: optimization parameters calculated by ROTS

ROTS.stat<- function(simdata, opt.para = c(0, 1)){
  group <- simdata$group
  gi <- unique(group)
  x <- simdata$data[group==gi[1]]
  y <- simdata$data[group==gi[2]]
  n1 <- length(x)
  n2 <- length(y)
  pooledSD <- sqrt(((n1+n2)/n1*n2)*((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2))
  stat <- abs(mean(x) - mean(y)) / (opt.para[1] + opt.para[2]*pooledSD)
  return(stat)
}
