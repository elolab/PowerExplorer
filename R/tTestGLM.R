# T Test based on Generalized Linear Model Fit
#'@import stats
# Author: Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modifed: 5th, Dec, 2017
tTestGLM <- function(simdata, formula=abundance ~ group) {
  GLMfit <- glm(data=simdata, formula=formula, family=gaussian(link="log"))
  statistics <- summary(GLMfit)$coefficients[2, "t value"]
  abs(statistics)
}
