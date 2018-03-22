# T Test based on Generalized Linear Model Fit
# Author: Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modifed: 5th, Dec, 2017
#'@importFrom stats glm gaussian
tTestGLM <- function(simdata, formula=data ~ group) {
  GLMfit <- glm(data=simdata, formula=formula, family=gaussian(link="log"))
  statistics <- summary(GLMfit)$coefficients[2, "t value"]
  abs(statistics)
}
