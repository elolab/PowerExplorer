# PowerExplorer
## Abstract
This is a power and sample size estimation tool for RNA-Seq and quantitative proteomics data.

`PowerExplorer` contains following main features:
- Estimation of power based on the current setting
- Prediction of power according to the future settings (e.g. increasing sample size)
- Visualizations of estimation and prediction results

## Introduction
Power and sample size estimation is still one of the important principles in designing next-generation sequencing experiments to discover differential expressions, a few methods on power estimation for RNA-Seq data have been studied, while the one specialized for proteomics data has not yet been developed. `PowerExplorer`is a power estimation and prediction tool currently applicable to RNA-Seq and quantitative proteomics experiments. The calculation procedure starts with estimating the distribution parameters of each gene or protein (following referred as entry for simplicity) accordingly, with the obtained prior distribution of each entry, a specified amount of simulations are executed to generate data (read counts for RNA-Seq and peptide abundance for proteomics) repetitively for each entry based on null and alternative hypotheses. Furthermore, the corresponding statistical tests (t-test or Wald-test) are performed and the test statistics are collected, eventually the result statistics will be summarized to calculate the statistical power.
