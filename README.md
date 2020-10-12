# mutSigMapper

R package mutSigMapper aims to resolve a critical shortcoming of existing software for mutational 
signature analysis, namely that of finding parsimonious and biologically plausible exposures. 
By implementing a shot-noise-based model to generate spectral ensembles, this package addresses 
this gap and provides a quantitative, non-parametric assessment of statistical significance for 
the association between mutational signatures and observed spectra.

To install: install.packages("devtools") # if you have not installed the "devtools" package 
devtools::install_github("juliancandia/mutSigMapper")
