\name{plotSpectraCaterpillar}
\alias{plotSpectraCaterpillar}
\title{generates caterpillar plot of signature weights for individual samples}
\description{
For each chosen sample, this function generates a caterpillar plot of signature weights. By default, signatures are ordered by p-value and those below the significance threshold are marked with an asterisk. Weights from the input sample are shown as red crosses, whereas weights from the 
generated random spectra are shown as blue circles (median values) and bars (confidence interval ranges). Options are provided to customize signature selection and display.}
\usage{
plotSpectraCaterpillar(x, spectra.set=NULL, sig.set=NULL, order_sig_by_pval=TRUE, 
signif=TRUE, pval=0.05, CI=0.95, xlab=NULL, ylab=NULL, cexRow=NULL, cex.lab=0.95, 
main.title=NULL, col.title="black", line.title=1.5, cex.title=0.85, subtitle=NULL, 
col.subtitle="darkgray", line.subtitle=0.5, cex.subtitle=0.55, legend=TRUE)
}
\arguments{
  \item{x}{\code{mutSigMapper} object.}
  \item{spectra.set}{Vector of samples names. If \code{NULL}, all samples are included.}
  \item{sig.set}{Vector of signature names. If \code{NULL}, all signatures are included.}
  \item{order_sig_by_pval}{Logical to order signatures by p-value. Default is \code{TRUE}.}
  \item{signif}{Logical to display asterisks in statistically significant signature labels. Default is \code{TRUE}.}
  \item{pval}{P-value threshold for statistical significance. Default is 0.05.}
  \item{CI}{Confidence interval (normalized to 1) around the median for random spectra (shown as blue bars). Default is 0.95.}
  \item{xlab}{Label for x axis.}
  \item{ylab}{Label for y axis.}
  \item{cexRow}{Size of row labels.}
  \item{cex.lab}{Axis label size.}
  \item{main.title}{Custom title.}
  \item{col.title}{Title color.}
  \item{line.title}{Title line position.}
  \item{cex.title}{Title size.}
  \item{subtitle}{Custom subtitle.}
  \item{col.subtitle}{Subtitle color.}
  \item{line.subtitle}{Subtitle line position.}
  \item{cex.subtitle}{Subtitle size.}
  \item{legend}{Logical to display legend. Default is \code{TRUE}.}
}
\value{
None.
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{mutSigMapper}}, \code{\link{plotSpectraHeatmap}}
}
\examples{
\dontshow{
data(spectra)
map = mutSigMapper(spectra[,1:2],n_rdm=5)
plotSpectraCaterpillar(map, spectra.set="colon1")
}
\donttest{
data(spectra)
map = mutSigMapper(spectra)
plotSpectraCaterpillar(map, spectra.set="colon1")
}
}
