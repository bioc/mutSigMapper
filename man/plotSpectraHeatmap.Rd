\name{plotSpectraHeatmap}
\alias{plotSpectraHeatmap}
\title{generates heatmap plots of spectra vs mutational signatures}
\description{This function generates a heatmap plot of associations between spectra (columns) and mutational signatures (rows), which may include statistical significance annotations based on shot noise modeling. By default, all spectra and signatures are plotted, but the user can provide custom subsets. Options are provided to plot only samples and/or signatures that satisfy statistical significance. A variety of  association metrics can be chosen via \code{type}; numerical results are returned.}
\usage{
plotSpectraHeatmap(x, spectra.set=NULL, sig.set=NULL, 
type=c("weights","cosine","entropy","correlation"), 
cor.method = c("pearson","kendall","spearman"), signif=TRUE, pval=0.05, 
signif.sig.only=FALSE, signif.sample.only=FALSE, xlab=NULL, ylab=NULL, 
color=NULL, breaks=NULL, breaks.scale.adj=TRUE, 
breaks.scale.adj.max=c("q99","q95","max"), scale="none", Rowv=FALSE, Colv=TRUE, 
na.color=NULL, cexRow=NULL, srtRow=0, cexCol=NULL, srtCol=45, margins=c(5,5), 
key=TRUE, key.title=NA, dendrogram="none", trace="none", notecol="red", 
cor.notecol="white", notecex=1.5, main.title=NULL, col.title="black", 
line.title=1, cex.title=0.95, subtitle=NULL, col.subtitle="red", 
line.subtitle=-1, cex.subtitle=0.65, \dots)
}
\arguments{
  \item{x}{\code{mutSigMapper} object.}
  \item{spectra.set}{Vector of sample names. If \code{NULL}, all samples are included.}
  \item{sig.set}{Vector of signature names. If \code{NULL}, all signatures are included.}
  \item{type}{Type of signature vs spectra matrix to show: \code{"weights"} (signature contributions calculated for each sample), \code{"cosine"} (cosine similarity), \code{"entropy"} (Jensen-Shannon divergence) and \code{"correlation"}. Default is \code{"weights"}.}
  \item{cor.method}{Correlation method (if \code{type="correlation"}). Default is \code{"pearson"}.}
  \item{signif}{Logical to display statistical significance annotations. Default is \code{TRUE}.}
  \item{pval}{P-value threshold for statistical significance. Default is 0.05.}
  \item{signif.sig.only}{Logical to only display statistically significant signatures. Default is \code{FALSE}.}
  \item{signif.sample.only}{Logical to only display statistically significant samples. Default is \code{FALSE}.}
  \item{xlab}{Label for x axis.}
  \item{ylab}{Label for y axis.}
  \item{color}{Heatmap color vector. Length must be one less than number of breaks.}
  \item{breaks}{Color breaks vector.}
  \item{breaks.scale.adj}{(If \code{breaks} is \code{NULL}) Logical to automatically adjust the breaks scale to the data being plotted. Default is \code{TRUE}.}
  \item{breaks.scale.adj.max}{(If \code{breaks} is \code{NULL} and \code{breaks.scale.adj} is \code{TRUE}) Breaks maximum: \code{"q99"} (99-th percentile), \code{"q95"} (95-th percentile), or \code{"max"} (maximum). Default is \code{"q99"}.}
  \item{scale}{Logical to scale the data for heatmap in either the row or column direction. Default is \code{"none"}.}
  \item{Rowv}{Logical to reorder rows (signatures) by hierarchical clustering. Default is \code{FALSE}.}
  \item{Colv}{Logical to reorder columns (samples) by hierarchical clustering. Default is \code{TRUE}.}
  \item{na.color}{Color to use for missing values.}
  \item{cexRow}{Size of row labels.}
  \item{srtRow}{Angle of row labels, in degrees from horizontal.}
  \item{cexCol}{Size of column labels.}
  \item{srtCol}{Angle of column labels, in degrees from horizontal.}
  \item{margins}{Numeric vector of length 2 containing the margins for column and row names, respectively.}
  \item{key}{Logical to display key. Default is \code{TRUE}.}
  \item{key.title}{Main title of the color key.}
  \item{dendrogram}{To draw dendrograms. Default is \code{"none"}.}
  \item{trace}{To display trace lines. Default is \code{"none"}.}
  \item{notecol}{Color of statistical significance annotations for \code{type="weights", "cosine"} or \code{"entropy"}.}
  \item{cor.notecol}{Color of statistical significance annotations for \code{type="correlation"}.}
  \item{notecex}{Size of significance annotations.}
  \item{main.title}{Custom title.}
  \item{col.title}{Title color.}
  \item{line.title}{Title line position.}
  \item{cex.title}{Title size.}
  \item{subtitle}{Custom subtitle.}
  \item{col.subtitle}{Subtitle color.}
  \item{line.subtitle}{Subtitle line position.}
  \item{cex.subtitle}{Subtitle size.}
  \item{\dots}{Additional plotting parameters to be passed to \code{gplots::heatmap.2()}}
}
\value{
A list with:
\item{sig_spectra}{Numerical data frame with the values plotted in the heatmap across signatures (rows) and samples (columns).}
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{mutSigMapper}}, \code{\link{plotSpectraCaterpillar}}
}
\examples{
\dontshow{
data(spectra)
map = mutSigMapper(spectra[,1:2],n_rdm=5)
weights = plotSpectraHeatmap(map)
}
\donttest{
data(spectra)
map = mutSigMapper(spectra)
weights = plotSpectraHeatmap(map)
cos_sim = plotSpectraHeatmap(map,type="cosine")
}
}
