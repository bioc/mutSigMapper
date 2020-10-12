\name{mutSigMapper-package}
\alias{mutSigMapper-package}
\docType{package}
\title{maps spectra to mutational signatures}
\description{
Provides a quantitative, non-parametric assessment of statistical significance for the association between mutational signatures and observed spectra.
}
\details{
\tabular{ll}{
Package: \tab mutSigMapper\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2020-10-05\cr
License: \tab GPL-3\cr
}
This package provides a full pipeline of analysis: 
\code{\link{mutSigMapper}} takes in \code{spectra} and \code{sig} data to assess the statistical significance of sample/signature associations via shot-noise modeling; 
\code{\link{plotSpectraCaterpillar}} generates caterpillar plots of signature weights for individual samples;
\code{\link{plotSpectraHeatmap}} generates heatmap plots of spectra vs mutational signatures based on weights, cosine similarity, Jensen-Shannon divergence, or correlation;
and \code{\link{export}} creates plain text output files for downstream
processing. Example spectra from published sources are provided. Curated built-in signature compendia extracted de novo from large-scale analysis of cancer genomes and exomes (COSMIC) as well as from cell lines treated with known mutagens (Kucab et al, Cell 2019) are also provided.  
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\references{
COSMIC (\emph{Catalogue Of Somatic Mutations In Cancer}), https://cancer.sanger.ac.uk/cosmic/signatures/index.tt

Kucab JE et al. 
\emph{A Compendium of Mutational Signatures of Environmental Agents}, Cell 177, 821-836.e16 (2019).

Candia J.
  \emph{mutSigMapper: an R package to map spectra to mutational signatures based on shot-noise modeling}, bioRxiv (2020).
}
\keyword{ package }
\seealso{
\code{\link{mutSigMapper}}, \code{\link{plotSpectraCaterpillar}}, \code{\link{plotSpectraHeatmap}}, \code{\link{export}}.
}
\examples{
\donttest{
data(spectra)
map = mutSigMapper(spectra)
plotSpectraCaterpillar(map,spectra.set="colon1")
cos_sim = plotSpectraHeatmap(map,type="cosine")
export(x=map,dest_dir=tempdir())
}
}
