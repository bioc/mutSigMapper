\name{mutSigMapper}
\alias{mutSigMapper}
\title{maps spectra to mutational signatures}
\description{
This method maps spectra to mutational signatures based on a robust, non-parametric statistical approach. For each spectrum (provided as the number of mutations observed across channels, e.g. single-substitution trinucleotide contexts), non-negative least square regression procedures are used to calculate the contributions of each mutational signature. Shot noise (implemented as a Poisson or negative binomial process) generates an ensemble of random spectra that enables the assessment of statistical significance for those contributions. The resulting object can be plotted and further analyzed using built-in functions provided in the package.
}
\usage{
mutSigMapper(spectra, sig=NULL, ref=c("cosmic_v2","cosmic_v3",
"cosmic_v3_exome","cosmic_v3.1","mutagen53"), 
method=c("MutationalPatterns","deconstructSigs"),
sig.bkg.adj=c("none","1/exome","1/genome","exome/genome",
"genome/exome","custom"), sig.bkg.adj.custom=NULL,
signature.cutoff=0.06, noise=c("poisson","neg.binom"), 
neg.binom.size=NULL, n_rdm=200,
save_obj=FALSE, dest_dir=getwd(), dest_dir_create=TRUE, 
dest_dir_create_recur=FALSE, dest_obj="mutSigMap.Robj")
}
\arguments{
\item{spectra}{Input numerical matrix or data frame with channels as rows and samples as columns. Data must be absolute counts (i.e. not normalized spectra). Channel and sample labels should be provided as row and column names, respectively.}
\item{sig}{Input numerical matrix or data frame with channels as rows and signatures as columns. Channel and signature labels should be provided as row and column names, respectively. Channel labels must agree with those of the input spectra.}
\item{ref}{If a \code{sig} data frame is not provided, a reference built-in compendium can be chosen. Options are \code{"cosmic_v2"} (Cosmic v2, March 2015), \code{"cosmic_v3"} and \code{"cosmic_v3_exome"} (the genome and exome releases from Cosmic v3, May 2019), \code{"cosmic_v3.1"} (Cosmic v3.1, June 2020), and \code{"mutagen53"} (mutagen compendium by Kucab et al, Cell (2019)).}
\item{method}{Non-negative least squares regression procedure. Options are \code{"MutationalPatterns"} (default) and \code{"deconstructSigs"}.}
\item{sig.bkg.adj}{Distribution of background counts to adjust the signatures to match the spectra data source. See \bold{Details} section below.}
\item{sig.bkg.adj.custom}{For \code{sig.bkg.adj="custom"}, data frame to adjust signatures. It should be provided in the same format as datasets \code{ref_bkg_exome} and \code{ref_bkg_exome}.} 
\item{signature.cutoff}{Signature contributions (weights) must be above this threshold.}
\item{noise}{Statistical noise model to generate random spectra. Options are \code{"poisson"} and \code{"neg.binom"}.}
\item{neg.binom.size}{For negative binomial noise, this is the size parameter. Must be a positive number. As the size is increased, the distribution's variance decreases. Size values much larger than the distribution's mean yield negative binomial distributions similar to Poisson.}
\item{n_rdm}{Number of random spectra generated for each sample. Default is 200.}
\item{save_obj}{Logical to save the \code{mutSigMapper} object. Default is \code{FALSE}.}
\item{dest_dir}{Destination directory. Default is the working directory.}
\item{dest_dir_create}{Creates destination directory if it does not exist already. Default is \code{TRUE}.}
\item{dest_dir_create_recur}{Creates destination directory recursively if it does not exist already. Default is \code{FALSE}.}
\item{dest_obj}{Name for output \code{mutSigMapper} object.}
}
\details{
Input spectra must be provided as unnormalized mutation counts. An input signature matrix or data frame may be provided, but the package offers a variety of built-in reference compendia. Two different built-in procedures for non-negative least squares regression are included. Signatures may be adjusted to match the spectra data source. Signatures based on human genome abundance (such as \code{"cosmic_v2"}, \code{"cosmic_v3"}, \code{"cosmic_v3.1"} and \code{"mutagen53"}) used to analyze spectra extracted from whole genome sequencing (WGS) do not require adjustment; for spectra extracted from whole exome sequencing (WES), however, one needs to adjust genome-based signatures using \code{sig.bkg.adj="exome/genome"}. Conversely, WES-based signatures (such as \code{"cosmic_v3_exome"}) assessed against WGS spectra need to be adjusted using \code{sig.bkg.adj="genome/exome"}. 
Other potential mismatches between sequencing sources of spectra and signatures may arise from analyzing targeted sequencing data, as well as when analyzing data from different organisms. A \code{sig.bkg.adj="custom"} option is provided to allow users to adjust signatures based on a custom data frame of scale factors across channels via \code{sig.bkg.adj.custom}.
}
\value{
An object with S3 class \code{"mutSigMapper"}.
\item{spectra}{Input data frame.}
\item{sig}{Input data frame.}
\item{method}{Input parameter.}  
\item{tri.counts.method}{Input parameter.}
\item{signature.cutoff}{Input parameter.}
\item{noise}{Input parameter.}
\item{neg.binom.size}{Input parameter.}
\item{n_rdm}{Input parameter.}
\item{weights}{List of weights (i.e. signature contributions) for each input sample.}
\item{weights_rdm}{List of random weight matrices (rows=random spectra, columns=signatures) for each input sample.}
\item{map_pval}{Matrix of significance p-values (rows=signatures, columns=samples).}
\item{counts}{For \code{method="MutationalPatterns"}, list of counts contributed across signatures for each input sample. This is the unnormalized version of \code{weights}.}
\item{counts_rdm}{For \code{method="MutationalPatterns"}, list of random count matrices (rows=random spectra, columns=signatures) for each input sample. This is the unnormalized version of \code{weights_rdm}.}
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\references{Rosenthal R et al.
\emph{DeconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution}, Genome Biol. 17, 31 (2016).

Blokzijl F et al. 
\emph{MutationalPatterns: comprehensive genome-wide analysis of mutational processes}, 
Genome Med. 10, 33 (2018).

COSMIC (\emph{Catalogue Of Somatic Mutations In Cancer}), https://cancer.sanger.ac.uk/cosmic/signatures/index.tt

Kucab JE et al. 
\emph{A Compendium of Mutational Signatures of Environmental Agents}, Cell 177, 821-836.e16 (2019).
}
\seealso{\code{\link{ref_bkg_exome}}, \code{\link{ref_bkg_genome}}, \code{\link{plotSpectraCaterpillar}}, \code{\link{plotSpectraHeatmap}}, \code{\link{export}}
}
\examples{
\dontshow{
data(spectra)
map = mutSigMapper(spectra[,1:2],n_rdm=5)
}
\donttest{
# Spectra from Blokzijl et al:
data(spectra)
map = mutSigMapper(spectra)
}
}

