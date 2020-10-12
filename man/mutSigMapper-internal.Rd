\name{mutSigMapper-internal}
\title{Internal mutSigMapper functions}
\alias{mutSigMapper_deconstructSigs}
\alias{mutSigMapper_MutationalPatterns}
\alias{print.mutSigMapper}
\description{Internal mutSigMapper functions}
\usage{
mutSigMapper_deconstructSigs(spectra, sig, n_rdm, noise, neg.binom.size, 
signature.cutoff, sig.bkg.adj)

mutSigMapper_MutationalPatterns(spectra, sig, n_rdm, noise, neg.binom.size, 
signature.cutoff, sig.bkg.adj)
}
\value{
See \code{\link{mutSigMapper}}.
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\details{These functions are not part of the package application programming interface (API) and are not intended for use by end users. Functions 
\code{mutSigMapper_deconstructSigs} and \code{mutSigMapper_MutationalPatterns} are called by function \code{mutSigMapper} to implement non-negative least squares regression procedures from packages 
\code{deconstructSigs} and \code{MutationalPatterns}, respectively.
}
\keyword{internal}
