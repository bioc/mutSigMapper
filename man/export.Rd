\name{export}
\alias{export}
\title{
generates plain text files from mutSigMapper object
}
\description{
This function enables the extraction of data  
from a \code{mutSigMapper} object. Plain text data files are produced with tab- or comma-separated-value formats.
}
\usage{
export(x, dest_dir=getwd(), dest_dir_create=TRUE, delim=c("tsv","csv"), 
input.data=TRUE)
}
\arguments{
  \item{x}{\code{mutSigMapper} object.}
  \item{dest_dir}{Destination directory. Default is the working directory.}
  \item{dest_dir_create}{Creates destination directory if it does not exist already. Default is \code{TRUE}.}
  \item{delim}{Delimiter for the generated files, either tab-separated (\code{"tsv"}) or comma-separated (\code{"csv"}) values. Default is \code{"tsv"}.}
  \item{input.data}{Logical variable indicating whether to generate files for the input signature and spectra data frames. Default is \code{TRUE}.}
}
\value{
None.
}
\author{Julian Candia\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{mutSigMapper}}}
\examples{
\dontshow{
data(spectra)
map = mutSigMapper(spectra[,1:2],n_rdm=5)
}
\donttest{
data(spectra)
map = mutSigMapper(spectra)
export(x=map,dest_dir=tempdir())
}
}

