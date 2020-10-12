\name{WGS_PCAWG}
\alias{WGS_PCAWG}
\docType{data}
\title{Whole genome sequencing datasets from PCAWG}
\description{
Data frame that contains the mutational counts of 2780 pan-cancer samples from 37 different cancer types (columns) along the 96 single-nucleotide-substitution trinucleotide contexts (rows) based on whole genome sequencing from the ICGC/TCGA Pan-Cancer Analysis of Whole Genomes (PCAWG) Network. 
}
\usage{data(WGS_PCAWG)}
\format{
A list that contains:
\itemize{
\item{\code{spectra}: Data frame with mutational counts for all samples (columns) across 96 trinucleotide contexts (rows).}
\item{\code{sample_study}: Data frame with sample names and the study they belong to. Rows in this data frame correspond to (and appear in the same order as) columns in the \code{spectra} data frame.}
}
}
\references{The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Network. Pan-cancer analysis of whole genomes. Nature 578, 82-93 (2019).

The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium, https://dcc.icgc.org/releases/PCAWG
}

\keyword{datasets}

