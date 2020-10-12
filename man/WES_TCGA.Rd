\name{WES_TCGA}
\alias{WES_TCGA}
\docType{data}
\title{Whole exome sequencing datasets from TCGA}
\description{
Data frame that contains the mutational counts of 9493 pan-cancer samples from 32 different cancer types (columns) along the 96 single-nucleotide-substitution trinucleotide contexts (rows) based on whole exome sequencing from the Cancer Genome Project (TCGA).
}
\usage{data(WES_TCGA)}
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

