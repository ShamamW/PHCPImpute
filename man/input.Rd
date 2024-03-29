% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation_functions.R
\name{input}
\alias{input}
\title{Preparation of input files for imputation}
\usage{
input(haps, genetic_map, ancient_tped, chr)
}
\arguments{
\item{haps}{Phased reference genomes in .haps format, which is the output of SHAPEIT2.}

\item{genetic_map}{genetic map for the reference genome with 3 columns: chromosome number, genetic distance in Morgans, physical position. The columns header should be: chr, genetic_distance, bp}

\item{ancient_tped}{pseudo-haplotype to impute in Plink tped format}

\item{chr}{the number of chromosome to impute}
}
\value{
an input for the imputation
}
\description{
This function is used as a nested function in phcp_impute()
}
\seealso{
\link{phcp_impute} which is the imputation function that uses this function
}
