% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IBT.R
\name{IBT}
\alias{IBT}
\title{Isomerization of breanched tetraethers index}
\usage{
IBT(gdgt, complete = TRUE, out = "all", na.ignore = FALSE)
}
\arguments{
\item{gdgt}{Data frame with the peak areas or fractional abundances of brGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.}

\item{complete}{Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.}

\item{na.ignore}{Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.}
}
\value{
The function returns the BIT value for each sample provided. Any value that blows up (<-10000 or >10000) is turned into NA.
}
\description{
Calculates the IBT index (Ding, 2015) for brGDGTs presented as either peak areas or fractional abundances.
}
\references{
Ding, S., Xu, Y., Wang, Y., He, Y., Hou, J., Chen, L., & He, J. S. (2015). Distribution of branched glycerol dialkyl glycerol tetraethers in surface soils of the Qinghai–Tibetan Plateau: implications of brGDGTs-based proxies in cold and dry regions. Biogeosciences, 12(11), 3141-3151.
}
