% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEX_forward.R
\name{TEX_forward}
\alias{TEX_forward}
\title{BaySPR TEX86 forward model}
\usage{
TEX_forward(lat, lon, t, type = "SST", mode = "standard", tolerance = NA)
}
\arguments{
\item{lat}{Latitude from -90 to 90}

\item{lon}{Longitude from -180 to 180}

\item{t}{A vector (N x 1) of SST or subT}

\item{type}{Selection of either "SST" (default) or "subT" for the model to use.}

\item{mode}{A string specifying "standard" or "analog".}

\item{tolerance}{Search tolerance in temperature units}
}
\value{
A N x 1000 ensemble of possible temperature values.
}
\description{
This function forward models TEX86 from either SSTs or subTs.
}
\references{
Tierney, J. E., & Tingley, M. P. (2014). A Bayesian, spatially-varying calibration model for the TEX86 proxy. Geochimica et Cosmochimica Acta, 127, 83-106.
}
