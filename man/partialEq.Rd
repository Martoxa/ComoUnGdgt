% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partialEq.R
\name{partialEq}
\alias{partialEq}
\title{Generates partial equations for GDGT indices}
\usage{
partialEq(gdgt, eq)
}
\arguments{
\item{gdgt}{Data frame of peak areas or fractional abundances of GDGTs where each row is a sample and each column a GDGT with the appropriate name format.}

\item{eq}{Equation for the index that wants to be calculated. This is provided by the index function as a text string.}
}
\value{
The function returns the modified equation for the index.
}
\description{
Compares the provided list of GDGTs to the required list of GDGTs for the index that one wants to calculate. Any GDGTs not provided in the dataset are removed from the equation and a warning is generated.
}
