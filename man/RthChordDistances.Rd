% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RthChordDistances.R
\name{RthChordDistances}
\alias{RthChordDistances}
\title{Earth Chord Distances}
\usage{
RthChordDistances(llPoints1, llPoints2)
}
\arguments{
\item{llPoints1}{Nx2 matrix where each row is a longitude-latitude pair (-180 to 180, and -90 to 90, respectively).}

\item{llPoints2}{Mx2 matrix where each row is a longitude-latitude pair (-180 to 180, and -90 to 90, respectively).}
}
\value{
A NxM matrix of chordal distance between the ith row of llPoint1 and the jth row of llPoint2. Radius of the earth is taken as 6378.137
}
\description{
Calculates the chordal distances (in km) between points 1 and 2.
}
\author{
Jessica Tierney https://github.com/jesstierney/BAYSPAR
}
