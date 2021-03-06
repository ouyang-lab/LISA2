% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{climbMount}
\alias{climbMount}
\title{Find nearest peak for each valley}
\usage{
climbMount(valley, peaks, dst, nid)
}
\arguments{
\item{dst}{vector, the density vector of each point}

\item{nid}{matrix, the knn matrix}
}
\value{
the vector of nearest peak for each valley
}
\description{
Find nearest peak for each valley
}
