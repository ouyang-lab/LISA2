% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainMethods.R
\name{getCluCellTh}
\alias{getCluCellTh}
\title{getCluCellTh: get the main cell types in the cluster}
\usage{
getCluCellTh(cluid, cell_labels, th = 0.1, subN = 4)
}
\arguments{
\item{cluid}{cluster id}

\item{cell_labels}{cell annotations}

\item{th}{threshold for the proportion of the cell type in the cluster 0.1,}

\item{subN}{threshold of the number of character in the string}
}
\value{
main cell types in each cluster
}
\description{
getCluCellTh: get the main cell types in the cluster
}
\examples{

}
