% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.r
\name{phenograph_clustering}
\alias{phenograph_clustering}
\title{PhenoGraph wrapper function}
\usage{
phenograph_clustering(data = NULL, k = 30, seed = 999)
}
\arguments{
\item{data}{A data frame: samples in rows, features in columns}

\item{k}{An integer of the number of neighbors}

\item{seed}{An integer of seed for reproducibility}
}
\value{
A data frame which appends the cluster assignment to input data \code{data}
}
\description{
PhenoGraph wrapper function
}
