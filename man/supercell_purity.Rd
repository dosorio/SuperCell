% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_purity.R
\name{supercell_purity}
\alias{supercell_purity}
\title{Compute purity of super-cells}
\usage{
supercell_purity(clusters, supercell_membership)
}
\arguments{
\item{clusters}{vector of clustering assignment (reference assignment)}

\item{supercell_membership}{vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)}
}
\value{
a vector of super-cell purity, which is defined as a proportion of the most abundant cluster within super-cell.
With 1 meaning that super-cell consists of single cells from one cluster (reference assignment)
}
\description{
Compute purity of super-cells
}
