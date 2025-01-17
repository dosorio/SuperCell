% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_knn_graph.R
\name{build_knn_graph}
\alias{build_knn_graph}
\title{Build kNN graph}
\usage{
build_knn_graph(
  X,
  k = 5,
  from = c("dist", "coordinates"),
  use.nn2 = TRUE,
  return_neighbors_order = F,
  dist_method = "euclidean",
  cor_method = "pearson",
  p = 2
)
}
\arguments{
\item{X}{either distance or matrix of coordinates (rows are samples and cols are coordinates)}

\item{k}{kNN parameter}

\item{from}{from which data type to build kNN network: "dist" if X is a distance (dissimilarity) or "coordinates" if X is a matrix with coordinates as cols and cells as rows}

\item{use.nn2}{whether use RANN::nn2 method to buid kNN network faster (avaivable only for "coordinates" option)}

\item{return_neighbors_order}{whether return order of neighbors (not available for nn2 option)}

\item{dist_method}{method to compute dist (if X is a matrix of coordinates) available: c("cor", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")}

\item{cor_method}{if distance is computed as correlation (dist_method == "cor), which type of correlation to use (available: "pearson", "kendall", "spearman")}

\item{p}{p param in \code{"dist"} function}
}
\value{
a list with components
\itemize{
  \item graph.knn - igraph object
  \item order - Nxk matrix with indecies of k nearest neighbors ordered by relevance (from 1st to k-th)
}
}
\description{
Build kNN graph either from distance (from == "dist") or from coordinates (from == "coordinates")
}
