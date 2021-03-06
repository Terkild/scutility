% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color_functions.R
\name{color_subcluster}
\alias{color_subcluster}
\title{Expand color vector for subclusters}
\usage{
color_subcluster(
  subcluster,
  cluster,
  cluster_colors = c(),
  order = TRUE,
  brightness_change = NA,
  max_change = 1
)
}
\arguments{
\item{subcluster}{vector of subcluster assignments (can be duplicated across cell assignments)}

\item{cluster}{vector of cluster assignment for each subcluster (vector of same length as subcluster)}

\item{cluster_colors}{(named) vector of colors for each cluster}

\item{order}{should colors be ordered by subcluster size? Set to FALSE to use factor levels}

\item{brightness_change}{how big a change in color brightness should be applied at each step (between 0 and 1)}
}
\value{
vector of colors
}
\description{
Make a vector of colors for a vector of subclusters grouped by their
assignment to a cluster
}
