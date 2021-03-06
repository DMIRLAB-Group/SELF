% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{randomGraph}
\alias{randomGraph}
\title{Generate a random graph}
\usage{
randomGraph(dim, indegree, maxite = 10000)
}
\arguments{
\item{dim}{The dimension of the random graph}

\item{indegree}{The average indegree of random graph for each nodes}

\item{maxite}{The maximum iterations to find the random graph}
}
\value{
Return a random graph
}
\description{
Generate a random graph based on the given dimension size and average indegree
}
\examples{
randomGraph(dim=10,indegree=1)
}
