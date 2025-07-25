% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_utils.R
\name{dist.group.sigma2}
\alias{dist.group.sigma2}
\title{Calculate Group-wise Sigma Squared}
\usage{
dist.group.sigma2(dm, f)
}
\arguments{
\item{dm}{Distance matrix.}

\item{f}{Factor variable for group definition.}
}
\value{
A diagonal matrix of group-wise sigma squared values.
}
\description{
This function calculates the sigma squared statistic for each group defined
by a factor variable in a given distance matrix.
}
\examples{
\dontrun{
dm <- as.dist(matrix(runif(100), nrow = 10))
f <- factor(c(rep("A", 5), rep("B", 5)))
dist.group.sigma2(dm, f)
}
}
