% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{Tw2}
\alias{Tw2}
\title{Calculate the Tw2 Statistic for Heteroscedastic Test}
\usage{
Tw2(dm, f)
}
\arguments{
\item{dm}{A distance matrix, representing dissimilarity between observations.}

\item{f}{A factor variable indicating the group for each observation.}
}
\value{
The calculated Tw2 statistic.
}
\description{
This function calculates the Tw2 statistic to compare means of k populations
with potentially unequal observations. It is suitable for microbiome data
and utilizes permutation testing for significance estimation.
}
\details{
The function first checks if the factor variable has exactly two levels.
The Tw2 statistic is a modification of Hotelling's T-square statistic
adapted for heteroscedasticity and specifically suitable for microbiome data.
It calculates the sum of squares within each group and then computes the
Tw2 statistic based on these sum of squares.
}
\examples{
\dontrun{
# Generate a synthetic distance matrix and a factor variable
dm <- matrix(runif(100), nrow = 10)
f <- factor(rep(1:2, each = 5))
# Calculate the Tw2 statistic
Tw2_stat <- Tw2(dm, f)
}

}
\references{
Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
}
