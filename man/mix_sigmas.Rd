% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post.R
\name{mix_sigmas}
\alias{mix_sigmas}
\title{Extract Standard Deviations of Mixture Components}
\usage{
mix_sigmas(x)
}
\arguments{
\item{x}{A mixture distributional object}
}
\value{
numeric or list object
}
\description{
Extract Standard Deviations of Mixture Components
}
\details{
If a distributional object that is a mixture of two or more normal
distributions is read in, the function will return a numeric object with
the standard deviations of each normal component. If the distributional
object is a mixture of two or more multivariate normal distributions, the
function will return a list with the covariance matrices of each multivariate
normal component.
}
\examples{
library(distributional)
mix_norm <- dist_mixture(comp1 = dist_normal(1, 10),
                         comp2 = dist_normal(1.5, 12),
                         weights = c(.5, .5))
mix_sigmas(mix_norm)
}
