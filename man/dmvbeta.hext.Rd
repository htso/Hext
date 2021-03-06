% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Density_fun.r
\name{dmvbeta.hext}
\alias{dmvbeta.hext}
\title{Multivariate beta emission density for hidden state j in an HMM}
\usage{
dmvbeta.hext(x, j, model)
}
\arguments{
\item{x}{numeric matrix of size T x p, where T is number of obs, p the dimension of x}

\item{j}{hidden state index}

\item{model}{HMM model}
}
\value{
vector of proabilities, length equal the number of rows in x
}
\description{
Multivariate beta emission density for hidden state j in an HMM
}

