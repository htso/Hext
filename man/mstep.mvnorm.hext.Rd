% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mstep.r
\name{mstep.mvnorm.hext}
\alias{mstep.mvnorm.hext}
\title{M-step re-estimation function for multivariate gaussian emission distribution, no regularization}
\usage{
mstep.mvnorm.hext(x, wt)
}
\arguments{
\item{x}{numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data}

\item{w}{weight matrix, gamma from EStep C function}

\item{model}{HMM model as defined in hmmspec.gen(..)}
}
\value{
emission slot of Hext object
}
\description{
Re-estimates the parameters of the multivariate gaussian emission distribution as part of the EM algorithm for HMM. No regularizer is used.
}
\author{
Horace W Tso
}

