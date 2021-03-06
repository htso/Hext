% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mstep.r
\name{mstep.gamma}
\alias{mstep.gamma}
\title{M-step re-estimation function for univariate gamma emission distribution}
\usage{
mstep.gamma(x, w)
}
\arguments{
\item{x}{numeric vector}

\item{w}{weight matrix, size T x K}

\item{model}{HMM model as defined in hmmspec.gen(..)}
}
\value{
emission slot of Hext object
}
\description{
Re-estimates the parameters of univariate gamma emission distribution as part of the EM algorithm for HMM. No regularizer is used.
}
\details{
Needs the solnp function in the Rsolnp package.
}
\author{
Horace W Tso
}

