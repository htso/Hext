% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mstep.r
\name{mstep.mvnorm.wish}
\alias{mstep.mvnorm.wish}
\title{M-step re-estimation function for multivariate gaussian emission distribution with inverse Wishart regularizer}
\usage{
mstep.mvnorm.wish(x, w, wdf = 1)
}
\arguments{
\item{x}{numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data (columns are features)}

\item{w}{weight matrix, size T x K, where each row is the weight for a time step, and columns are states.}

\item{model}{HMM model as defined in hmmspec.gen(..)}
}
\value{
emission slot of Hext object
}
\description{
Re-estimates the parameters of the multivariate normal emission distribution with inverse Wishart regularization.
}
\details{
Ref : Snoussi, Mohammad-Djafari (2005), Degeneracy and likelihood penalization
     in multivariate gaussian mixture models, preprint.  
     wdf : Wishart degree of freedom, i.e. wishart df = wdf + p, where p is the dimension of the data. Thus, wdf > -1.
}
\author{
Horace W Tso
}

