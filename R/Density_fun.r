#' Generic density function for HMM's emission probability
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data
#' @param j hidden state index
#' @param model current model as defined in hmmspec.gen(..)
#' @return the emssion density of hidden state j
#' @details see example code for usage
#' @export
generic.density = function( x, j, model, NegInf=-1000 ) {
   require(corpcor)
   p = NCOL(x)
   Ttot = NROW(x)
   D = model$D
   M = model$M
   if ( length(model$mtype) != p )
      stop("no. of columns does not equal length of mtype.")

   # 1 = discrete
   ix1 = which(model$mtype == 1)
   if ( length(ix1) == 1 ) {  
     # prob.v is a vector of length D, the number of multinomial outcomes
     prob.v = model$parms.emission$M1$prob.l[[j]]
     # x is coded so that it just take a simple table look-up to get the emission probability
     if ( model$Log ) {
        logb1 = log(prob.v[x])   # table look up
     } else {
        logb1 = prob.v[x]   # table look up     
     }
   } else if ( length(ix1) > 1 ) {
     # prob.m is a matrix max(D) x p
     prob.m = model$parms.emission$M1$prob.l[[j]]
     bb = matrix(0, nrow=Ttot, ncol=length(ix1))
     for ( i in 1:length(ix1) ) {
       if ( model$Log ) {
         bb[,i] = log(prob.m[x[, ix1[i]], i])  # <-- make sure usage matches definition
       } else {
         bb[,i] = prob.m[x[, ix1[i]], i]  # <-- make sure usage matches definition
       }
     }
     if ( model$Log ) {
       logb1 = apply(bb, 1, sum)
     } else {
       logb1 = apply(bb, 1, prod)
     }
   } else {
     if ( model$Log ) {
       logb1 = rep(0, Ttot)
     } else {
       logb1 = rep(1, Ttot)
     }
   }

   # 2 = normal/multivariate
   ix2 = which(model$mtype == 2)
   if ( length(ix2) >= 1 ) {
     # for multivariate normal, M2$mu is a list of vector where the ith element
     # of the list is a vector for state i, and each element of the vector
     # is the mu for a gaussian density;
     # M2$sigma is a list of matrices where the ith element is the cov matrix for
     # state i.
     mu.v = model$parms.emission$M2$mu.l[[j]]
     sigma.m = model$parms.emission$M2$sigma.l[[j]]
     #cat("In generic.density, i see sigma :", as.vector(sigma.m), "\n")
     #ev = eigen(sigma.m)$values
     #cat("In generic.density, i see eigenvalues :", ev, "\n")
     if ( model$Log == FALSE & model$Intv == TRUE ) {
       logb2 = pmvnorm.intv(x[,ix2, drop=FALSE], mean=mu.v, sigma=sigma.m)
     } else {
       # Note that Dmvnorm will return a log value if log is set to TRUE.
       logb2 = Dmvnorm(x[,ix2, drop=FALSE], mean=mu.v, sigma=sigma.m, log=model$Log)
     }
   } else {
     if ( model$Log ) {
       logb2 = rep(0, Ttot)
     } else {
       logb2 = rep(1, Ttot)
     }
   }

   # 5 = gaussian mixture
   ix5 = which(model$mtype == 5)
   #cat("ix5", ix5, "\n")
   if ( length(ix5) >= 1 ) {
     # M5$mu is a list of lists where each element of the lists is a list, and
     # each element of this sublist is a vector corresponding to the mu for each
     # component of the mixture gaussian density; likewise, M5$sigma is a list
     # of lists, each list corresponds to a state, and each element of the sublist
     # corresponds to the cov matrix of a mixture component;
     # M5$c.l is a list of vectors, each element is the state i's mixture proportion
     # NOTATION :   x.ll = list of list     x.l = list
     mu.l = model$parms.emission$M5$mu.ll[[j]]
     sigma.l = model$parms.emission$M5$sigma.ll[[j]]
     c.v = model$parms.emission$M5$c.l[[j]]
     Ttot = NROW(x)
     # check for consistency
     #if ( length(c.v) != M ) {
      # cat("c.v length, M :", length(c.v), M, "\n")
      # stop("length mismatch in mixture parameters.")
     #}
     # logbb : each row is a log probability of observing one observation (a row in x[,])
     # conditioned on mu and sigma of the hidden state
     logbb = matrix(NA, nrow=Ttot, ncol=M)
     for ( m in 1:M ) {
       tmp = x[, ix5, drop=FALSE]
       mu.v = mu.l[[m]]
       sigma.m = sigma.l[[m]]
       logbb[,m] = log(c.v[m]) + Dmvnorm(tmp, mean=mu.v, sigma=sigma.m, log=TRUE)
     }
     sum.p = rep(NA, Ttot)
     # since the columns of logbb are the gaussian mixture components, we need to sum them.
     for ( tt in 1:Ttot) {
       sum.p[tt] = Elnsum.vec(logbb[tt,]) # in the non-log version: sum.pr = rowSums(bb)
     }
     
     if ( model$Log == TRUE ) {
       logb5 = sum.p        
     } else {
       logb5 = exp(sum.p)
     }
   } else {
     if ( model$Log ) {
       logb5 = rep(0, Ttot)
     } else {
       logb5 = rep(1, Ttot)
     }
   }

   # 3 = gamma
   ix3 = which(model$mtype == 3)
   if ( length(ix3) >= 1 ) {
     shape.v = model$parms.emission$M3$shape.l[[j]]
     scale.v = model$parms.emission$M3$scale.l[[j]]
     bb = matrix(NA, ncol=length(ix3), nrow=Ttot)
     for ( i in 1:length(ix3) ) {
       tmp = dgamma(x[,ix3[i]], shape=shape.v[i], scale=scale.v[i])       
       #iz = which(tmp < 1e-323)
       #tmp[iz] = 1e-323
       #cat("counting small values in dgamma:", length(iz), "\n")
       if ( model$Log ) {
         bb[,i] = log(tmp)
       } else {
         bb[,i] = tmp
       }
     }
     if ( model$Log ) {
       logb3 = apply(bb, 1, sum)
     } else {
       logb3 = apply(bb, 1, prod)
     }
   } else {
     if ( model$Log ) {
       logb3 = rep(0, Ttot)
     } else {
       logb3 = rep(1, Ttot)
     }
   }

   # 4 = beta
   ix4 = which(model$mtype == 4)
   if ( length(ix4) >= 1 ) {
     shape1.v = model$parms.emission$M4$shape1.l[[j]]
     shape2.v = model$parms.emission$M4$shape2.l[[j]]
     bb = matrix(NA, ncol=length(ix4), nrow=Ttot)
     for ( i in 1:length(ix4) ) {
        tmp = dbeta(x[,ix4[i]], shape1=shape1.v[i], shape2=shape2.v[i])
        #iz = which(tmp < 1e-323)
        #tmp[iz] = 1e-323
        #cat("counting small values in dgamma:", length(iz), "\n")
        if ( model$Log ) {
          bb[,i] = log(tmp)
        } else {
          bb[,i] = tmp
        }
     }
     if ( model$Log ) {
       logb4 = apply(bb, 1, sum)
     } else {
       logb4 = apply(bb, 1, prod)
     }
   } else {
     if ( model$Log ) {
       logb4 = rep(0, Ttot)
     } else {
       logb4 = rep(1, Ttot)
     }
   }
   if ( model$Log ) {
     b = logb1 + logb2 + logb3 + logb4 + logb5
     # NOTE : b may have -Inf, which can't be passed into C, so need to make them very small
     b[which(b < NegInf)] = NegInf
   } else {
     b = logb1 * logb2 * logb3 * logb4 * logb5
   }   
   return(b)
}

#' Multivariable normal distribution at a point, calculate probability for a small interval around a point
#'
#' @param x numeric matrix
#' @param mean mean vector of the multivariate normal 
#' @param sigma covariance matrix
#' @param eps small numeric neighborhood around x
#' @return vector of proabilities, length equal the number of rows in x
#' @details used in generic.density
#' @export
pmvnorm.intv = function (x, mean, sigma, eps=0.1) {
  require(mvtnorm)
  ll = x - eps/2
  uu = x + eps/2
  pp = double(NROW(x))
  for ( i in 1:NROW(x)) {
    pp[i] = pmvnorm(lower=ll[i,], upper=uu[i,], mean=mean, sigma=sigma)    
  }
  return(pp)
}

#' Multivariable normal density for a hidden state evaluated at x
#'
#' @param x numeric matrix
#' @param model HMM model
#' @param j hidden state index
#' @param Log boolean whether log is used, default to FALSE
#' @return vector of proabilities, length equal the number of rows in x
#' @details used in generic.density
#' @export
dmvnorm.hext = function (x, j, model, Log=FALSE) {
  mu.v = model$parms.emission$M2$mu.l[[j]]
  sigma.m = model$parms.emission$M2$sigma.l[[j]]
  ans = Dmvnorm(x, mean=mu.v, sigma=sigma.m, log=Log)
  #ans[is.na(ans)] <- 1   # not sure what this is for?
  return(ans)
}

#' Multivariable normal density, similar to the similar named function in mvtnorm
#'
#' @param x numeric matrix
#' @param mean mean vector
#' @param sigma covariance matrix
#' @param log boolean whether log is used, default to FALSE
#' @param num.cutoff numeric cutoff value
#' @return vector of proabilities, length equal the number of rows in x
#' @details use ginv from MASS to compute inverse matrix
#' @export
Dmvnorm = function (x, mean, sigma, log=TRUE, num.cutoff=1e-100)  {
    require(MASS)  # MASS is needed because of ginv
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    #if (missing(mean)) {
    #    mean <- rep(0, length = NCOL(x))
    #}
    #if (missing(sigma)) {
    #    sigma <- diag(NCOL(x))
    #}
    #if (NCOL(x) != NCOL(sigma)) {
    #    stop("x and sigma have non-conforming size")
    #}
    #if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
    #   assign("SIGMA", sigma, envir = .GlobalEnv)
    #   stop("sigma must be a symmetric matrix")
    #}
    #if (length(mean) != NROW(sigma)) {
    #    stop("mean and sigma have non-conforming size")
    #}
    # Use this function from MASS, instead of the default solve() 
    inv.sigma = ginv(sigma)    
    
    distval = mahalanobis(x, center=mean, cov=inv.sigma, inverted=TRUE)
    ev = eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
    logdet = sum(log(ev))
    
    logretval = -(NCOL(x) * log(2 * pi) + logdet + distval)/2
    ix = which(logretval > 0) # because exp(logretval) would be > 1
    #if ( length(ix) > 0 ) {
    #  cat("NCOL(x):", NCOL(x), " ")
    #  cat("logdet:", logdet, " ")
    #  cat("distval[ix]:", distval[ix], "\n")
    #}
    
    if (log == TRUE) { 
      return(logretval)
    } else { 
      # put a floor on the probability value
      tmp = exp(logretval)
      ix = which( tmp < num.cutoff )
      if ( length(ix) > 0 ) {
        tmp[ix] = num.cutoff
      } 
      return(tmp)
    }
}

#' Beta emission density for hidden state j in an HMM
#'
#' @param x numeric vector
#' @param j hidden state index
#' @param model HMM model
#' @return vector of proabilities, length equal the number of rows in x
#' @export
dbeta.hext = function(x, j, model) {
    b = dbeta(x, shape1=model$parms.emission$shape1[j],
                 shape2=model$parms.emission$shape2[j])
    b[is.na(b)] = 1
    return(b)
}

#' Multivariate beta emission density for hidden state j in an HMM
#'
#' @param x numeric matrix of size T x p, where T is number of obs, p the dimension of x
#' @param j hidden state index
#' @param model HMM model
#' @return vector of proabilities, length equal the number of rows in x
#' @export
dmvbeta.hext = function(x, j, model) {
    p = NCOL(x)
    n = NROW(x)
    Shape1 = model$parms.emission$shape1[[j]]
    Shape2 = model$parms.emission$shape2[[j]]
    # shape and scale should be vector of length p
    if ( length(Shape1) != p | length(Shape2) != p )
      stop("shape1 or shape2 are of the wrong length!")
    b = matrix(NA, nrow=n, ncol=p)
    for ( i in 1:p ) {
       b[,i] = dbeta(x[,i], shape1=Shape1[i], shape2=Shape2[i])
    }
    ans = apply(b, 1, prod)
    return(ans)
}

#' Univariate gamma emission density for hidden state j in an HMM
#'
#' @param x numeric vector
#' @param j hidden state index
#' @param model HMM model
#' @return vector of proabilities, length equal the number of rows in x
#' @export
dgamma.hext = function(x, j, model) {
    b = dgamma(x, shape=model$parms.emission$shape[j],
                  scale=model$parms.emission$scale[j])
    b[is.na(b)] = 1
    return(b)
}

#' Multivariate gamma emission density for hidden state j in an HMM
#'
#' @param x numeric matrix of size T x p, where T is number of obs, p the dimension of x
#' @param j hidden state index
#' @param model HMM model
#' @return vector of proabilities, length equal the number of rows in x
#' @export
dmvgamma.hext = function(x, j, model) {
    p = NCOL(x)
    n = NROW(x)
    Shape = model$parms.emission$shape[[j]]
    Scale = model$parms.emission$scale[[j]]
    # shape and scale should be vector of length p
    if ( length(Shape) != p | length(Scale) != p )
      stop("shape or scale are of the wrong length!")
    b = matrix(NA, nrow=n, ncol=p)
    for ( i in 1:p ) {
       b[,i] = dgamma(x[,i], shape=Shape[i], scale=Scale[i])
    }
    ans = apply(b, 1, prod)
    return(ans)
}

#' Weighted covariance matrices, similar to cov.wt in hsmm
#'
#' @param x numeric matrix of size T x p, where T is number of obs, p the dimension of x
#' @param wt vector of weights for each obs
#' @param cor boolean whether correlation matrix to be returned
#' @param center centers to be used when computing covariances
#' @param method either "unbiased" or "ML
#' @return list with components cov, center, n.obs, wt, cor
#' @export
cov.wt.hext = function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE, method = c("unbiased", "ML")) {
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    else if (!is.matrix(x)) 
        stop("'x' must be a matrix or a data frame")
    if (!all(is.finite(x))) 
        stop("'x' must contain finite values only")
    n <- nrow(x)
    if (with.wt <- !missing(wt)) {
        if (length(wt) != n) 
            stop("length of 'wt' must equal the number of rows in 'x'")
        s = sum(wt)
        if (any(wt < 0) || s == 0) 
            stop("weights must be non-negative and not all zero")
        wt <- wt/s
    
    }
    if (is.logical(center)) {
        center <- if (center) 
            colSums(wt * x)
        else 0
    }
    else {
        if (length(center) != ncol(x)) 
            stop("length of 'center' must equal the number of columns in 'x'")
    }
    x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
    denom = 1 - sum(wt^2)
    if ( denom < 1e-300 ) {
      stop("1 - sum(wt^2) is very close to zero!\n")
    }
    cov <- switch(match.arg(method), unbiased=crossprod(x)/denom, ML = crossprod(x))
    y <- list(cov = cov, center = center, n.obs = n)
    if (with.wt) 
        y$wt <- wt
    if (cor) {
        Is <- 1/sqrt(diag(cov))
        R <- cov
        R[] <- Is * cov * rep(Is, each = nrow(cov))
        y$cor <- R
    }
    y
}

