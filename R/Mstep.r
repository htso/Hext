#' Generic M-step function to re-estimate the paramters of an emission probability distribution.
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data
#' @param wt weight matrix
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return the five types of M-step results
#' @description Re-estimate the parameters of an emission distribution as part of the EM algorithm for HMM. 
#'    This is called by the hext function.
#' @usage mstep.generic(x, wt, model)    
#' @details see example code for usage
#' @author Horace W Tso
#' @export
mstep.generic = function(x, wt, model) {
  # type codes (mtype):
  #   discrete/binomial = 1
  #   multivariate normal = 2
  #   gamma = 3
  #   beta = 4
  #   gaussian mixture = 5
  D = model$D
  mtype = model$mtype
  ix1 = which(mtype == 1)
  if ( length(ix1) == 1 ) {
    # if there is only one column in obs, then no need to use ix1
    if ( length(mtype) == 1 ) {
      M1 = mstep.disc(x, wt, D)
    } else { 
      # this is the case where there is only one discrete column
      # and one or more other types of observation, thus need ix1.
      M1 = mstep.disc(x[,ix1, drop=FALSE], wt, D)
    }
  } else if ( length(ix1) > 1 ) {
    M1 = mstep.mvdisc(x[,ix1], wt, D)
  } else {
    M1 = list(prob.l=list(), labels=NULL)
  }

  ix2 = which(mtype == 2)
  if ( length(ix2) == 1 ) {
    M2 = mstep.norm.hext(x[,ix2, drop=FALSE], wt)
  } else if ( length(ix2) > 1 ) {
    # =================================
    M2 = mstep.mvnorm.wish(x[,ix2], wt)
    # =================================
  } else {
    M2 = list(mu.l=list(), sigma.l=list())
  }

  ix3 = which(mtype == 3)
  if ( length(ix3) == 1 ) {
    M3 = mstep.gamma(x[,ix3, drop=FALSE], wt)
  } else if ( length(ix3) > 1 ) {
    M3 = mstep.mvgamma(x[,ix3], wt)
  } else {
    M3 = list(shape.l=list(), scale.l=list())
  }

  ix4 = which(mtype == 4)
  if ( length(ix4) == 1 ) {
    M4 = mstep.beta(x[,ix4, drop=FALSE], wt)
  } else if ( length(ix4) > 1 ) {
    M4 = mstep.mvbeta(x[,ix4], wt)
  } else {
    M4 = list(shape1.l=list(), shape2.l=list())
  }

  ix5 = which(mtype == 5)
  if ( length(ix5) > 0 ) {
    #cat("ix5:", ix5, "\n")
    # ====================================================================
    M5 = mstep.gaussian.mix(x[,ix5, drop=FALSE], wt, model)
    # ====================================================================
  } else {
    M5 = list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K))
  }
  return(list(M1=M1, M2=M2, M3=M3, M4=M4, M5=M5))
}

#' M-step re-estimation function for Gaussian Mixture emission distribution with inverse Wishart regularizer
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data
#' @param w weight matrix, gamma from EStep C function
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of the gaussian mixture emission distribution as part of the EM algorithm for HMM. 
#' @details Ref : Snoussi, Mohammad-Djafari (2005), Degeneracy and likelihood penalization
#'      in multivariate gaussian mixture models, preprint.  
#'      NOTE : gam.jk is *not* in log domain.      
#' @author Horace W Tso
#' @export
mstep.gaussian.mix = function (x, w, model) {
    idx = apply(is.na(x), 1, any)
    x = x[!idx, , drop = FALSE]
    w = w[!idx, , drop = FALSE]
    K = ncol(w)  # should agree with model$K
    Ttot = NROW(x) # = sum(model$TT)
    p = NCOL(x)  # = model$p
    M = model$M
    wdf = model$wdf
    
    if ( any(is.na(w)) | any(is.infinite(w))) {
      stop("na/infty found in w[].")
    }
    if ( all(w < 1e-100) == TRUE ) {
      stop("all w[] are near zero.")
    }
    w.tot = colSums(w)  # denominator for the reestimation equations, 1 x K
    
    emission = list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K))
    for ( i in 1:K ) {
      emission$c.l[[i]] = rep(NA, M)
    }
    
    gam.jk = array(NA, dim=c(Ttot, K, M))
    for ( i in 1:K )  {
       mu.l = model$parms.emission$M5$mu.ll[[i]]
       #if ( length(mu.l) != M ) cat("length of mu.l <> M\n")
       sigma.l = model$parms.emission$M5$sigma.ll[[i]]
       #if ( any(is.na(sigma.l))) {
        # stop("na in sigma.l:")
       #}
       #if ( length(sigma.l) != M ) cat("length of sigma.l <> M\n")
       c.v = model$parms.emission$M5$c.l[[i]]  # c.v is a vector of length M
       #if ( length(c.v) != M ) cat("length of c.v <> M\n")
       # Compute the conditional probabilities of the gaussian components in log domain :
       logbb = matrix(NA, nrow=Ttot, ncol=M)
       for ( m in 1:M ) {
         logbb[,m] = log(c.v[m]) + Dmvnorm(x, mean=mu.l[[m]], sigma=sigma.l[[m]], log=TRUE)
       }
       # Sum the contribution from the individual gaussian components
       # <==> sum by rows ==> a vector of length Ttot
       denom = double(Ttot)
       for ( tt in 1:Ttot) {
         denom[tt] = Elnsum.vec(logbb[tt,]) # in the non-log version: denom = rowSums(bb)
       }
       for ( m in 1:M ) {
          tmp = logbb[,m] + log(w[,i]) - denom
          gam.jk[,i,m] = exp(tmp)
       }
       if ( any(is.na(gam.jk[,i,m])) | any(is.infinite(gam.jk[,i,m]))) {
         cat("na/inf in gam.jk: i, m ", i, m)
         stop("...")
       }
    }
    
    for ( i in 1:K ) {
      # sum over mixture state k and time t, for every state i
      denom.j = sum(gam.jk[,i,])  # sum over mixture and sum over time, a scalar
      # cat("denom.j:", denom.j, "\n")
      for ( k in 1:M ) {
        if ( all(gam.jk[,i,k] < 1e-100) == FALSE ) {
          # c_jk =======================================================================
          emission$c.l[[i]][k] = sum(gam.jk[,i,k]) / denom.j  # sum over t; scalar / scalar = scalar
          #if ( any(is.na(emission$c.l[[i]][k])) ) {
          #  warning("na in c.l at m-step.")
          #}
          # mu_jk ====================================================================================
          
          xx = gam.jk[,i,k] * x   # (1 x T) (*) (T x p) = T x p, component by component multiplication
          denom.jk = sum(gam.jk[,i,k])  # summing over t for every i and k, result is a scalar
          mu.bar = colSums(xx) / denom.jk  # (1 x p) / scalar = 1 x p
          #if ( any(is.na(mu.bar)) ) {
          #  warning("na in mu.bar at m-step.")
          #}
          emission$mu.ll[[i]][[k]] = mu.bar
          
          # sigma_jk =============================================================================================
          x.centered = sweep(x, 2, mu.bar)  # T x p, 'sweep'' subtracts mu.bar ( 1 x p vector) from each row of x,
                                            # in effect, it is subtracting the mean from each column of
                                            # the multivariate observation vector (matrix), or (O_t - mu_jk)
          mat.sum = matrix(0, nrow=p, ncol=p)                                 
          for ( tt in 1:Ttot ) {
            crossp = x.centered[tt,] %*% t(x.centered[tt,]) # x %*% t(x) = p x p matrix
            mat.w = crossp * gam.jk[tt,i,k]  # p x p matrix * scalar = p x p
            mat.sum = mat.sum + mat.w   # element by element addition : p x p
          }
          # Regularization with Inverse Wishart distribution ======================
          # Fraser, p. 55, regularized covariance matrix.
          # see p. 8 of Snoussi(2005), where nu = p + wdf, in numerator
          # nu + n + 1 = (p + wdf) + p + 1 = wdf + 2*p + 1, in denominator
          cov.mat = ( mat.sum + (p + wdf)*diag(p) ) / ( denom.jk + wdf + 2*p + 1 )
          # =======================================================================
          # cat("na in cov.mat:", sum(is.na(cov.mat)), "\n")
          # =============================================================
          # matrix's columns and rows must have the same names, otherwise
          # isSymmetric in dmvnorm complains
          rownames(cov.mat) = paste("X", 1:p, sep="")
          colnames(cov.mat) = rownames(cov.mat)
          #if ( any(is.na(cov.mat)) ) {
          #  warning("na in cov.mat at m-step.")
          #}
          emission$sigma.ll[[i]][[k]] = cov.mat 
        } else {
          emission$c.l[[i]][k] = 0
          emission$mu.ll[[i]][[k]] = rep(0, p)
          emission$sigma.ll[[i]][[k]] = diag(p)
        }
      }
    }
    emission
}

#' M-step re-estimation function for multivariate gaussian emission distribution with inverse Wishart regularizer
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data (columns are features)
#' @param w weight matrix, size T x K, where each row is the weight for a time step, and columns are states.
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of the multivariate normal emission distribution with inverse Wishart regularization.
#' @details Ref : Snoussi, Mohammad-Djafari (2005), Degeneracy and likelihood penalization
#'      in multivariate gaussian mixture models, preprint.  
#'      wdf : Wishart degree of freedom, i.e. wishart df = wdf + p, where p is the dimension of the data. Thus, wdf > -1.
#' @author Horace W Tso
#' @export
mstep.mvnorm.wish = function (x, w, wdf=1) {
    ix = apply(is.na(x), 1, any)
    x = x[!ix, , drop = FALSE]
    w = w[!ix, , drop = FALSE]
    K = ncol(w)
    TT = NROW(x)
    p = NCOL(x)
    
    w.tot = colSums(w)  # denominator for the reestimation equations, 1 x K
    #if ( any(w.tot < 1e-100) ) {
      # cat("w.tot :", w.tot, "\n") 
      #stop("w.tot too small!")
    #}
    emission = list(mu.l=list(), sigma.l=list())
    mat.sum = matrix(0, nrow=p, ncol=p)
    for ( i in 1:K ) {
        xx = w[,i] * x   # (1 x T) x (T x p) = T x p
        mu.bar = colSums(xx) / w.tot[i]  # (1 x p) x scalar = 1 x p
        x.centered = sweep(x, 2, mu.bar)  # T x p        
        for ( tt in 1:TT ) {
          crossp = x.centered[tt,] %*% t(x.centered[tt,]) 
          mat.w = crossp * w[tt, i]  # p x p matrix times a scalar  = p x p
          mat.sum = mat.sum + mat.w   # element by element addition : p x p
        }
        # Regularization using the Inverse Wishart distribution ==============
        # Ref : Snoussi (2005), p. 8, where the inverse scale matrix is to be taken
        # as a diagonal matrix, wish.df is nu, the degree of freedom, p is the 
        # dimension of the data n.
        cov.mat = (mat.sum + (p + wdf)*diag(p)) / ( w.tot[i] + wdf + 2*p + 1 )
        # ====================================================================
        # matrix's columns and rows must have the same names, otherwise
        # isSymmetric() complains 
        rownames(cov.mat) = paste("X", 1:p, sep="")
        colnames(cov.mat) = rownames(cov.mat)
        emission$mu.l[[i]] = mu.bar
        emission$sigma.l[[i]] = cov.mat 
    }
    emission
}

#' M-step re-estimation function for multivariate gaussian emission distribution, no regularization
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data
#' @param w weight matrix, gamma from EStep C function
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of the multivariate gaussian emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @author Horace W Tso
#' @export
mstep.mvnorm.hext = function (x, wt) {
    ix = apply(is.na(x), 1, any)
    x = x[!ix, , drop = FALSE]
    wt = wt[!ix, , drop = FALSE]
    emission = list(mu.l=list(), sigma.l=list())
    for (i in 1:ncol(wt)) {
        tmp <- cov.wt(x, wt[,i])
        emission$mu.l[[i]] <- tmp$center
        emission$sigma.l[[i]] <- tmp$cov
    }
    emission
}

#' M-step re-estimation function for univarite gaussian emission distribution, no regularization
#'
#' @param x numeric matrix, size T x p, where T is the number of observations, p the dimension of the obs data
#' @param w weight matrix, gamma from EStep C function
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of the univariate gaussian emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @author Horace W Tso
#' @export
mstep.norm.hext = function (x, wt) {
    k = ncol(wt)
    mu = numeric(k)
    sigma = numeric(k)
    for (i in 1:k) {
        tmp = cov.wt(data.frame(x[!is.na(x)]), wt[!is.na(x),i])
        mu[i] = tmp$center
        sigma[i] = tmp$cov
    }
    list(mu.l = mu, sigma.l = sigma)
}

#' M-step re-estimation function for univariate discrete emission distribution
#'
#' @param x numeric matrix, size sum(T) x 1. All sequences stacked onto one column
#' @param w weight matrix, size sum(T) x K, gamma(i, t)
#' @param D number of alphabets in x
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of discrete emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @author Horace W Tso
#' @export
mstep.disc = function(x, w, D) {
  K = NCOL(wt)
  TT = NROW(x)
  w.sum = colSums(w)  # a vector of length K
  if ( length(which(w.sum < 1e-300)) > 0 ) {
    stop("some of the wt column sums are too close to zero.")
  }
  # prob.l is a list of matrices
  prob.l = list()
  # In this univariate case, b.vec is really just a vector; but to be consistent with
  # the usage in generic.density, i have to make it a single column matrix.
  b.vec = matrix(0, nrow=D, ncol=1)
  for ( j in 1:K ) {
    for ( fa in 1:D ) { 
      ix = which( x == fa )
      b.vec[fa, 1] = sum(w[ix, j])
    }  
    prob.l[[j]] = b.vec / w.sum[j]
  }
  list(prob.l=prob.l, labels=NULL)
}

#' M-step re-estimation function for multivariate discrete emission distribution
#'
#' @param x numeric matrix, size sum(T) x p
#' @param w weight matrix, size T x K
#' @param D number of alphabets in x, a vector of length equal to ncol(x)
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of discrete emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @author Horace W Tso
#' @export
mstep.mvdisc = function(x, w, D)  {
  K = NCOL(wt)
  TT = NROW(x)
  p = NCOL(x)
  wt.sum = colSums(wt)  # a vector of length K
  if ( length(which(wt.sum < 1e-300)) > 0 ) {
    stop("some of the wt column sums are too close to zero.")
  }
  prob.l = list()
  lookup.mat = matrix(0, nrow=max(D), ncol=p)
  for ( j in 1:K ) {
    for ( k in 1:p ) {
      for ( fa in 1:D[k] ) {
        # IMPORTANT NOTE : here i assume the discrete obs take on numeric value
        # from 1 to D[k] for column k of the obs. The next line looks for all
        # occurences of the factor (fa) in obs x[,k],
        ix = which( x[,k] == fa )
        # then sum up the soft weights, i.e. probability of the occurence of fa
        # in x for state j.
        lookup.mat[fa, k] = sum(wt[ix, j])
        # RECALL : wt = P( x_t=fa | s_t=j )
        # it has as many rows as x and as many columns as the no of states
      }
    }
    prob.l[[j]] = lookup.mat / wt.sum[j]
  }   
  list(prob.l=prob.l, labels=NULL)
}

#' M-step re-estimation function for univariate beta emission distribution
#'
#' @param x numeric matrix, size sum(T) x p
#' @param w weight matrix, size T x K
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of univariate beta emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @details Needs the solnp function in the Rsolnp package.
#' @author Horace W Tso
#' @export
mstep.beta = function(x, w) {
  require(Rsolnp)
  ix = which(is.na(x))
  if ( length(ix) > 0 ) {
    x = x[!ix]
    w = w[!ix, ,drop = FALSE]
  }
  para.start = c(runif(1,min=0.5,max=7), runif(1,min=0.5,max=7))   # c(shape1, shape2)
  lower.b = c(0.01, 0.01)
  upper.b = c(Inf, Inf)
  K = ncol(wt)
  Shape1 = list()
  Shape2 = list()
  i = 1
  nfail = 0
  while ( i <= K ) {
    sol = try(solnp(para.start, Qaux.beta, x=x, wt=w[,i],
             LB=lower.b, UB=upper.b,
             control=list(trace=FALSE, inner.iter=1000, outer.iter=1000,
             tol=1e-13, delta=1e-13)), FALSE)
    if ( !inherits(sol, "try-error") ) {
      Shape1[[i]] = sol$par[1]
      Shape2[[i]] = sol$par[2]
      i = i + 1
    } else {
      nfail = nfail + 1
      if ( nfail >= 5 ) {
        stop("Too many problems with solnp!")
      }
    }
  }
  list(shape1.l=Shape1, shape2.l=Shape2)
}

#' M-step re-estimation function for multivariate beta emission distribution
#'
#' @param x numeric matrix, size sum(T) x p
#' @param w weight matrix, size T x K
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of multivariate beta emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @details Needs the solnp function in the Rsolnp package.
#' @author Horace W Tso
#' @export
mstep.mvbeta = function(x, wt) {
  require(Rsolnp)
  p = NCOL(x)
  K = NCOL(wt)
  ix = apply(is.na(x), 1, any)
  if ( length(ix) > 0 ) {
    x = x[!ix, ,drop = FALSE]
    wt = wt[!ix, ,drop = FALSE]
  }
  emission = list(shape1.l=list(), shape2.l=list())
  lower.b = c(0.01, 0.01)
  upper.b = c(Inf, Inf)
  h1.mat = matrix(NA, ncol=K, nrow=p)
  h2.mat = matrix(NA, ncol=K, nrow=p)
  for ( j in 1:p ) {
    i = 1
    nfail = 0
    while ( i <= K ) {
      para.start = c(runif(1,min=0.5,max=7), runif(1,min=0.5,max=7))
      sol <- try(solnp(para.start, Qaux.beta, x=x[,j], wt=wt[,i],
               LB=lower.b, UB=upper.b,
               control=list(trace=FALSE, inner.iter=1000, outer.iter=1000,
               tol=1e-13, delta=1e-13)), FALSE)
      if (!inherits(sol,"try-error")) {
        h1.mat[j,i] = sol$par[1]
        h2.mat[j,i] = sol$par[2]
        i = i + 1
      } else {
        nfail = nfail + 1
        if ( nfail >= 5 ) {
          stop("Too many problems with solnp!")
        }
      }
    }
  }
  for ( i in 1:K ) {
    # columns are states
    emission$shape1.l[[i]] = h1.mat[,i]
    emission$shape2.l[[i]] = h2.mat[,i]
  }
  return(emission)
}

#' M-step re-estimation function for univariate gamma emission distribution
#'
#' @param x numeric vector
#' @param w weight matrix, size T x K
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of univariate gamma emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @details Needs the solnp function in the Rsolnp package.
#' @author Horace W Tso
#' @export
mstep.gamma = function(x, w) {
  require(Rsolnp)
  ix = which(is.na(x))
  if ( length(ix) > 0 ) {
    x = x[!ix]
    w = w[!ix, ,drop = FALSE]
  }
  lower.b = c(0.01, 0.01)
  upper.b = c(Inf, Inf)
  K = ncol(w)
  Shape = list()
  Scale = list()
  i = 1
  nfail = 0
  while( i <= K ) {
    para.start = c(runif(1,min=0.5,max=10), runif(1,min=0.1,max=2.5))   # c(shape, scale)
    sol <- try(solnp(para.start, Qaux.gamma, x=x, wt=w[,i],
             LB=lower.b, UB=upper.b,
             control=list(trace=FALSE, inner.iter=1000, outer.iter=1000,
             tol=1e-13, delta=1e-13)), FALSE)
    if ( !inherits(sol, "try-error")) {
       Shape[[i]] = sol$par[1]
       Scale[[i]] = sol$par[2]
       i = i + 1
    } else {
      nfail = nfail + 1
      if ( nfail >= 5 ) {
        stop("Too many problems with solnp!")
      }
    }
  }
  list(shape.l=Shape, scale.l=Scale)
}

#' M-step re-estimation function for multivariate gamma emission distribution
#'
#' @param x numeric vector
#' @param w weight matrix, size T x K
#' @param model HMM model as defined in hmmspec.gen(..)
#' @return emission slot of Hext object
#' @description Re-estimates the parameters of multivariate gamma emission distribution as part of the EM algorithm for HMM. No regularizer is used.
#' @details Needs the solnp function in the Rsolnp package.
#' @author Horace W Tso
#' @export
mstep.mvgamma = function(x, w) {
  require(Rsolnp)
  p = NCOL(x)
  K = NCOL(w)
  ix = apply(is.na(x), 1, any)
  if ( length(ix) > 0 ) {
    x = x[!ix, ,drop = FALSE]
    w = w[!ix, ,drop = FALSE]
  }
  emission = list(shape.l=list(), scale.l=list())
  lower.b = c(0.01, 0.01)
  upper.b = c(Inf, Inf)
  a.mat = matrix(NA, ncol=K, nrow=p)
  s.mat = matrix(NA, ncol=K, nrow=p)
  for ( j in 1:p ) {
    i = 1
    nfail = 0
    while ( i <= K ) {
      para.start = c(runif(1,min=0.5,max=10), runif(1,min=0.1,max=2.5))
      sol <- try(solnp(para.start, Qaux.gamma, x=x[,j], wt=w[,i],
               LB=lower.b, UB=upper.b,
               control=list(trace=FALSE, inner.iter=1000, outer.iter=1000,
               tol=1e-13, delta=1e-13)), FALSE)
      if (!inherits(sol,"try-error")) {
        a.mat[j,i] = sol$par[1]
        s.mat[j,i] = sol$par[2]
        i = i + 1
      } else {
        nfail = nfail + 1
        if ( nfail >= 5 ) {
          stop("Too many problems with solnp!")
        }
      }
    }
  }
  for ( i in 1:K ) {
    # columns are states
    emission$shape.l[[i]] = a.mat[,i]
    emission$scale.l[[i]] = s.mat[,i]
  }
  return(emission)
}


