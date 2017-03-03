#' Simulate gaussian mixture data. Internal function called by simulate.gaussian.mix
#'
#' @param N number of data points to generate
#' @param mu a list of mean vectors, length equal to the number of mixtures
#' @param sigma a list of covariance matrices, length = no. of mixtures 
#' @param cj a vector with the fraction to draw from each component, must sum to 1.0
#' @description See example
#' @return matrix of size N x p, where p is the dimension of the data
#' @export
sim.gaussmix = function(N, mu, sigma, cj) {
  require(mvtnorm)
  M = length(cj)
  p = ncol(sigma[[1]])
  dat = rep(list(list()), M)
  for ( i in 1:M ) {
    # rmvnorm returns a N x p matrix
    dat[[i]] = rmvnorm(N, mean=mu[[i]], sigma=sigma[[i]])
  }
  # ==============================
  draws = rmultinom(N, 1, prob=cj)
  # draws is a M x N matrix of zeros and ones, where each 
  # column has only one 1, where the number of 1s in each
  # row is roughly the proportion given in cj.
  # ======================================================
  i.pick = apply(draws, 2, function(xx) which(xx==1))
  # i.picks is a vector consisting of integers from 1 to M.
  # this is the random picks according to probability vector cj 
  # from the M components
  simdat = matrix(NA, ncol=p, nrow=N)
  for ( i in 1:N ) {
    tmp = dat[[i.pick[i]]]
    # tmp is a N x p matrix
    simdat[i,] = tmp[i,]
  }
  return(simdat)
}

#' Simulate gaussian mixture data (wrapper function)
#'
#' @param N number of data points to generate
#' @param mix.st mixture structure M5 as defined in gen.init of the Hext package
#' @param state integer indicating one of the K states to use for this simulation
#' @description wrapper function that calls sim.gaussmix(..)
#' @return matrix of size N x p, where p is the dimension of the data
#' @export
simulate.gaussian.mix = function(N, mix.st, state ) {
  mu = mix.st$mu.ll[[state]]
  sigma = mix.st$sigma.ll[[state]]
  cj = mix.st$c.l[[state]]
  dat = sim.gaussmix(N, mu, sigma, cj)
  return(dat)
}

#' Simulate gaussian mixture data, returns a matrix (wrapper function)
#'
#' @param Len vector of lengths of the sequences to generate
#' @param mix.st mixture structure M5 as defined in gen.init of the Hext package
#' @param SS vector of the length of each state, e.g. {100, 150, 130, 200, 180}
#' @param st.progression integer vector for the state progression, e.g. {1, 3, 2, 1, 4} 
#' @description wrapper function that calls simulate.gaussian.mix(..)
#' @details Length of TT must equal the length of st.progression.
#' @return matrix of size N x p, where p is the dimension of the data
#' @export
simulate.gauss.mix.matrix = function(Len, mix.st, st.progression) {
  dat = NULL
  n = length(Len)
  for ( i in 1:n ) {
    tmp = simulate.gaussian.mix(Len[i], mix.st, st.progression[i])
    dat = rbind(dat, tmp)
  }
  return(dat)
}

#' Generate multiple sequences of multivariate gaussian mixture data 
#'
#' @param TT vector of the lengths of each sequence
#' @param LL vector of the lengths of each state in a sequence
#' @param mix.st mixture structure with mu, sigma, cj for each state and each mixture component
#' @param st.progression vector of integer indicating the state progression, e.g. {1, 3, 2, 1, 4}
#' @description generate multiple sequences of multivariate gaussian mixture data. Calls simulate.gauss.mix.matrix()
#' @return list with two elements : dat, the data matrix; ss : the state sequences
#' @export
simulate.GMM.Nseq = function(TT, LL, mix.st, st.progression) {
  if ( length(LL) != length(st.progression)) {
    error("length of SS not equal to st.progression.")
  }
  Nseq = length(TT)
  n = sum(LL)
  frac = LL/n
  dat = NULL
  ss = NULL
  for ( i in 1:Nseq) {
    len = ceiling(TT[i]*frac)
    tmp = rep(st.progression, len)
    ss = c(ss, tmp[1:TT[i]])
    tmp1 = simulate.gauss.mix.matrix(len, mix.st, st.progression)
    dat = rbind(dat, tmp1[1:TT[i],])
  }
  list(dat=dat, ss=ss)
}

#' Generate multinomial sequence data
#'
#' @param TT vector of the lengths of each state
#' @param D 
#' @param st.progression vector of integer indicating the state progression, e.g. {1, 3, 2, 1, 4}
#' @description Note : TT and st.progress should have the same length
#' @return matrix of size n x p
#' @export
simulate.multinom = function(TT, D, st.progress) {
  K = max(st.progress)
  n = length(TT)
  pr = matrix(0, ncol=D, nrow=K)
  for ( k in 1:K ) {
    pr[k,] = as.vector(rdirichlet(1, rep(1,D)))
  }
  res = NULL
  for ( i in 1:n ) {
    tmp = rmultinom(TT[i], size=1, prob=pr[st.progress[i],])
    vec = apply(tmp, 2, function(d)which(d==1))
    res = c(res, vec)
  }
  return(res)
}


#' Generate multivariate binary data
#'
#' @param TT vector of sequence lengths 
#' @param p dimension of observation
#' @param K number of states
#' @description generate multivariate binary sequence
#' @return matrix of size n x p
#' @export
simulate.mv.binary = function(TT, p, K) {
  Nseq = length(TT)
  Ttot = sum(TT)
  theta = runif(K, min=0.05, max=0.95)
  
  mat = NULL
  for ( i in 1:Nseq ) {
    n = ceiling(TT[i] / K)
    for ( k in 1:K ) {
      tmp = matrix(rbinom(n*p, 1, prob=theta[k]), ncol=p)
      mat = rbind(mat, tmp)
    }
  }
  mat = mat[1:Ttot,]
  # following four lines convert the 0/1 matrix into 1/2 matrix 
  if ( p == 1) {
    tmp1 = as.character(mat)
  } else {
    tmp1 = apply(mat, 2, as.character)  
  }
  tmp.f = sapply(tmp1, factor)  # this collapses the matrix into a vector of factors
  tmp.lvl = as.integer(tmp.f)  # remains a vector
  x.lvl = matrix(tmp.lvl, ncol=p)  # turns it back to a matrix for use below
  return(x.lvl)
}

