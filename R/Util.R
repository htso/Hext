#' Computes the euclidean distances for n pairs of points
#'
#' @param tgt target, a vector of length p 
#' @param pts n x p matrix, where each row is a point in p-dimensional space 
#' @description Computes the euclidean distances for n pairs of points, e.g. pt0 vs pts[1,], pts[2,], pts[3,], etc. 
#'              It could be slow for large number of points.
#' @return vector of distances 
#' @export
dist.one = function(tgt, pts) {
  del2 = (t(pts) - tgt)^2
  d = sqrt(apply(del2, 2, sum))
  # R runs by columns, so when a matrix is subtracted by a vector
  # it takes the first column of the matrix, subtract the vector, then 2nd column, so on.
  #twopts.dist = function(x) {sqrt(sum((x - pt0)^2))}
  #d = apply(pts, 1, twopts.dist)
  return(d)
}

#' Computes the euclidean distances for n pairs of points using faster C codes
#'
#' @param tgt target, a vector of length p 
#' @param pts n x p matrix, where each row is a point in p-dimensional space 
#' @description Computes the euclidean distances for *all* pairs of points, e.g. pt0 vs pts[1,], pts[2,], pts[3,], etc.
#'              Similar to dist.one(). This version uses the build-in dist() function
#' @return matrix of distance
#' @export
dist.one0 = function(tgt, pts) {
  mat = rbind(tgt, pts)
  d = as.matrix(dist(mat, method="euclidean"))
  return(d[,1][-1])
}

#' Compute the density of points 
#'
#' @param dat numeric matrix of size n x p, where n is the nubmer of points in p-dimensional space
#' @description compute the density in an ellipse
#' @return numeric value which is the density in p-dim ellipsoid
#' @export
density.ndim = function(dat) {
  p = NCOL(dat)
  range.vec = apply(dat, 2, function(x){max(x)-min(x)})
  # formula for n-dimensional ellipsoid volume from wiki (search keyword : ellipsoid)
  radii = range.vec / 2  # vector of length p
  ave.radii = mean(radii)
  C.n = pi^(p/2) / gamma(p/2 + 1)
  ellipsoid.volume = prod(radii) * C.n  
  pt.density = n / ellipsoid.volume
  return(pt.density)
}

#' Compute the expected number of neighbors in a ball of given radius
#'
#' @param dat numeric matrix of size n x p, where n is the nubmer of points in p dimensional space
#' @param ball.r radius of a p-dimensional sphere
#' @description compute the expected number of neighbors
#' @return number of points in the sphere 
#' @export
expected.neighbors = function(dat, ball.r) {
  p = NCOL(dat)
  n = NROW(dat)
  den = density.ndim(dat)
  C.n = pi^(p/2) / gamma(p/2 + 1)
  ball.volume = C.n * (ball.r)^p 
  expected.neighbors = ball.volume * den
  return(expected.neighbors)
}


#' Check each covariance matrix in the gaussian mixture specification
#'
#' @param model Hext model structure with an non-empty M5 slot
#' @description Display the sum of log eigenvalues of each covariance matrix
#' @return None
#' @export
cov.mat.check = function(model) {
  M = model$M
  K = model$K
  cat("----- determinant of cov matrices ------\n")
  for ( i in 1:K ) {
    cat("i:", i, "\t")
    for ( k in 1:M ) {
      cat("k:", k, "\t")
      sigma = model$parms.emission$M5$sigma.ll[[i]][[k]]
      ev = eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
      logdet = sum(log(ev))
      cat(exp(logdet), "\t")
    }
    cat("\n")
  }
  cat("----------------------------------------")
}

#' Check each covariance matrix in the gaussian mixture specification
#'
#' @param K number of viterbi states 
#' @param ss.vitb vector of the viterbi-decoded states in the same time order as the observation sequence
#' @param dat data frame having forward return columns, whose names start with "F"
#' @description See example
#' @return list of Ret.by.state, Sd.by.state, event.counts, state.counts, enum
#' @export
Transition.events.by.return = function(K, ss.vitb, dat) {
  # enumerate all possible transitions
  # note that the matrix is not symmetric, i.e i-> j is not the same as j-> i
  enum.tran = expand.grid(x=1:K, y=1:K)
  # eliminate self transitions 
  enum.tran = enum.tran[-which(enum.tran[,1]==enum.tran[,2]),]
  colnames(enum.tran) = c("from", "to")
  rownames(enum.tran) = 1:nrow(enum.tran)
  # find all transition events
  # ss1 is one period lag of ss
  ss = ss.vitb
  ss1 = c(NA, ss[-length(ss)])
  # the event vector keeps track of transition events
  event = rep(NA, length(ss))
  # a non-NA entry in the vector tells you that a transition has occurred at that time point
  # ie. it went from state i yesterday to state j today
  for ( i in 2:length(ss) ) {
    for ( j in 1:nrow(enum.tran)) {
      # yesterday was in state i and today state j, then an event is recorded
      if ( ss1[i]==enum.tran[j,1] && ss[i]==enum.tran[j,2])
        event[i] = j
    }
  }
  dat = cbind(dat, event)
  ic = grep(pat="F", colnames(dat), ignore.case=FALSE)
  Ret.by.state = aggregate(dat[,ic], list(event=event), mean)
  Sd.by.state = aggregate(dat[,ic], list(event=event), sd)
  event.counts = table(dat[,"event"])
  state.counts = table(ss)
  list(Ret.by.state=Ret.by.state, Sd.by.state=Sd.by.state, event.counts=event.counts, state.counts=state.counts, enum=enum.tran)
}

#' Internal function for testing the sensitivity of HMM parameters
#'
#' @param dat : matrix with rows as observation and columns as features
#' @param TT : vector of sequence lengths, one entry per sequence
#' @param min.length : the minimum length of the *last* sequence
#' @param shift : the number of observations in the first sequence to discard so that the whole continuous observation is shifted forward
#' @description Purpose of this function is to test how sensitive an hmm model is to 
#' different beginning of the sequence division. Its use is very limited
#' to the approach taken in DiscretRet.R. Some analysis needs to be made
#' when modeling long sequences that evaluates parameter stability when a 
#' few data points is added or taken away at the beginning and end of 
#' the sequence. 
#' @return list of x, N
#' @export
shift.sequence = function(dat, TT, min.length=10, shift=5) {
  n = length(TT)
  len = nrow(dat)
  tmp = dat
  if ( TT[n] > min.length + shift ) {
    tmp = tmp[-(1:shift),]
    TT[n] = TT[n] - shift  
  } else if ( TT[n] > shift ) {
    tmp = tmp[(shift+1):(len - TT[n] + shift),]
    TT = TT[-n]
  } else {
    error("length requirement not met.")
  }
  return(list(x=tmp, N=TT))
}

#' Inspect vector or matrix for number of NA, NaNs, and infinities
#'
#' @param x vector or matrix
#' @description Count the number of NAs, NaNs, and infinities in a matrix or vector
#' @return list of ina, number of NAs, inan, number of NaN, and iinf, number of infinities
#' @export
weirdness = function(x) {
  if (NCOL(x) == 1 ) {
    ina = sum(is.na(x))
    inan = sum(is.nan(x))
    iinf = sum(is.infinite(x))
    return(list(ina=ina, inan=inan, iinf=iinf))
  } else if ( NCOL(x) > 1 ) {
    ina = apply(x, 2, function(y) sum(is.na(y)))
    inan = apply(x, 2, function(y) sum(is.nan(y)))
    iinf = apply(x, 2, function(y) sum(is.infinite(y)))
    return(list(ina=ina, inan=inan, iinf=iinf))
  } else {
    return(list(ina=NA, inan=NA, iinf=NA))
  }
}

#' Compute log(X+Y) given log(X) and log(Y)
#'
#' @param lX log(X), a scalar
#' @param lY log(Y), a scalar
#' @description Calculate log(X+Y), Scalar version
#' @return numeric value
#' @export
Elnsum = function(lX, lY) {	
  myeLnX = max(lX, lY)
	myeLnY = min(lX, lY)
  # Based on R 2.15.0 under 64-bit Windows 7,
  # log(-745) = 0
	if ( myeLnY <= -745 ) {
		return( myeLnX )
	} else {
		return(myeLnX + Eln(1.0 + exp(myeLnY - myeLnX)))
	}
}

#' Compute log(X+Y) given log(X) and log(Y), vectorized version
#'
#' @param lx vector of log(x1), log(x2), .... log(xn)
#' @description Calculate log(x1 + x2 + ... = xn), vectorized version
#' @return numeric value
#' @export
Elnsum.vec = function(lx) {
  cum = -Inf
  n = length(lx)
  for ( i in 1:n ) {
    cum = Elnsum(cum, lx[i])
  }
  return(cum)
}

#' Compute the logarithm of a value
#' 
#' @param x numeric value
#' @description Internal function to compute the log of a value
#' @return log(x) if x > 0, otherwise it gives -Inf
#' @export
Eln = function(x) {
  if (x > 0.0) {
    return(log(x))
  } else {
    return(-Inf)
  }
}


