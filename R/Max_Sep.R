


#' Max Sep algorithm that find M points in the data cloud that are maximally separated
#'
#' @param dat data matrix
#' @param M number of maximally separate points
#' @description see my paper for detailed description
#' @return matrix
#' @export
max.separate.N = function(dat, M) {
  require(fields)
  N = NROW(dat)
  Dmat = rdist(dat)
  first.pair = which(Dmat == max(Dmat), arr.ind=TRUE)
  Smax = c(first.pair[1,1], first.pair[2,1])
  for ( m in 1:(M-2) ) {
    sub.i = (1:N)[-Smax]
    Fi = rep(1e300, N)
    for ( i in 1:length(sub.i) ) {
      q.i = rep(1, length(Smax))
      Fi[sub.i[i]] = F.potential(dat[Smax,], dat[sub.i[i],], q.i)
    }
    # find the point with the smallest field potential
    ix = which.min(Fi)
    Smax = c(Smax, ix)
  }
  return(Smax)
}

#' Field potential
#'
#' @param src sources
#' @param pt point vector
#' @param q.i electric charges for each source point
#' @description see my paper for detailed description
#' @return potential matrix
#' @export
F.potential = function(src, pt, q.i) {
  N = NROW(src)
  p = NCOL(src)
  r = rep(0, N)
  pot = 0
  for ( i in 1:N ) {
    r[i] = sqrt(sum((src[i,] - pt)^2) )
    pot = pot + q.i[i]/r[i]
  }
  return(pot)
}
