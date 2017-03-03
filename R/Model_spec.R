#' Model specificaiton in one object for hextfit call
#'
#' @param linit initial probability vector of hidden states, in log domain
#' @param ltrans transition probability matrix, in log domain 
#' @param parms.emission paramters for emission distribution
#' @param dens.emission density function for emission distribution
#' @param rand.emission function used to generate obs from emission distribution
#' @param mstep density function for the M-step
#' @param mtype integer code to indicate the type of emission distribution
#' @param K number of hidden states
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param Log boolean, whether log is used
#' @param Intv interal for guassian distribution
#' @param wdf degree of freedom for Wishart distribution  
#' @return list with all these elements for hexfit
#' @details see code for details
#' @export
hmmspec.gen = function (linit, ltrans, parms.emission, 
                        dens.emission=generic.density, 
                        rand.emission=NULL,
                        mstep=mstep.generic, mtype=mtype, 
                        K, M, D, Log=TRUE, Intv=FALSE, wdf=1) {
  ans <- list(K=K, M=M, D=D,
              linit=linit, ltransition=ltrans,
              parms.emission=parms.emission, 
              dens.emission=dens.emission,
              rand.emission=rand.emission, 
              mstep=mstep, 
              mtype=mtype, 
              Log=Log, 
              Intv=Intv, 
              wdf=wdf)
  class(ans) <- "Hext"
  return(ans)
}


#' Initialization of the HMM parameters
#'
#' @param p the entries of this vector give the dimension of the ith component of the obs
#' @param K number of hidden states  
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param M number of gaussian mixtures
#' @param A transition probability matrix, if NULL, the init code is used to generate one
#' @param Pi initial state probabiilty, if NULL, the init code is used to generate one
#' @param Log boolean, whether log is used
#' @return list consists of the initial state probability Pi, the transition matrix A, and emission probability b
#' @details see code for details
#' @export
gen.init = function(p, K, D=0, M=0, A=NULL, Pi=NULL, Log=FALSE) {
  require(MCMCpack)
  # p is a vector of length 5. 
  # K is a numeric value
  # D is a vector for multivariate discrete obs, otherwise just a value
  # M is a value
  # A, Pi are a list
  
  #  p = the entries of this vector give the dimension of the ith component of the obs
  #  K = no. of hidden states 
  #  D = a vector of the size of the alphabets in each column of the multivariate discrete obs matrix;
  #      one can think of it as number of multinomial outcomes for each component of the multivariate 
  #      discrete observation vector. D=2 means binary.
  #  M = no. of gaussian mixtures
  #  TT = no. of observation, or length of (single) obs sequence
  
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  # M1 : initialize prob table -------------------------------------------------
  # p[1] : dimension of the binary observations
  # if p[1] is non-zero, D must be non-zero
  if ( p[1] > 0 ) {
    for ( i in 1:K ) {
      lookup.mat = matrix(0, nrow=max(D), ncol=p[1])
      # dim : max(D) x N, where N is the dimension of the obs vector
      rownames(lookup.mat) = paste("f=",1:max(D),sep="")
      for ( j in 1:p[1] ) {
        dir.parm = rep(1, D[j])
        trial = as.vector(rdirichlet(1, dir.parm))
        while( any(trial < (0.05/D[j]) ) ) {
          trial = as.vector(rdirichlet(1, dir.parm))
        }
        lookup.mat[1:D[j], j] = trial
      }
      # as a check, the columns of this matrix should sum to 1.
      b$M1$prob.l[[i]] = lookup.mat
    }
  }
  
  # M2 : multivariate normal ---------------------------------------------------
  # p[2] : dimension of the univariate observations
  if ( p[2] > 0 ) {
    V = crossprod(matrix(rnorm(p[2]*p[2], sd=3), p[2], p[2]))
    for ( i in 1:K ) {
      b$M2$mu.l[[i]] = runif(p[2], min=-2, max=2)
      b$M2$sigma.l[[i]] = rwish(p[2]+i,V)
    }
  }
  
  # M3 : gamma -----------------------------------------------------------------
  if ( p[3] > 0 ) {
    for ( i in 1:K ) {
      b$M3$shape.l[[i]] = runif(p[3], min=1, max=10)
      b$M3$scale.l[[i]] = runif(p[3], min=0.1, max=3)
    }
  }
  
  # M4 : beta ----------------------------------------
  if ( p[4] > 0 ) {
    for ( i in 1:K ) {
      b$M4$shape1.l[[i]] = runif(p[4], min=0.5, max=7)
      b$M4$shape2.l[[i]] = runif(p[4], min=0.5, max=7)
    }
  }
  
  # M5 : gaussian mixture --------------------------------------------------------
  # p[5] : dimension of the multivariate observations
  # if p[5] is non-zero, M must be greater than 1.
  # There is a specialized function to initialize gaussian mixture starting values
  # ==============================================================================
  if ( p[5] > 0 ) {
    V5 = crossprod(matrix(rnorm(p[5]*p[5], sd=2), p[5], p[5]))
    # cat("M:", M, "\n")
    dir.parm = rep(1, M)
    for ( i in 1:K ) {
      b$M5$c.l[[i]] = as.vector(rdirichlet(1, dir.parm))
      for ( m in 1:M ) {
        b$M5$mu.ll[[i]][[m]] = rnorm(p[5])
        b$M5$sigma.ll[[i]][[m]] = rwish(p[5]+1,V5)  
      }
    }
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha =============
  if ( is.null(Pi) ) {
    dir.parm = rep(1, K)
    if ( Log) { 
      Pi = log(as.vector(rdirichlet(1, dir.parm)))
    } else {
      Pi = as.vector(rdirichlet(1, dir.parm))
    }
  } 
  # initial transition prob generated by Dirichlet distribution with equal alpha ==============
  if ( is.null(A) ) {
    dir.parm = rep(1, K)
    if ( Log ) {
      A = log(rdirichlet(K, dir.parm))
    } else {
      A = rdirichlet(K, dir.parm)
    }
  }
  list(Pi=Pi, A=A, b=b)
}



