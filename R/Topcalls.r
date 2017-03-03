#' Train HMM with model specificaiton provided in `model`
#'
#' @param x vector or matrix of data
#' @param model model specification created by hmmspec.gen(..)
#' @param lock.transition boolean, whether to fix the transition matrix
#' @param tol numeric tolerance for convergence evaluation
#' @param maxit maximum number of EM iterations
#' @param verbose boolean, whether to print diagnostics
#' @return list with elements that describe the trained model
#' @details See example
#' @export
hextfit = function(x, model, lock.transition=FALSE, tol=1e-01, maxit=1000, verbose=TRUE) {
  K = model$K
  D = model$D  # size of the alphabet sets, this is a vector!
  M = model$M  # no. of gaussian mixtures
  if ( class(x)=="numeric" | class(x)=="integer") {
     warning('x is a primitive vector. Assuming single sequence.')
     TT = NROW(x)
     Nseq = 1
  } else {
    TT = x$TT
    Nseq = length(TT)
    x = x$x
  }
  mtype = model$mtype
  if ( length(mtype) != NCOL(x) ) stop("length of mtype differs from no of columns in x.") 
  if ( K < 2 ) stop("K must be larger than one.")
  if( any(dim(model$ltransition)!=K) )  stop("dimensions of ltransition incorrect!")
  if( length(model$linit)!=K ) stop("dimensions of linit incorrect!")
  if( NROW(x)!=sum(TT) )  stop("no. of rows in x not equal sum(TT) !")

  mstep = model$mstep
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  Emit.Prob.Func = model$dens.emission
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  loglik1 = numeric(maxit)
  loglik1[] = NA
  loglik1[1] = -1000
  loglik2 = numeric(maxit)
  loglik2[] = NA
  loglik2[1] = -1000

  # initialize gamma ===========
  lgam = rep(-1.0e100, K*sum(TT))
  # ============================
  lAlpha = rep(-1.0e300, (K+1)*sum(TT))
  lBeta = rep(-1.0e300, K*sum(TT))

  ll = c(0.0, 0.0)
  mu.it = rep(list(list()), maxit)
  
  for(i in 1:maxit) {
    if (verbose == TRUE) cat("\niter:",i," ====")
    # calculate emission prob by state given model >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    emit.p = sapply(1:K, fn <- function(state) Emit.Prob.Func(x, state, model))
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # E-step <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    result = .C("mo_estep_hmm_ln",
        lA=as.double(t(model$ltransition)),
        lPi=as.double(t(model$linit)),
        lP=as.double(t(emit.p)),
        TT=as.integer(TT),
        nsequences=as.integer(Nseq),
        K=as.integer(K),
        forward_p=lAlpha,
        backward_p=lBeta,
        lgam=lgam,
        ll=as.double(ll),
        PACKAGE='Hext')

    loglik1[i] = result$ll[1]
    loglik2[i] = result$ll[2]
    if (verbose == TRUE) { 
      cat("loglikelihood(alpha) ", result$ll[1])
      cat(".....(beta)", result$ll[2])
    }
    if ( i > 1 ) {
      if(abs(loglik1[i]-loglik1[i-1]) < tol) {
        if (verbose == TRUE) cat("\nConverged.\n")
        break
      }
    }
    # M-step >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gam = matrix(exp(result$lgam), ncol=K)
    # simple check on gamma
    gamsum = rowSums(gam)
    
    model$parms.emission = mstep(x, gam, model)

    # write the optimized parameters back to model
    if(!lock.transition) {
      model$ltransition = matrix(result$lA, nrow=K, byrow=TRUE)
      model$linit = result$lPi
    }
    # save the means vector for output
    mu.it[[i]] = model$parms.emission$M5$mu.ll
  }
  
  LL = loglik1[i]    # last iteration
  # Calculate no of parameters in the model >>>>>>>>>>>>>>>>>>>>>>>>>
  # 1 = discrete
  ix1 = which(mtype == 1)
  if ( length(ix1) > 0 ) {
    n.par1 = K*length(ix1)*(D - 1)
  } else {
    n.par1 = 0
  }
  # 2 = normal/multivariate
  ix2 = which(mtype == 2)
  if ( length(ix2) > 0 ) {  
    # Zucchini & MacDonald (2009), p.186  
    # K*p(p+3)/2, where p is the dimension of the observation, and is equal
    # the number of 2s in mtype.
    n.par2 = K*length(ix2)*(length(ix2) + 3)/2
  } else {
    n.par2 = 0
  }
  # 3 = gamma
  ix3 = which(mtype == 3)
  if ( length(ix3) > 0 ) {  
    # two free parameters for each component obeying the gamma distribution
    n.par3 = K * length(ix3) * 2
  } else {
    n.par3 = 0
  }
  # 4 = beta
  ix4 = which(mtype == 4)
  if ( length(ix4) > 0 ) {  
    # two free parameters for each component obeying the beta distribution     
    n.par4 = K * length(ix4) * 2
  } else {
    n.par4 = 0
  }
  # 5 = gaussian mixture
  ix5 = which(mtype == 5)
  if ( length(ix5) > 0 ) {  
    # for each multivariate component, there are M mixtures, and each mixture 
    # has p(p+3)/2 parameters, thus total :
    n.par5 = K * length(ix5)* M * (length(ix5) + 3) / 2
  } else {
    n.par5 = 0
  }
  np.trans = K * ( K - 1 )
  # Total number of parameters in the HMM
  N.param = np.trans + n.par1 + n.par2 + n.par3 + n.par4 + n.par5
  # Akaike information criterion (Ref : wiki on AIC)
  Aic = - 2*LL + 2*N.param 
  # AIC with finite sample correction (see wiki on AIC)
  # Burnham, Anderson (2002) Model selection and multimodel inference: A Practical
  #   Information-Theoretic Approach, 2nd ed, Springer.
  Aicc = Aic + 2*N.param*(N.param + 1)/(sum(TT) - N.param - 1)
  Bic = - 2*LL + N.param*log(sum(TT))
  
  ret = list(model=model,mstep=mstep, gam=gam, lgam=result$lgam, last.it=i,
         loglik1=loglik1, loglik2=loglik2, TT=TT, yhat=apply(gam,1,which.max),
         Aic=Aic, Aicc=Aicc, Bic=Bic, mu.it=mu.it)
  class(ret) <- "Hext"
  return(ret)
}

