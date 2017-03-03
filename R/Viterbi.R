#' Viterbi decoding algorithm
#'
#' @param model trained HMM model, output of hextfit
#' @param newdata data to decode with the trained model
#' @return list with elements that describe the trained model
#' @details the decoded sequence, likelihood, see code for details
#' @export
Viterbi.hext = function (model, newdata) {
  if (missing(newdata)) 
    stop("no data passed to predict!")
  else 
    x = newdata
  
  if (class(x) == "numeric" | class(x) == "integer") {
    warning("x is a primitive vector.  Assuming single sequence.")
    Nseq = 1
    TT = NROW(x)
    Tcum = c(0, TT)
    x = list(x=x, TT=TT)
  } else {
    Nseq = length(x$TT)
    TT = x$TT
    Tcum = cumsum(c(0, TT))
  }
  
  K = model$K
  lb = sapply(1:K, fn <- function(state) model$dens.emission(x$x, state, model))
  
  if ( model$Log == FALSE ) {    
    p[p == 0] = 1e-300
    tmp = model$ltransition
    tmp[!tmp > 0] = 1e-300
    logtrans = as.double(log(t(tmp)))
    tmp1 = model$linit
    tmp1[!tmp1 > 0] = 1e-300
    logpi = as.double(log(t(tmp1)))
  } else {
    logtrans = model$ltransition
    logpi = model$linit
  }
  
  state.seq = integer(sum(x$TT))
  loglik = 0
  for ( i in 1:Nseq ) {
    emit.prob = t(lb[(Tcum[i] + 1):Tcum[i + 1],])
    res <- .C("viterbi_hmm", 
              la = logtrans, 
              lpi = logpi, 
              lb = as.double(emit.prob), 
              TT = as.integer(x$TT[i]), 
              Nseq = as.integer(c(1)), 
              K = as.integer(K), 
              ss = as.integer(rep(-1,x$TT[i])), 
              loglik = as.double(c(0)), 
              PACKAGE = "Hext")
    loglik = loglik + res$loglik
    state.seq[(Tcum[i] + 1):Tcum[i + 1]] = res$ss + 1
  }
  ans <- list(s=state.seq, x=x$x, TT=x$TT, loglik=loglik)
  class(ans) <- "Hext.data"
  ans
}



