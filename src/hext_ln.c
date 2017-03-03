#include "stdio.h"
#include <R.h>
#include "stdlib.h"

#define int8 unsigned char
//#define LOGZERO -DBL_MAX
#define LOGZERO  (-1.0E15)   /* ~log(0) */
#define LSMALL (-1.0E15)   /* log values < LSMALL are set to LZERO */
#define MINEARG (-708.3)   /* lowest exp() arg  = log(MINLARG) */
#define MINLARG 2.45E-308  /* lowest log() arg  = exp(MINEARG) */

#ifndef MIN
	#define MIN(p,q) ((p) < (q) ? (p) : (q))
#endif //MIN
#ifndef MAX
	#define MAX(p,q) ((p) > (q) ? (p) : (q))
#endif //MAX

/* ==========================================================================*/
double eexp(const double theX)
{
	if (theX <= LOGZERO)
		return(0.0L) ;
	else
		return(exp(theX)) ;
}

double eln(const double theX)
{
	if (theX > 0.0L)
		return(log(theX)) ;
	else
		return(LOGZERO) ;
}

double elnsum1(const double theX, const double theY)
{
double	myeLnX = eln(theX),
		myeLnY = eln(theY) ;

	if ( (myeLnX <= LOGZERO) || (myeLnY <= LOGZERO) )
	{	if (myeLnX <= LOGZERO)
			return(myeLnY) ;
		else
			return(myeLnX) ;
	} else{
	    if (myeLnX > myeLnY)
			return(myeLnX + eln(1.0L+exp(myeLnY-myeLnX))) ;
		else
			return(myeLnY + eln(1.0L+exp(myeLnX-myeLnY))) ;
	}
}

double elnsum(const double theeLnX, const double theeLnY)
{
// elnsum(eln(x), eln(y)) = eln(x+y) pour x, y > LOGZERO
// elnsum(LOGZERO, eln(y)) = eln(y)
// elnsum(eln(x), LOGZERO) = eln(x)
double	myeLnX = MAX(theeLnX, theeLnY),
		myeLnY = MIN(theeLnX, theeLnY) ;

	if (myeLnY <= LOGZERO)
		return(myeLnX) ;
	else
		return(myeLnX + eln(1.0L+exp(myeLnY-myeLnX))) ;
}

double elnproduct1(const double theX, const double theY)
{
double	myeLnX = eln(theX),
		myeLnY = eln(theY) ;

	if ( (myeLnX <= LOGZERO) || (myeLnY <= LOGZERO) )
		return(LOGZERO) ;
	else
		return(myeLnX + myeLnY) ;
}

double elnproduct(const double theeLnX, const double theeLnY)
// elnproduct(eln(x), eln(y)) = eln(x) + eln(y) pour x, y > 0
// elnproduct(LOGZERO, eln(y)) = elnproduct(eln(x), LOGZERO) = LOGZERO
{
	if ( (theeLnX <= LOGZERO) || (theeLnY <= LOGZERO) )
		return(LOGZERO) ;
	else
		return(theeLnX + theeLnY) ;
}
// =============================================================================

void checkmem(void *x) {
  if(x==NULL) error("Out of memory.");
}

int min(int a,int b) {
  if(a<b) return(a);
  else return(b);
}

void **alloc_matrix(int nrow,int ncol,int size) {
  int i;
  void **x = malloc(sizeof(void *)*nrow);
  checkmem(x);
  for(i=0;i<nrow;i++)	{
    x[i]=malloc(size*ncol);
    checkmem(x[i]);
  }
  return(x);
}

void free_matrix(int nrow,int ncol,void **x) {
  int i;
  for(i=0;i<nrow;i++)	free(x[i]);
  free(x);
}

void print_matrix(int nrow,int ncol,double *x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)
      Rprintf("%.3g\t",x[i*ncol+j]);
    //			Rprintf("%.3f\t",x[j*nrow+i]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_matrix2(int nrow,int ncol,double **x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)
      Rprintf("%.3g\t",x[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_imatrix2(int nrow,int ncol,int **x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)
      Rprintf("%d\t\t",x[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

// ================================================================================================
void forward_ln(double *la, double *lpi, double *lb, int offset, int *timelength,
     int *nstates, double ***output) {
	 
  int K = *nstates;
  int T = *timelength;
  int i,k,t;
  double **lalpha = *output;
  double cum = LOGZERO;

  for(i=0;i<K;i++) {
    lalpha[i][0] = lpi[i] + lb[offset*K+i];
  }

  for(t=1; t<T; t++) {
    for(k=0; k<K; k++) {
      cum = LOGZERO;
      for(i=0; i<K; i++) {
        cum = elnsum(cum, lalpha[i][t-1] + la[i*K+k]);
      }
      lalpha[k][t] = cum + lb[offset*K+t*K+k];
      //Rprintf("lalpha[k:%d,t:%d] = %.5g\n", k,t, lalpha[k][t]);
    }
  }
}

void backward_ln(double *la, double *lb, int offset, int *timelength, int *nstates, double ***output) {

  int K = *nstates;
  int T = *timelength;
  int i,j,t;
  double cum = LOGZERO;
  double **lbeta = *output;

  for(i=0;i<K;i++) {
    lbeta[i][T-1] = 0.0L;
  }
  //Notice that t recurses back from T-1 to 1.
  for(t=T-2; t>=0; t--) {
    for(i=0; i<K; i++) {
      cum = LOGZERO;
	  for(j=0; j<K; j++) {
          cum = elnsum(cum, la[i*K+j] + lb[offset*K+(t+1)*K+j] + lbeta[j][t+1]);
	  }
	  lbeta[i][t] = cum;
    }
  }
}

// ===== PRODUCTION VERSION ===========================================================================
void mo_estep_hmm_ln(double *la, double *lpi, double *lb, int *T, int *nsequences, int *nstates,
					   double *forward_p, double *backward_p, double *lgam, double *loglik) {
  int K = *nstates;
  int Nseq = *nsequences;
  double ***lalpha, ***lbeta, *lxi;
  double normalizer=LOGZERO, cum=LOGZERO, num=LOGZERO, denom=LOGZERO;
  double val=0.0L, asum=0.0L;
  int i, j, t, n;
  int *Tcum = (int *)malloc((Nseq+1)*sizeof(int));
  double *LLa = (double *)malloc(Nseq*sizeof(double));
  double *LLb = (double *)malloc(Nseq*sizeof(double));
  double tmp, tmp1, tmp2;

  Tcum[0]=0;
  for( n=1; n < Nseq+1; n++) {
    Tcum[n] = T[n-1] + Tcum[n-1];
  }
  // K x K x sum(T[1], T[2], ... T[N])
  lxi = (double *)malloc(K*K*Tcum[Nseq]*sizeof(double));

  // lalpha and lbeta have the same size in this log version of EM
  lalpha = (double ***)malloc(Nseq*sizeof(double **));
  lbeta  = (double ***)malloc(Nseq*sizeof(double **));

  for( n=0; n < Nseq; n++ ) {
    lalpha[n] = (double **)malloc(K*sizeof(double *));
    lbeta[n]  = (double **)malloc(K*sizeof(double *));

     // set pointers to the memory chunks pass down from R
    for(i=0; i<K; i++) {
        lalpha[n][i] = forward_p  + i*Tcum[Nseq] + Tcum[n];
        lbeta[n][i]  = backward_p + i*Tcum[Nseq] + Tcum[n];
    }
	 // calling Alpha, Beta to perform forward and backward recursions
    forward_ln( la, lpi, lb, Tcum[n], T+n, nstates, &lalpha[n]);
    backward_ln(la,      lb, Tcum[n], T+n, nstates, &lbeta[n]);

     // Gamma ---------------------------------------------------------------------------
    for( t=0; t < T[n]; t++) {
        normalizer = LOGZERO;
        for (i=0; i<K; i++) {
		    lgam[i*Tcum[Nseq] + Tcum[n] + t] = lalpha[n][i][t] + lbeta[n][i][t];
		    normalizer = elnsum(normalizer, lgam[i*Tcum[Nseq] + Tcum[n] + t]);
        }
	    for (i=0; i<K; i++) {
		    lgam[i*Tcum[Nseq] + Tcum[n] + t] = lgam[i*Tcum[Nseq] + Tcum[n] + t] - normalizer;
		      //Rprintf("lgam[n:%d,i:%d,t:%d] = %.5g\n",n,i,t,lgam[i*Tcum[Nseq] + Tcum[n] + t]);
        }
    }

	 // Xi ------------------------------------------------------------------------------------
	 // Mann's Algo 8, Rabiner, p.351-352 ---------------------------------------------------------------
  	for(t=0; t < T[n]-1; t++) {
		 normalizer = LOGZERO;
    	 for (i=0; i<K; i++) {
    		for (j=0; j<K; j++) {
              //a = la[i*K+j];
              lxi[Tcum[n]*K*K + t*K*K + i*K + j] = lalpha[n][i][t] + la[i*K+j] + lb[Tcum[n]*K+(t+1)*K+j] + lbeta[n][j][t+1];
  			  normalizer = elnsum(normalizer, lxi[Tcum[n]*K*K + t*K*K + i*K + j]);
    		}
    	}
    	for (i=0; i<K; i++) {
    		for (j=0; j<K; j++) {
    			lxi[Tcum[n]*K*K + t*K*K + i*K + j] = lxi[Tcum[n]*K*K + t*K*K + i*K + j] - normalizer;
    		}
    	} //i-for
	} // t-for
  } // n-for

  // Calculate loglikelihood using alpha[i][T]
  asum = 0.0L;
  for(n=0; n < Nseq; n++) {
    LLa[n] = LOGZERO; // need to initialize each of the cum[]
    for(i=0; i < K; i++) {
      LLa[n] = elnsum(LLa[n], lalpha[n][i][T[n]-1]);
      //Rprintf("ll cum[i:%d,n:%d] = %.5g\n",i,n,cum);
    } // i-for
    asum += LLa[n];
  } // n-for
  //Rprintf("ll(alpha) = %.5g\n", asum);
  loglik[0] = asum;

  // Calculate loglikelihood using beta[i][1]
  asum = 0.0L;
  for(n=0; n < Nseq; n++) {
    LLb[n] = LOGZERO;
    for(i=0; i < K; i++) {
      LLb[n] = elnsum(LLb[n], lbeta[n][i][0] + lpi[i] + lb[Tcum[n]*K+i]);
    } // i-for
    asum += LLb[n];
  } // n-for
  //Rprintf("ll(beta) = %.5g\n", asum);
  loglik[1] = asum;

  // Pi, or initial prob ----------------------------------------------------------
  // see p. 10, Bilmes, Algo 9 of Mann --------------------------------------------
  for(i=0; i<K; i++) {
      cum = LOGZERO;
	  for(n=0; n < Nseq; n++) {
	      // lgam evaluated at t=0
		  cum = elnsum(cum, lgam[i*Tcum[Nseq] + Tcum[n]]);
	  } // n-for
      lpi[i] = cum - eln(Nseq);
  } // i-for

  // a(i,j) -----------------------------------------------------------------------
  // Rabiner Eq (109), p.273 multiple sequence reestimation formulas --------------------------------------------
  // Unlike Rabiner, i calculate the denominator from gamma, which already has the
  // P(O|lambda) built in. Therefore, no need to have the 1/P_k factor.
  for(i=0; i < K; i++) {
    denom = LOGZERO;
    for(n=0; n < Nseq; n++) {
        for(t=0; t < T[n] - 1; t++) {
		  denom = elnsum(denom, lgam[i*Tcum[Nseq] + Tcum[n] + t]);
		  //Rprintf("-----> denom(t=%d, i=%d) : %.5g\n", t, i, denom);
		}
    }
	// Similary, unlike Rabiner, I calculate the numerator from xi, which
	// is normalized already, thus no need to include the 1/P_k factor in
	// the sum over sequence.
    for(j=0; j < K; j++) {
      num = LOGZERO;
      for(n=0; n < Nseq; n++) {
		for(t=0; t < T[n]-1; t++) {
           num = elnsum(num, lxi[Tcum[n]*K*K + t*K*K + i*K + j]);
           //Rprintf("num(t=%d, j=%d, i=%d) : %.5g\n", t, j, i, num);
		}
	  }
	  la[i*K+j] = num - denom;
	}
  }

  // release memory
  for(n=0; n<Nseq; n++){
    free(lbeta[n]);
    free(lalpha[n]);
  }
  free(lalpha);
  free(lbeta);
  free(lxi);
  free(Tcum);
}

void viterbi_hmm(double *la, double *lpi, double *lb, int *T, int *nsequences, int *nstates, int *q, double *loglik) {

  int K = *nstates;
  int Nseq = *nsequences;
  int i,j,t,n;
  double P;
  int **psi;
  double **delta;
  double *tmp1;
  int maxind;
  int *Tcum = (int *)malloc((Nseq+1)*sizeof(int));

  Tcum[0]=0;
  for(n=1;n < Nseq+1;n++) Tcum[n]= T[n-1]+Tcum[n-1];

  psi = (int **)alloc_matrix(K,Tcum[Nseq],sizeof(int));
  delta = (double **)alloc_matrix(K,Tcum[Nseq],sizeof(double));
  tmp1 = (double *)malloc(K*sizeof(double));

  for(n=0;n<Nseq;n++) {
    t=Tcum[n];
    for(i=0;i<K;i++) {
      delta[i][t] = lpi[i] + lb[i];
      psi[i][t]=0;
    }
    for(t=Tcum[n]+1; t < Tcum[n+1];t++) {
      for(j=0;j<K;j++) {
	i=0;
	maxind = i;
	tmp1[i] = delta[i][t-1] + la[i*K+j];
	for(i=1;i<K;i++) {
	  tmp1[i] = delta[i][t-1] + la[i*K+j];
	  if(tmp1[i] > tmp1[maxind]) maxind = i;
	}
	delta[j][t] = tmp1[maxind] + lb[t*K+j];
	psi[j][t] = maxind;
      }
    }
  }

  P = 1;
  *loglik = 0.0;
  for(n=1;n<=Nseq;n++) {
    maxind = 0;
    for(i=1;i<K;i++)
      if ( delta[i][Tcum[n]-1] > delta[maxind][Tcum[n]-1] ) 
	    maxind=i;
    *loglik += delta[maxind][Tcum[n]-1];
    q[Tcum[n]-1] = maxind;
    P += delta[maxind][Tcum[n]-1];
  }

  for(n=0;n < Nseq;n++) {
    for(t=Tcum[n+1]-2; t >= Tcum[n]; t--) {
      if(q[t+1]<0) {
	     error("Invalid state at n = %d and t = %d\n",n,t+1);
	     free_matrix(K,Tcum[Nseq],(void **)psi);
	     free_matrix(K,Tcum[Nseq],(void **)delta);
      }
      else q[t] = psi[q[t+1]][t+1];
    }
  }
  //	print_imatrix2(K,T,psi);
  //	print_matrix2(K,3050,delta);
  free_matrix(K,Tcum[Nseq],(void **)psi);
  free_matrix(K,Tcum[Nseq],(void **)delta);
  free(tmp1);
  free(Tcum);
}
