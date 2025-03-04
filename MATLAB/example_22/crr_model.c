/* Conceptual rainfall-runoff model with adaptive-step explicit 
 * Runge-Kutta mwSizeegrator 
 ---written by G Schoups, Nov 2012
 ---modified by J.A. Vrugt, Jan 2013
*/

#include "mex.h"                                                                                                                                                 
#include <math.h>
#include <stdlib.h>                                                                                                                                              
#include <stdio.h>       
#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

static void runge_kutta(mwSize nvar, mwSize nt, double *tout, double *y0, const mxArray *data, const mxArray *options, double *y);
static void rk2(const mxArray *data, mwSize nvar, mwSize s, double t, double h, double *u, double *LTE);
static mwSize  fRhs(mwSize s, double t, double *u, double *udot, const mxArray *data);
static double expFlux(double Sr, double a);
static double exponen(double x);

void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{ 
  double *tout, *y0, *y;
  const mxArray *options, *data;
  mwSize nvar, nt;
  
  /* Create C pomwSizeers to input arguments */
  tout    = mxGetPr(prhs[0]);
  y0      = mxGetPr(prhs[1]);
  data    = prhs[2];
  options = prhs[3];
    
  /* Allocate memory for output arguments */
  nvar = mxGetM(prhs[1]);
  nt   = mxGetN(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(nvar,nt,mxREAL);

  /* Create C pomwSizeers to output arguments */
  y = mxGetPr(plhs[0]);

  /* Call the C subroutine */
  runge_kutta(nvar,nt,tout,y0,data,options,y);
  
  return;
}

static void runge_kutta(mwSize nvar, mwSize nt, double *tout, double *y0, const mxArray *data, const mxArray *options, double *y)
{
  double hin, hmax_, hmin_, reltol, *abstol, order;
  double h, t, t1, t2, *LTE, *ytmp, wrms, *w;
  mwSize i, s, ns, accept;

  /* Initialize */
  ns = nt - 1;
  for (i=1; i<=nvar; i++) y[i-1] = y0[i-1];
  
  /* Memory allocation for temporary vectors */
  LTE    = malloc(nvar*sizeof(double));
  ytmp   = malloc(nvar*sizeof(double));
  w      = malloc(nvar*sizeof(double));
    
  /* Extract mwSizeegration options */
  hin    = *mxGetPr(mxGetField(options,0,"InitialStep"));
  hmax_   = *mxGetPr(mxGetField(options,0,"MaxStep"));
  hmin_   = *mxGetPr(mxGetField(options,0,"MinStep"));
  reltol = *mxGetPr(mxGetField(options,0,"RelTol"));
  abstol =  mxGetPr(mxGetField(options,0,"AbsTol"));
  order  = *mxGetPr(mxGetField(options,0,"Order"));
  
  /* mwSizeegrate from tout[0] to tout[end] */
  for (s=1; s<=ns; s++) {
      /* Set start and end times */
      t1 = tout[s-1];
      t2 = tout[s];
      /* Set initial step */
      h = hin;
      h = max(hmin_,min(h,hmax_));
      h = min(h,t2-t1);
      /* Set initial y */
      for (i=1; i<=nvar; i++) y[i-1+nvar*s] = y[i-1+nvar*(s-1)];
      /* mwSizeegrate from t1 to t2 */
      t = t1;
      while (t < t2) {
        /* Advance solution by step h */
        for (i=1; i<=nvar; i++) ytmp[i-1] = y[i-1+nvar*s];
        rk2(data,nvar,s,t,h,ytmp,LTE);
        /* Decide whether to accept current step */
    	accept = 0;
        wrms = 0;
        for (i=1; i<=nvar; i++) {
            w[i-1] = 1.0/(reltol*fabs(ytmp[i-1]) + abstol[i-1]);
            wrms = wrms + pow(w[i-1]*LTE[i-1],2);
        }
        wrms = pow(wrms/nvar,0.5);
        if (wrms <= 1) accept = 1;
    	/* If accepted, update y and t */
        if (accept > 0) {
    		for (i=1; i<=nvar; i++) y[i-1+nvar*s] = ytmp[i-1];
            t = t + h;
    	}
        /* Compute new step */
        h = h*max(0.2,min(5.0,0.9*pow(wrms,-1.0/order)));
        h = max(hmin_,min(h,hmax_));
        h = min(h,t2-t);
      }
  }
  
  /* Free memory for temporary vectors */
  free(LTE);
  free(ytmp);
  free(w);

  return;
}

static void rk2(const mxArray *data, mwSize nvar, mwSize s, double t, double h, double *u, double *LTE)
{
	mwSize i, flag;
    double *udotE, *uE, *udot;

	/* Memory allocation for temporary vectors */
    udotE = malloc(nvar*sizeof(double));
    uE    = malloc(nvar*sizeof(double));
    udot  = malloc(nvar*sizeof(double));
    
	/* Euler solution */
	flag = fRhs(s,t,u,udotE,data);
    for (i=1; i<=nvar; i++) uE[i-1] = u[i-1] + h*udotE[i-1];
	/* Heun solution */
	flag = fRhs(s,t+h,uE,udot,data);
    for (i=1; i<=nvar; i++) u[i-1] = u[i-1] + 0.5*h*(udotE[i-1] + udot[i-1]);    
	/* Compute estimate of LTE */
    for (i=1; i<=nvar; i++) LTE[i-1] = fabs(uE[i-1]-u[i-1]);

	/* Free memory for temporary vectors */
    free(udotE);
    free(uE);
    free(udot);
   
    return;
}

/* Conceptual rainfall-runoff model */
static mwSize fRhs(mwSize s, double t, double *u, double *udot, const mxArray *data)
{
    double *P, *Ep, Imax, Sumax, Qsmax, aE, aF, aS, Kf, Ks;
    double Si, Su, Sf, Ss;
    double Precip, Evap, Perc, Runoff, FastQ, SlowQ;
    double EvapI, P_e, Ep_e;

    P     =  mxGetPr(mxGetField(data,0,"P"));
    Ep    =  mxGetPr(mxGetField(data,0,"Ep"));
    Imax  = *mxGetPr(mxGetField(data,0,"Imax"));
    Sumax = *mxGetPr(mxGetField(data,0,"Sumax"));
    Qsmax = *mxGetPr(mxGetField(data,0,"Qsmax"));
    aE    = *mxGetPr(mxGetField(data,0,"aE"));
    aF    = *mxGetPr(mxGetField(data,0,"aF"));
    aS    = *mxGetPr(mxGetField(data,0,"aS"));
    Kf    = *mxGetPr(mxGetField(data,0,"Kf"));
    Ks    = *mxGetPr(mxGetField(data,0,"Ks"));

    Si = u[0];
    Su = u[1];
    Sf = u[2];
    Ss = u[3];

    Precip = P[s-1];
    if(Imax > 0.0) {
        EvapI  = Ep[s-1]*expFlux(Si/Imax,50.0);
        P_e    = P[s-1]*expFlux(Si/Imax,-50.0);
        Ep_e   = max(0.0,Ep[s-1]-EvapI);
    } else {
        EvapI = 0.0;
        P_e  = P[s-1];
        Ep_e = Ep[s-1];
    }
    Evap   = Ep_e * expFlux(Su/Sumax,aE);
    Perc   = Qsmax * expFlux(Su/Sumax,aS);
    Runoff = P_e * expFlux(Su/Sumax,aF);
    FastQ  = Sf/Kf;
    SlowQ  = Ss/Ks;

    udot[0] = Precip - EvapI - P_e;
    udot[1] = P_e - Evap - Perc - Runoff;
    udot[2] = Runoff - FastQ;
    udot[3] = Perc - SlowQ;
    udot[4] = FastQ + SlowQ;
/*
    if (!(FastQ >= 0) || !(SlowQ >= 0)){
	mexPrmwSizef("%s%f\n", "Precip ",Precip);
	mexPrmwSizef("%s%f\n", "Runoff ", Runoff);
	mexPrmwSizef("%s%f\n", "Evap ",Evap);
	mexPrmwSizef("%s%f\n", "Perc ",Perc);
	mexPrmwSizef("%s%f\n","aF ",aF);
	mexPrmwSizef("%s%f\n","Su ",Su);
	mexPrmwSizef("%s%f\n","Ss ",Ss);
    mexPrmwSizef("%s%f\n","Sf ",Sf);
    mexPrmwSizef("%s%f\n","Kf ",Kf);
	mexPrmwSizef("%s%f\n","FastQ ",FastQ);
	mexPrmwSizef("%s\n", "----");
    return(0) ;
	}
*/
return(0);
}

/* Relative flux from exponential storage-flux relation */
static double expFlux(double Sr, double a)
{
    double Qr;
    
    Sr = max(0.0,min(1.0,Sr));
    if(fabs(a) < 1e-6) {
        Qr = Sr;    /* approximately linear */
    } else {
        Qr = (1.-exponen(-a*Sr))/(1.-exponen(-a));
    }
            
    return Qr;
}

/* Exponential function with protection against overflow */
static double exponen(double x)
{
    double f;
    
    f = exp(min(300.0,x));
            
    return f;
}
