#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/types.h>

#include "image.h"
#include "mkl_bm_form.h"


    /*
      assume signal is of the form |s|=A/(t^2+w^2);

      estimate_amp(acf, sigma, lags, glags, &amp, &wid_l);

    */
    
int estimate_amp(MKL_Complex16 *acf, double *sigma, int *lags, int glags, double *amp, double *wid_l){

  int jl;
  double *pwr;
  pwr=(double *)calloc(sizeof(double),glags);
  for( jl=0; jl<glags; jl++ )
    pwr[jl]=acf[jl].real*acf[jl].real+acf[jl].imag*acf[jl].imag;

  *wid_l=44;
  
  double ybar=0;
  double xbar=0;
  int cnt=0;
  double alpha;
  double beta;
  for( jl=0; jl<glags/2; jl++ ){
    ybar+=pwr[jl];
    xbar+=lags[jl];
    cnt++;
  }
  ybar/=(double)cnt;
  xbar/=(double)cnt;

  double num=0;
  double denom=0;
  for( jl=0; jl<glags/2; jl++ ){
    num+=(lags[jl]-xbar)*(pwr[jl]-ybar);
    denom+=(lags[jl]-xbar)*(lags[jl]-xbar);
  }
  beta=num/denom;
  alpha=ybar-beta*xbar;

  if( beta>1. || isnan(alpha) || isnan(beta) )return(-1);
  //  fprintf(stderr,"%lf %lf\n",beta,alpha);

  double p0=alpha/2;
  
  int bad_cnt=0;
  int bthresh=4;
  for( jl=3; jl<glags-2; jl++ ){
    if( pwr[jl]>p0*2.4 )bad_cnt++;
  }
  if( bad_cnt>bthresh)return(-1);
  
  for( jl=3; jl<glags-2; jl++ ){
    if( (pwr[jl-1]+pwr[jl]+pwr[jl+1])/3 < p0) break;
  }					       
  *wid_l=(double)lags[jl];
  *amp=sqrt((pwr[0]+pwr[1])/2.);
  *amp=sqrt(alpha);  
    
  return(0);
}
