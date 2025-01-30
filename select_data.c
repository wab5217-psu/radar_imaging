#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/types.h>
#include "image.h"

#define MIN_RATIO 5.0
#define MIN_PT_RATIO 2.0
#define MIN_POWER 1.0
#define MIN_GOOD 10

#define VERBOSE 0

int select_data(MKL_Complex16 **lags_v_ang, int nlags, int nangs, int *selections){

  double pwr;
  int jl,ja;

  int buffer=3;
  double mxpwr;
  double avgpwr;
  double goodcnt;
  
  for( ja=0; ja<nangs; ja++ )selections[ja]=0;
  goodcnt=0;
  avgpwr=0;
  mxpwr=0;
  
  if( isgood(lags_v_ang[0][(int)(nangs/2)].real) ){
    jl=0;
  }else{
    jl=1;
  }
  
  for( ja=buffer; ja<nangs-buffer; ja++ )
    if( isgood(lags_v_ang[jl][ja].real) ){
      pwr=lags_v_ang[jl][ja].real*lags_v_ang[jl][ja].real+
	lags_v_ang[jl][ja].imag*lags_v_ang[jl][ja].imag;
      if( pwr>mxpwr )mxpwr=pwr;
    }

  for( ja=buffer; ja<nangs-buffer; ja++ )
    if( isgood(lags_v_ang[jl][ja].real) ){
      pwr=lags_v_ang[jl][ja].real*lags_v_ang[jl][ja].real+
	lags_v_ang[jl][ja].imag*lags_v_ang[jl][ja].imag;
      if( pwr < mxpwr/MIN_RATIO ){
	avgpwr+=pwr;
	goodcnt++;
      }
    }
  

  if( goodcnt < MIN_GOOD ){
    if( VERBOSE )fprintf(stderr,"min good failed\n");
    return(0);
  }
    
  avgpwr/=goodcnt;

  if( VERBOSE )fprintf(stderr,"mx: %lf avg: %lf\n",mxpwr,avgpwr);
  
  if( mxpwr<MIN_POWER ){
    if( VERBOSE )fprintf(stderr,"mxpwr failed\n");
    return(0);
  }
  if(mxpwr/avgpwr < MIN_RATIO ){
    if( VERBOSE )fprintf(stderr,"ratio failed\n");
    return(0);
  }
  
  /* fprintf(stderr,"select data: %lf %lf %lf\n",mxpwr,avgpwr,mxpwr/avgpwr); */
  
  for( ja=0; ja<nangs; ja++ ){
    if( isgood(lags_v_ang[jl][ja].real) && lags_v_ang[jl][ja].real>0. ){
      pwr=lags_v_ang[jl][ja].real*lags_v_ang[jl][ja].real+
	lags_v_ang[jl][ja].imag*lags_v_ang[jl][ja].imag;
      if((pwr>MIN_POWER) && (pwr/mxpwr>1.e-2) && (pwr/avgpwr > MIN_PT_RATIO)){
	  selections[ja]=1;
	}
    }
  }
  
  
  return(0);
}
