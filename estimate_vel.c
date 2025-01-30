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

int estimate_vel(MKL_Complex16 *acf, double *sigma, int *lags, int glags, double amp, double wid_l, int tfreq, int mpinc, double *vel, double *err_min){

  double v_min=-1500.;
  double v_max=1500.;
  double dv=2.5;
  double lag;
  double err=0;
  int nv=(v_max-v_min)/dv;
  int jv,jlag;
  
  double m_acf_r;
  double m_acf_i;
  double arg;
  double amp_l;

  *err_min=99999999.;
  *vel=BAD_VALUE;
  for( jv=0; jv<nv; jv++ ){
    arg=-4.*PI*(v_min+(double)jv*dv)*((double)tfreq/3.e5)*((double)mpinc/1.e6);
    err=0;
    for( jlag=0; jlag<glags; jlag++ ){
      lag=(double)lags[jlag];
      amp_l=amp*sqrt(wid_l/(lag*lag+wid_l*wid_l));
      m_acf_r=amp_l*cos(lag*arg);
      m_acf_i=amp_l*sin(lag*arg);
      err+=(pow((acf[jlag].real-m_acf_r),2)+pow((acf[jlag].imag-m_acf_i),2))/sigma[jlag];
      //      err+=(fabs(acf[jlag].real-m_acf_r)+fabs(acf[jlag].imag-m_acf_i))*sigma[jlag];
    }
    //        fprintf(stdout,"%d %d %lf %lf %lf\n",tfreq,mpinc,v_min+(double)jv*dv,err,*err_min);
    if( err<*err_min ){
      *err_min=err;
      *vel=v_min+(double)jv*dv;
    }
  }
  return(0);
}
