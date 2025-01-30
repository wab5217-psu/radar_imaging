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

int calculate_bm_matrix(double ang_min, double ang_max, double dang, double asep, double k_r, MKL_Complex16 bm_mat[] ){
  double dphi;
  lapack_int nangs=(int)((ang_max-ang_min)/dang);
  lapack_int nant=MAIN_ANT;
  int ja0,ja1,jang;
  
  dang*=(double)DTOR;
  ang_min*=(double)DTOR;
  ang_max*=(double)DTOR;

  double *stheta;
  stheta=(double *)malloc(sizeof(double)*nangs);
  for( jang=0; jang<nangs; jang++)stheta[jang]=sin(ang_min+(double)jang*dang+dang/2.);
  
  for(ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nant; ja1++ ){
      //      dphi=0.95*k_r*asep*(double)(ja0-ja1);
      dphi=k_r*asep*(double)(ja0-ja1);
      for( jang=0; jang<nangs; jang++ ){ 
	bm_mat[(ja0*nant+ja1)*nangs+jang].real=cos(dphi*stheta[jang]);
	bm_mat[(ja0*nant+ja1)*nangs+jang].imag=sin(dphi*stheta[jang]);
      }
    }

  free(stheta);
  
  return(0);
}
