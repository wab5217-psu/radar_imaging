#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>

#include "image.h"


int pinv_c(lapack_complex_double *in_mat, MKL_INT m, MKL_INT n, lapack_complex_double *ret_mat,int print){

  int jm,jn;
  
  MKL_INT lda=n;
  MKL_INT ldv=n;
  MKL_INT ldu=m;

  FILE *svfp;
    
  lapack_complex_double alpha;
  alpha.real=1;
  alpha.imag=0;
  lapack_complex_double beta;
  beta.real=0;
  beta.imag=0;
    
  lapack_complex_double *a;
  a=(lapack_complex_double *)mkl_calloc(lda*m,sizeof(lapack_complex_double),64);
  for( jm=0; jm<m; jm++ )for( jn=0; jn<n; jn++ )a[jm*n+jn]=in_mat[jm*n+jn];

  double *sva;
  sva=(double *)calloc(n,sizeof(double));
    
  lapack_complex_double *u;
  u=(lapack_complex_double *)mkl_calloc(ldv*m,sizeof(lapack_complex_double),64);

  lapack_complex_double *vt;
  vt=(lapack_complex_double *)mkl_calloc(ldv*n,sizeof(lapack_complex_double),64);
  
  int info=LAPACKE_zgesdd(LAPACK_ROW_MAJOR,'A',m,n,a,lda,sva,u,ldu,vt,ldv);

  if( info!=0 ){
    fprintf(stderr,"INFO -- nonzero\n");
    return(-1);
  }

  /* double sv_thresh=5e-2*sva[0]; */
  double sv_thresh=1e-2*sva[0];

  print=0;
  if(print==1)for( jn=0; jn<n; jn++ )fprintf(stderr,"%d %le\n",jn,sva[jn]);

  if( sva[0] > 1.0 ){
    svfp=fopen("singVals","w");
    fprintf(svfp,"%d\n",MIN(m,n));
    for(jn=0; jn<MIN(m,n); jn++)fprintf(svfp," %8.4le \n",sva[jn]);
    fclose(svfp);
  }
  
  /* for( jn=0; jn<n; jn++ )if( sva[jn] > sv_thresh )sva[jn]=1./sva[jn]; */
  /*   else sva[jn]=0; */

  for( jn=0; jn<n; jn++ )sva[jn]=sva[jn]/(sva[jn]*sva[jn]+sv_thresh*sv_thresh);

  int indx;
  for( jm=0; jm<m; jm++ )for( jn=0; jn<n; jn++ ){
      indx=jn*n+jm;
      vt[indx].real*=sva[jn];
      vt[indx].imag*=sva[jn];
    }

  cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasConjTrans,n,m,n,&alpha,vt,n,u,n,&beta,ret_mat,n);
  
  mkl_free(a);
  mkl_free(u);
  mkl_free(vt);
  free(sva);
  
  return(0);
}
