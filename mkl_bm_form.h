#ifndef C_BM_MAT
#define C_BM_MAT
int calculate_beam_matrix(double a0, double a1, double dang, double asep, double k_r, int beam, lapack_complex_double *bm_mat);
#endif

#ifndef MKL_ZGESVJ
#define MKL_ZGESVJ
lapack_int LAPACKE_zgesvj (int matrix_layout, char joba, char jobu, char jobv, lapack_int m, lapack_int n, lapack_complex_double * a, lapack_int lda, double * sva, lapack_int mv, lapack_complex_double * v, lapack_int ldv, double * stat);
#endif
  
#ifndef EST_AMP
#define EST_AMP
int estimate_amp(lapack_complex_double *acf, double *sigma, int *lags, int glags, double *amp, double *wid_l);
#endif

#ifndef EST_VEL
#define EST_VEL
int estimate_vel(lapack_complex_double *acf, double *sigma, int *lags, int glags, double amp, double wid_l, int tfreq, int mpinc, double *vel, double *err);
#endif
