#define NREC 20
#define N_ANT 20
#define MAIN_ANT 16
#define MAX_RANGE 600
#define MAX_SAMP 1000
#define PI acos(-1.)
#define DTOR PI/180.
#define C0 3.e5
#define MAX_SEQ 50

#define ANG_MIN -30.
#define ANG_MAX 30.
#define DANG 1.

#define VEL_MIN -1500.
#define VEL_MAX 1500.
#define DVEL 25.



#define BAD_VALUE -99999.

#ifndef _MKL_H_
#include "mkl.h"
#endif

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))



struct IMAGE_STR{
  float version;
  int stid;
  int channel;
  int year;
  int month;
  int day;
  int hour;
  int minut;
  int sec;
  int nsec;
  int nrang;
  int mpinc;
  int smsep;
  int lagfr;
  int pplen;
  int beam;
  int tfreq;
  int mppul;
  int *ppat;
  int nbaud;
  int *pcode;
  int samples;
  int sequences;
  int n_antennas;
  int *antennas;
  int seq_t_sec[MAX_SEQ];
  int seq_t_msec[MAX_SEQ];
  double **i_array[N_ANT];
  double **q_array[N_ANT];
};

struct CORR_MTX{
  double **i_array;
  double **q_array;
  double **i_var;
  double **q_var;
};
  
struct img_header{
  int year;
  int month;
  int day;
  int hour;
  int minut;
  int minut0;
  int sec;
  int sec0;
  int beam;
  int freq;
  int nf;
  int nsamp;
};

#define PATH_LEN 512
#define F_NAME_LEN 64

#define IS_LEAPYEAR(year) ((year%4)&&(!(year%100))||(year%400))


typedef struct finfo{
  char dir_path[PATH_LEN];
  char fname[F_NAME_LEN];
}FILE_INFO;

/* #ifndef GSL */
/* #define GSL */
/* #include <gsl/gsl_math.h>  */
/* #include <gsl/gsl_eigen.h> */
/* #include <gsl/gsl_complex.h> */
/* #include <gsl/gsl_complex_math.h> */
/* #include <gsl/gsl_rng.h> */
/* #endif */

#ifndef TIME
#define TIME
#include <time.h>
#endif

#ifndef SEL_F
#define SEL_F
FILE_INFO *select_file(char *radar, time_t time);
#endif

#ifndef RD_IMAGE
#define RD_IMAGE
int rd_image(FILE *fp, char *radar, int not_beam, int stid,  struct IMAGE_STR *dstr);
#endif

  
#ifndef DECODE
#define DECODE
int decode(double* decoded, double* coded, int len, int nbaud, int brkr[]);
#endif

/* #ifndef CORRMTX */
/* #define CORRMTX */
/* int corrmtx(gsl_vector_complex *s,int s_dim, gsl_matrix_complex *corm); */
/* #endif */

#ifndef CORR_STR
#define CORR_STR
int corr_mtx(struct IMAGE_STR *dstr, int range, struct CORR_MTX *corm);
#endif

#ifndef PINV_STR 
#define PINV_STR
int pinv_c(lapack_complex_double *in_mat, lapack_int m, lapack_int n, lapack_complex_double *ret_mat,int print);
#endif

#ifndef CALC_STR
#define CALC_STR
int calculate_bm_matrix(double ang_min, double ang_max, double dang, double asep, double k_r, lapack_complex_double bm_mat[] );
#endif

#ifndef CORR_LAG
#define CORR_LAG
int corr_lag_mtx(struct IMAGE_STR *dstr,int lag, int range, struct CORR_MTX *corm);
#endif

#ifndef XCORR_LAG
#define XCORR_LAG
int x_corr_lag_mtx(struct IMAGE_STR *dstr,int lag, int range, struct CORR_MTX *x_corm);
#endif

#ifndef SEL_DATA
#define SEL_DATA
int select_data(lapack_complex_double **lags_v_ang, int nlags, int nangs, int *selections);
#endif

#ifndef EST_AMP
#define EST_AMP
int estimate_amp(lapack_complex_double *acf, double *sigma, int *lags, int glags, double *amp, double *wid_l);
#endif

#ifndef EST_VEL
#define EST_VEL
int estimate_vel(lapack_complex_double *acf, double *sigma, int *lags, int glags, double amp, double wid_l, int tfreq, int mpinc, double *vel, double *err_min);
#endif

#ifndef IS_GOOD
#define ID_GOOD
int isgood(double val);
#endif


