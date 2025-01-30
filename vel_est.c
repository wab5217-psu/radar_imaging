#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include "image.h"

#define VERBOSE 0

lapack_complex_double *dvec;

int corr_to_lag_v_ang( struct CORR_MTX corm, int nangs, int jl, lapack_complex_double **lags_v_ang, int print, lapack_complex_double *bm_mat ){

  int stat=0;
  int ja0,ja1,jang,jang1;
  int nant=MAIN_ANT;

  MKL_INT m=nangs;
  MKL_INT n=nangs;
  MKL_INT k=nant*nant;

  MKL_INT lda=k;
  MKL_INT ldb=m;
  MKL_INT ldc=m;
  
  lapack_complex_double alpha;
  alpha.real=1;
  alpha.imag=0;
  lapack_complex_double beta;
  beta.real=0;
  beta.imag=0;
  double var;

  lapack_complex_double *bm_mat_T;
  lapack_complex_double *aa;
  
  aa=(lapack_complex_double *)mkl_malloc(sizeof(lapack_complex_double)*nangs*(nangs+1),64); 
  bm_mat_T=(lapack_complex_double *)mkl_malloc(sizeof(lapack_complex_double)*nangs*(nant*nant+1),64);

  for( ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nant; ja1++ ){
      for( jang=0; jang<nangs; jang++ ){
  	bm_mat_T[jang*nant*nant+ja0+ja1*nant].real=0;
  	bm_mat_T[jang*nant*nant+ja0+ja1*nant].imag=0;
  	var=sqrt(pow(corm.i_var[ja0][ja1],2)+pow(corm.q_var[ja0][ja1],2));
  	if( var!=0. ){
  	  bm_mat_T[jang*nant*nant+ja0+ja1*nant].real=bm_mat[(ja0*nant+ja1)*nangs+jang].real/var;
  	  bm_mat_T[jang*nant*nant+ja0+ja1*nant].imag=-bm_mat[(ja0*nant+ja1)*nangs+jang].imag/var;
  	}
      }
    }

  if( bm_mat == NULL )fprintf(stderr,"bm_mat null\n");
  if( bm_mat_T == NULL )fprintf(stderr,"bm_mat_t null\n");
  if( aa == NULL )fprintf(stderr,"aa null\n");

  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,&alpha,bm_mat_T,lda,bm_mat,ldb,&beta,aa,ldc);
  
  fflush(stdout);
  stat=pinv_c(aa,n,n,aa,print);
  if( stat==-1 ) fprintf(stdout,"stat -1 pinv\n");
  if( stat==-1 ) return(-1);      
  
  for( jang=0; jang<nangs; jang++ ){
    dvec[jang].real=0;
    dvec[jang].imag=0;
    for( ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nant; ja1++ ){
    	dvec[jang].real+=bm_mat_T[jang*nant*nant+ja0*nant+ja1].real*corm.i_array[ja0][ja1] -
    	  bm_mat_T[jang*nant*nant+ja0*nant+ja1].imag*corm.q_array[ja0][ja1];
    	dvec[jang].imag+=bm_mat_T[jang*nant*nant+ja0*nant+ja1].real*corm.q_array[ja0][ja1] +
    	  bm_mat_T[jang*nant*nant+ja0*nant+ja1].imag*corm.i_array[ja0][ja1];
      }
        
  }
  
  for( jang=0; jang<nangs; jang++ ){
    lags_v_ang[jl][jang].real=0;
    lags_v_ang[jl][jang].imag=0;
    for( jang1=0; jang1<nangs; jang1++ ){
      lags_v_ang[jl][jang].real+=aa[jang*nangs+jang1].real*dvec[jang1].real -
	aa[jang*nangs+jang1].imag*dvec[jang1].imag;
      lags_v_ang[jl][jang].imag+=aa[jang*nangs+jang1].real*dvec[jang1].imag +
	aa[jang*nangs+jang1].imag*dvec[jang1].real;
    }
  }
  mkl_free(aa);
  mkl_free(bm_mat_T);
  return(0);
}


int main(int argc, char *argv[]){
  
  FILE *fp;
  FILE *fp_info;
  char fname[1024];
  char info_file[64];
  char r_chan;
  struct IMAGE_STR dstr;
  int stat;
  int j,jj;
  int beam=0;  
  int glags;
  int stid=-1;
  int printinfo=0;

  int c;
  extern char *optarg;
  extern int optind;
  
  while ((c = getopt(argc, argv, "b:is")) != -1) {
    switch (c) {
    case 'b':
      beam = atoi(optarg);
      break;

    case 'i':
      printinfo=1;
      break;

    case 's':
      stid=1;
      break;
    }
  }
      
      
  strcpy(fname,argv[optind]);
  fprintf(stderr,"%s",fname);
  fp=fopen(fname,"r");

  r_chan=fname[(strlen(fname)-1)];
  fprintf(stderr," %c\n",r_chan);

  char radar[16];
  strcpy(radar,argv[optind+1]);

  if( printinfo ){
    sprintf(info_file,"info_%d",beam);
    fp_info=fopen(info_file,"w");
  }

  
  dstr.ppat=NULL;
  dstr.pcode=NULL;
  dstr.antennas=NULL;

  for(j=0; j<N_ANT; j++){
    dstr.i_array[j]=(double**)calloc(sizeof(double*),MAX_SEQ);
    for(jj=0; jj<MAX_SEQ; jj++){
      dstr.i_array[j][jj]=(double*)calloc(sizeof(double),MAX_SAMP);
    }
    dstr.q_array[j]=(double**)calloc(sizeof(double*),MAX_SEQ);
    for(jj=0; jj<MAX_SEQ; jj++){
      dstr.q_array[j][jj]=(double*)calloc(sizeof(double),MAX_SAMP);
    }
    if( (dstr.i_array[j] == NULL) || (dstr.q_array[j] == NULL)){
      fprintf(stderr,"Can't allocate I&Q arrays\n"); return(-1);
    }
  }


  int not_beam=1;
  for( j=0; j<=beam; j++){
    if( j==beam )not_beam=0;
    if((stat=rd_image(fp, radar, not_beam, stid, &dstr))<0){
      if( j==beam ){
	fprintf(stderr,"FILE READ ERROR\n");
	exit(-1);
      }
    }
  }

  int jr,jant,jl,jang,seq;
  int nseq=dstr.sequences;    

  int nant=MAIN_ANT;
  double asep;

  if( dstr.nbaud>1 ){
    fprintf(stderr,"coded:\n");
    for( jl=0; jl<dstr.nbaud; jl++ )fprintf(stderr," %d ",dstr.pcode[jl]);
    fprintf(stderr,"\n samples: %d\n",dstr.samples);
    double *decoded;
    decoded=(double *)calloc(sizeof(double),dstr.samples);
    for( seq=0; seq<nseq; seq++ )for( jant=0; jant<nant; jant++){
      for( j=dstr.nbaud/2; j<dstr.samples-dstr.nbaud; j++ ){
  	decoded[j]=0.;
  	for( jl=0; jl<dstr.nbaud; jl++ ){
	  if(dstr.i_array[jant][seq][j+jl-dstr.nbaud/2]==BAD_VALUE) continue;	  
	  decoded[j]+=dstr.i_array[jant][seq][j+jl-dstr.nbaud/2]*dstr.pcode[jl];
	}
      }
      for( j=dstr.nbaud/2; j<dstr.samples-dstr.nbaud/2; j++ )
  	/* dstr.i_array[jant][seq][j]=decoded[j]; */
  	dstr.i_array[jant][seq][j]=decoded[j]/(double)dstr.nbaud;

      for( j=dstr.nbaud/2; j<dstr.samples-dstr.nbaud; j++ ){
  	decoded[j]=0.;
  	for( jl=0; jl<dstr.nbaud; jl++ ){
	  if(dstr.q_array[jant][seq][j+jl-dstr.nbaud/2]==BAD_VALUE) continue;	  
	  decoded[j]+=dstr.q_array[jant][seq][j+jl-dstr.nbaud/2]*dstr.pcode[jl];
	}
      }
      for( j=dstr.nbaud/2; j<dstr.samples-dstr.nbaud/2; j++ )
  	/* dstr.q_array[jant][seq][j]=decoded[j]; */
  	dstr.q_array[jant][seq][j]=decoded[j]/(double)dstr.nbaud;
    }
    free(decoded);
  }

  
  if( strcmp(radar,"kod")==0 ){
    asep=15.24;
    if(printinfo)fprintf(fp_info,"kod %d %d %d %d %d %d %d\n",dstr.year,dstr.month,dstr.day,dstr.hour,dstr.minut,dstr.tfreq,dstr.beam);
  }else if( strcmp(radar,"sps")==0 ){
    asep=15.24;
    if(printinfo)fprintf(fp_info,"sps %d %d %d %d %d %d %d\n",dstr.year,dstr.month,dstr.day,dstr.hour,dstr.minut,dstr.tfreq,dstr.beam);
  }else if( strcmp(radar,"mcm")==0 ){
    asep=12.8016;
    if(printinfo)fprintf(fp_info,"mcm %d %d %d %d %d %d %d\n",dstr.year,dstr.month,dstr.day,dstr.hour,dstr.minut,dstr.tfreq,dstr.beam);
  }else if( strcmp(radar,"bks")==0 ){
    asep=12.8016;
    if(printinfo)fprintf(fp_info,"bks %d %d %d %d %d %d %d\n",dstr.year,dstr.month,dstr.day,dstr.hour,dstr.minut,dstr.tfreq,dstr.beam);
  }else{
    fprintf(stderr,"Must specify radar: kod, mcm, sps, or bks\n");
    exit(-1);
  }
  

  long start_time=(long)((double)dstr.year*1.e8+(double)dstr.month*1.e6+(double)dstr.day*1.e4+(double)dstr.hour*100+dstr.minut);
  fprintf(stderr,"start time: %ld\n",start_time);


  double dang=DANG;
  double ang_min=ANG_MIN;
  double ang_max=ANG_MAX;
  int nangs=(int)((ang_max-ang_min)/dang);
  int nlags=dstr.ppat[dstr.mppul-1]+1;

  int *selections;
  selections=(int *)calloc(sizeof(int),nangs);
    
  double *vel_v_ang[MAX_RANGE];
  for( jr=0; jr<MAX_RANGE; jr++ )vel_v_ang[jr]=(double *)calloc(sizeof(double),nangs);

  double k_r=2.*PI*(double)dstr.tfreq/3.e5;
  
  
  struct CORR_MTX corm;
  corm.i_array=(double **)calloc(sizeof(double *),nant);
  corm.q_array=(double **)calloc(sizeof(double *),nant);
  for( jl=0; jl<nant; jl++){
    corm.i_array[jl]=(double *)calloc(sizeof(double),nant);
    corm.q_array[jl]=(double *)calloc(sizeof(double),nant);
  }
  corm.i_var=(double **)calloc(sizeof(double *),nant);
  corm.q_var=(double **)calloc(sizeof(double *),nant);
  for( jl=0; jl<nant; jl++){
    corm.i_var[jl]=(double *)calloc(sizeof(double),nant);
    corm.q_var[jl]=(double *)calloc(sizeof(double),nant);
  }

  int nback=4;
  struct CORR_MTX x_corm;
  x_corm.i_array=(double **)calloc(sizeof(double *),nant);
  x_corm.q_array=(double **)calloc(sizeof(double *),nant);
  for( jl=0; jl<nant; jl++){
    x_corm.i_array[jl]=(double *)calloc(sizeof(double),nback*nlags);
    x_corm.q_array[jl]=(double *)calloc(sizeof(double),nback*nlags);
  }
  x_corm.i_var=(double **)calloc(sizeof(double *),nant);
  x_corm.q_var=(double **)calloc(sizeof(double *),nant);
  for( jl=0; jl<nant; jl++){
    x_corm.i_var[jl]=(double *)calloc(sizeof(double),nback*nlags);
    x_corm.q_var[jl]=(double *)calloc(sizeof(double),nback*nlags);
  }
  
  k_r=2.*PI*(double)dstr.tfreq/3.e5;
  beam=dstr.beam;

  double amp;
  double vel;
  double err;
  double wid_l;
  
  lapack_complex_double **lags_v_ang;
  lags_v_ang=(lapack_complex_double **)mkl_calloc(nlags,sizeof(lapack_complex_double *),64);
  for( jl=0; jl<nlags; jl++ )
    lags_v_ang[jl]=(lapack_complex_double *)mkl_calloc(nangs,sizeof(lapack_complex_double),64);

  lapack_complex_double *acf;
  acf=(lapack_complex_double *)mkl_calloc(nlags,sizeof(lapack_complex_double),64);
  
  double *sigma;
  sigma=(double *)calloc(sizeof(double),nlags);

  double *sigma_l;
  sigma_l=(double *)calloc(sizeof(double),nlags);

  int *lags;
  lags=(int *)calloc(sizeof(int),nlags);
  
  double ssi,vi,vq;

  /* MKL_INT m=nant*nant; */
  /* MKL_INT n=nangs; */
  lapack_complex_double *bm_mat; 

  
  bm_mat=(lapack_complex_double *)mkl_malloc(sizeof(lapack_complex_double)*(nangs+1)*nant*nant,64); 
  dvec=(lapack_complex_double *)mkl_malloc(sizeof(lapack_complex_double)*nangs,64);

  int lag=1;
  int ja0,ja1;

  FILE *outf;
  FILE *out_v;
  FILE *outx;
  FILE *outc;
  
  double mean_sv,mean_count;
  
  char ang_file[64];
  char vel_ang_file[64];

  
  if( dstr.beam < 10){
    sprintf(ang_file,"%ld.0%d.angs",start_time,dstr.beam);
    sprintf(vel_ang_file,"%ld.0%d.vel_v_ang",start_time,dstr.beam);
  }else{
    sprintf(ang_file,"%ld.%d.angs",start_time,dstr.beam);
    sprintf(vel_ang_file,"%ld.%d.vel_v_ang",start_time,dstr.beam);
  }

  outf=fopen(ang_file,"w");
  out_v=fopen(vel_ang_file,"w");
  outx=fopen("xcorl","w");
  outc=fopen("corl_rng","w");

  fprintf(stderr,"beam: %d\n",dstr.beam);
  
  fprintf(out_v,"%d %d %d %d %d %d\n",dstr.year,dstr.month,dstr.day,dstr.hour,dstr.minut,dstr.sec);
  fprintf(out_v,"%d %d %d\n",dstr.nrang,nangs,dstr.smsep);
  if(printinfo){
    fprintf(fp_info,"nrang %d\n",dstr.nrang);
    fprintf(fp_info,"nang %d\n",nangs);
  }
  fprintf(out_v,"%d %d %d\n",dstr.beam,dstr.tfreq,dstr.mpinc);
  
  k_r=2.*PI*(double)dstr.tfreq/3.e5;
  beam=dstr.beam;

  stat=calculate_bm_matrix(ang_min, ang_max, dang, asep, k_r, bm_mat);
  int print=0;
    
  for( jr=0; jr<dstr.nrang-4; jr++ ){
    if(printinfo)fprintf(fp_info,"range %d\n",jr);
    print=0;
    
    glags=0;

    mean_sv=0;
    mean_count=0;
    for( jl=0; jl<nlags; jl++ ){
      for( jang=0; jang<nangs; jang++ ){
	lags_v_ang[jl][jang].real=BAD_VALUE;
	lags_v_ang[jl][jang].imag=BAD_VALUE;
      }
      
      stat=corr_lag_mtx(&dstr,jl,jr+4,&corm);

      stat=x_corr_lag_mtx(&dstr,jl,jr+4,&x_corm);

      
      ssi=0;
      glags=0;
      for( ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nant; ja1++ )if( corm.i_var[ja0][ja1]!=0. ){
	  vi=corm.i_var[ja0][ja1];
	  vq=corm.q_var[ja0][ja1];
	  ssi+=sqrt(vi*vi+vq*vq)/2.;
	  glags++;
	}
      if( glags>8 ){
	sigma_l[jl]=ssi/(double)glags;
	mean_sv+=sigma_l[jl];
	mean_count+=1;
      }else{ sigma_l[jl]=1.e10; }
      
      if( stat==-1 && printinfo==1)fprintf(fp_info,"stat -1 corr\n");
      if( stat==-1 ) continue;
      if(printinfo)fprintf(fp_info,"stat 0 corr lag %d %d\n",jl,lag);
      
      if( jl==lag && printinfo==1){
	fprintf(fp_info,"corm %d %d\n",jr,lag);
	for( ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nant; ja1++ ){
	      fprintf(fp_info,"%d %d %12.10le %12.10le %12.10le %12.10le\n",ja0,ja1,corm.i_array[ja0][ja1],corm.q_array[ja0][ja1],corm.i_var[ja0][ja1],corm.q_var[ja0][ja1]);
	    }
      }
      if( jr==20)print=1;
      stat= corr_to_lag_v_ang( corm, nangs, jl, lags_v_ang, print, bm_mat);      
      print=0;
    }

    mean_sv/=mean_count;

    if( jr==40 ){
      for( jl=0; jl<nlags; jl++ )for( jang=0; jang<nangs; jang++ ){	  
	  fprintf(outc,"%d %d %12.10le %12.10le %12.10le\n",jl,jang,lags_v_ang[jl][jang].real,lags_v_ang[jl][jang].imag,sigma_l[jl]);
	}
    }
    
    if(printinfo){
      fprintf(fp_info,"lag_v_ang %d\n",jr);
      for( jl=0; jl<nlags; jl++ )for( jang=0; jang<nangs; jang++ ){
	  fprintf(fp_info,"%d %d %12.10le %12.10le %12.10le\n",jl,jang,lags_v_ang[jl][jang].real,lags_v_ang[jl][jang].imag,sigma_l[jl]);
	}
    }
    
    if( VERBOSE ) fprintf(stderr,"range: %d ",jr);
    stat=select_data(lags_v_ang,nlags,nangs,selections);
    for( jang=0; jang<nangs; jang++ )if( selections[jang]==1 )fprintf(outf,"%d %d\n",jr,jang);

    
    for( jang=0; jang<nangs; jang++ )if( selections[jang]==1 ){       
	glags=0;
	/* for( jl=0; jl<nlags; jl++ )if(sigma_l[jl] != 0){ */
	
	for( jl=0; jl<nlags; jl++ ){
	  if( (isgood(lags_v_ang[jl][jang].real) == 1) && (sigma_l[jl] <= 1.25*mean_sv) && ( sigma_l[jl] != 0)){
	    acf[glags].real=lags_v_ang[jl][jang].real;
	    acf[glags].imag=lags_v_ang[jl][jang].imag;
	    sigma[glags]=sigma_l[jl];
	    lags[glags]=jl;
	    glags++;
	  }
	}
	stat=estimate_amp(acf, sigma, lags, glags, &amp, &wid_l);
	if( stat==-1 )continue;
	stat=estimate_vel(acf,sigma,lags,glags,amp,wid_l,dstr.tfreq,dstr.mpinc,&vel,&err);
	fprintf(out_v,"%d %d %lf %lf %lf %lf\n",jr,jang,(double)jang*DANG+ang_min,amp,3.e10/(wid_l*(double)dstr.tfreq*(double)dstr.mpinc),vel);
      }

    fprintf(outx,"range=%d\n",jr);
    for(jl=0; jl<nback*nlags; jl++){
      for( jant=0; jant<nant; jant++ )fprintf(outx,"%12.10le %12.10le ",x_corm.i_array[jant][jl],x_corm.q_array[jant][jl]);
      fprintf(outx,"\n");
    }
  }    

  fclose(outf);
  fclose(out_v);
  fclose(outx);
  fclose(outc);
  
}
