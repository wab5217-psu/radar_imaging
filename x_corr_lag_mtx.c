#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "image.h"

int x_corr_lag_mtx(struct IMAGE_STR *dstr,int lag, int range, struct CORR_MTX *x_corm){

  int j,j1,seq,ja0,ja1,ja_lag,jb;

  int nant=MAIN_ANT;
  int nback=4;
  int mppul=dstr->mppul;
  int nlags=dstr->ppat[mppul-1]+1;
  int nseq=dstr->sequences;    
  double samp_p_tau=(double)dstr->mpinc/(double)dstr->smsep;
  double ss1;
  double ss2;
  double taper;
  
  
  int *lag_table[2];
  for( j=0; j<2; j++)lag_table[j]=(int *)calloc(sizeof(int),nlags);
  for( j1=0; j1<nlags; j1++ )lag_table[0][j1]=-1;
  for( j1=0; j1<nlags; j1++ )lag_table[1][j1]=-1;

  lag_table[0][0]=0;
  lag_table[1][0]=0;
  for( j=0; j<mppul-1; j++) for( j1=j+1; j1<mppul; j1++){
      lag_table[0][abs(dstr->ppat[j]-dstr->ppat[j1])]=dstr->ppat[j];
      lag_table[1][abs(dstr->ppat[j]-dstr->ppat[j1])]=dstr->ppat[j1];
    }

  if( lag_table[0][lag] == -1 )return(-1);
  
  int s0=lag_table[0][lag]*samp_p_tau;
  int s1=lag_table[1][lag]*samp_p_tau;

  int buf0=1;
  int buf1=2;    
  if( dstr->nbaud>1 ){
    buf0=dstr->nbaud;
    buf1=2*dstr->nbaud;
  }
  int pulse_range=0;
  int bad=0;
  for( j=0; j<mppul; j++){
    pulse_range=dstr->ppat[j]*samp_p_tau;

    if( range+s0 >= pulse_range-buf0 && range+s0 <= pulse_range+buf1 ) bad=1;
    if( range+s1 >= pulse_range-buf0 && range+s1 <= pulse_range+buf1 ) bad=1;
  }

  if( bad == 1 )return(-1);

  long seq_cnt=0;
  double data_mx=1000.;
  double a0,a1;
  double taper_width=64.;
  
  for( ja0=0; ja0<nant; ja0++ )for(ja1=0; ja1<nback; ja1++){
      jb=nant+ja1;
      ja_lag=ja1+lag*nback;
      x_corm->i_array[ja0][ja_lag]=0;
      x_corm->q_array[ja0][ja_lag]=0;
      seq_cnt=0;
      for( seq=0; seq<nseq; seq++ ){	    
	if( isgood(dstr->i_array[ja0][seq][range+s0]) &&
	    isgood(dstr->q_array[ja0][seq][range+s0]) &&
	    isgood(dstr->i_array[jb][seq][range+s1]) &&
	    isgood(dstr->q_array[jb][seq][range+s1]) &&
	    fabs(dstr->i_array[ja0][seq][range+s0])<=data_mx &&
	    fabs(dstr->i_array[jb][seq][range+s1])<=data_mx &&
	    fabs(dstr->q_array[ja0][seq][range+s0])<=data_mx &&
	    fabs(dstr->q_array[jb][seq][range+s1])<=data_mx){

	  x_corm->i_array[ja0][ja_lag]+=
	    (dstr->i_array[ja0][seq][range+s0]*dstr->i_array[jb][seq][range+s1]
	     +dstr->q_array[ja0][seq][range+s0]*dstr->q_array[jb][seq][range+s1]);
	  x_corm->q_array[ja0][ja_lag]+=
	    (-dstr->i_array[ja0][seq][range+s0]*dstr->q_array[jb][seq][range+s1]
	     +dstr->q_array[ja0][seq][range+s0]*dstr->i_array[jb][seq][range+s1]);
	  seq_cnt++;
	  
	}
      }
      if( seq_cnt!=0 ){
	x_corm->i_array[ja0][ja_lag]/=(double)(seq_cnt);
	x_corm->q_array[ja0][ja_lag]/=(double)(seq_cnt);
      }      
    }	
  
   
  for( ja0=0; ja0<nant; ja0++ )for( ja1=0; ja1<nback; ja1++ ){
      jb=nant+ja1;
      ja_lag=ja1+lag*nback;
      x_corm->i_var[ja0][ja_lag]=0;
      x_corm->q_var[ja0][ja_lag]=0;
      a0=(double)ja0-7.5;
      a1=(double)ja_lag-7.5;
      taper=1.;
      taper=exp(-(a0*a0/taper_width+a1*a1/taper_width));
      seq_cnt=0;
      for( seq=0; seq<nseq; seq++ ){
	if( isgood(dstr->i_array[ja0][seq][range+s0]) &&
	    isgood(dstr->q_array[ja0][seq][range+s0]) &&
	    isgood(dstr->i_array[jb][seq][range+s1]) &&
	    isgood(dstr->q_array[jb][seq][range+s1]) &&
 	    fabs(dstr->i_array[ja0][seq][range+s0])<=data_mx &&
	    fabs(dstr->i_array[jb][seq][range+s1])<=data_mx &&
	    fabs(dstr->q_array[ja0][seq][range+s0])<=data_mx &&
	    fabs(dstr->q_array[jb][seq][range+s1])<=data_mx){

	  ss1=
	    (dstr->i_array[ja0][seq][range+s0]*dstr->i_array[jb][seq][range+s1]
	     +dstr->q_array[ja0][seq][range+s0]*dstr->q_array[jb][seq][range+s1])*taper;
	  ss2=
	    (-dstr->i_array[ja0][seq][range+s0]*dstr->q_array[jb][seq][range+s1]
	     +dstr->q_array[ja0][seq][range+s0]*dstr->i_array[jb][seq][range+s1])*taper;
	  x_corm->i_var[ja0][ja_lag]+=pow(ss1-x_corm->i_array[ja0][ja_lag],2);
	  x_corm->q_var[ja0][ja_lag]+=pow(ss2-x_corm->q_array[ja0][ja_lag],2);
	  seq_cnt++;
	}
      }
      if( seq_cnt!=0 ){
	x_corm->i_var[ja0][ja_lag]=x_corm->i_var[ja0][ja_lag]/(double)(seq_cnt-1);
	x_corm->q_var[ja0][ja_lag]=x_corm->q_var[ja0][ja_lag]/(double)(seq_cnt-1);
      }
    }
  free(lag_table[0]);
  free(lag_table[1]);
  return(0);
}
