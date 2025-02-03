#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "image.h"


int rd_image(FILE *fp, char *radar, int not_beam, int stid, struct IMAGE_STR *dstr){
  int jjj,jj,j,ja,js;
  size_t bytes;

  if((bytes=fread(&(dstr->version),1,sizeof(float),fp))<sizeof(float))return(-1);
  if( stid>=0 ){
    if((bytes=fread(&(dstr->stid),1,sizeof(int),fp))<sizeof(int))return(-1);
    if((bytes=fread(&(dstr->channel),1,sizeof(int),fp))<sizeof(int))return(-1);
  }
  if((bytes=fread(&(dstr->year),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->month),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->day),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->hour),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->minut),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->sec),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->nsec),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->nrang),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->mpinc),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->smsep),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->lagfr),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->pplen),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->beam),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->tfreq),1,sizeof(int),fp))<sizeof(int))return(-1);
  
  if((bytes=fread(&(dstr->mppul),1,sizeof(int),fp))<sizeof(int))return(-1);
  if( dstr->ppat != NULL )free(dstr->ppat);
  dstr->ppat=(int *)calloc(sizeof(int),dstr->mppul);
  for( jj=0; jj<dstr->mppul; jj++)if((bytes=fread(&(dstr->ppat[jj]),1,sizeof(int),fp))<sizeof(int))return(-1);
  
  if((bytes=fread(&(dstr->nbaud),1,sizeof(int),fp))<sizeof(int))return(-1);
  if( dstr->pcode != NULL )free(dstr->pcode);
  dstr->pcode=(int *)calloc(sizeof(int),dstr->nbaud);
  for( jj=0; jj<dstr->nbaud; jj++)if((bytes=fread(&(dstr->pcode[jj]),1,sizeof(int),fp))<sizeof(int))return(-1);
  
  if((bytes=fread(&(dstr->samples),1,sizeof(int),fp))<sizeof(int))return(-1);
  if((bytes=fread(&(dstr->sequences),1,sizeof(int),fp))<sizeof(int))return(-1);
  
  if((bytes=fread(&(dstr->n_antennas),1,sizeof(int),fp))<sizeof(int))return(-1);
  if( dstr->antennas != NULL )free(dstr->antennas);
  dstr->antennas=(int *)calloc(sizeof(int),dstr->n_antennas);
  for(jj=0; jj<dstr->n_antennas; jj++)if((bytes=fread(&(dstr->antennas[jj]),1,sizeof(int),fp))<sizeof(int))return(-1);
  
  float isamp;
  float qsamp;
  
  for( j=0; j<dstr->sequences; j++){
    if((bytes=fread(&(dstr->seq_t_sec[j]),1,sizeof(int),fp))<sizeof(int))return(-1);
    if((bytes=fread(&(dstr->seq_t_msec[j]),1,sizeof(int),fp))<sizeof(int))return(-1);
    for( jj=0; jj<dstr->n_antennas; jj++)for(jjj=0; jjj<dstr->samples; jjj++ ){
	
	if((bytes=fread(&isamp,1,sizeof(float),fp))<sizeof(float))return(-1);
	if((bytes=fread(&qsamp,1,sizeof(float),fp))<sizeof(float))return(-1);
	if(isnan(isamp) || isnan(qsamp)){
	  fprintf(stderr,"********BAD SAMPLES*******\n");
	  return(-1);
	}else{
	  dstr->i_array[dstr->antennas[jj]][j][jjj]=(double)isamp;
	  dstr->q_array[dstr->antennas[jj]][j][jjj]=(double)qsamp;
	}
      }
  }
  
  if( not_beam )return(0);
  
  double *mask;
  double samp_p_tau=(double)dstr->mpinc/(double)dstr->smsep;
  int s1;
  int sbuf;
  
  sbuf=3;

  fprintf(stderr,"nbaud: %d\n",dstr->nbaud);

  mask=(double *)calloc(sizeof(double),dstr->samples);
  for( j=0; j<dstr->samples; j++ )mask[j]=1.;
  for( j=0; j<dstr->mppul; j++ ){
    s1=dstr->ppat[j]*samp_p_tau;
    for( js=MAX(s1-2,0); js<MIN(s1+sbuf*dstr->nbaud,dstr->samples); js++ )mask[js]=0;
  }
  
  for( ja=0; ja<N_ANT; ja++ ){
    for( js=0; js<dstr->sequences; js++ )for( j=0; j<dstr->samples; j++ ){
	if( dstr->i_array[ja][js][j] == BAD_VALUE ) continue;
	if( dstr->q_array[ja][js][j] == BAD_VALUE ) continue;
	if( mask[j] == 0 ){
	  dstr->i_array[ja][js][j]=BAD_VALUE;
	  dstr->q_array[ja][js][j]=BAD_VALUE;
	  continue;
	}
      }
  }
  
  free(mask);
  return(0);
}
