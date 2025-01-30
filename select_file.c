#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "image.h"

#define CNV_TO_STR(num,str) if(num < 10){sprintf(str,"0%d",num);} else{sprintf(str,"%d",num);}


FILE_INFO *select_file(char *radar, time_t time)
{
  char yr_str[5], mo_str[3], dy_str[3];
  char *raid_path;
  char dir_path[PATH_LEN];
  time_t ftime;
  time_t diftime;
  time_t mindif=100000;
  DIR *dp;
  struct dirent *ep;
  char *tz;
  extern FILE_INFO *file;
  
  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  raid_path=getenv("RAID_PATH");
  struct tm *in_time;
  in_time=gmtime(&time);
  int yr=in_time->tm_year+1900;
  int mo=in_time->tm_mon+1;
  int dy=in_time->tm_mday;
  int hr=in_time->tm_hour;
  int mn=in_time->tm_min;
  int sc=in_time->tm_sec;

  CNV_TO_STR(yr,yr_str);
  CNV_TO_STR(mo,mo_str);
  CNV_TO_STR(dy,dy_str);
  
  sprintf(file->dir_path,"%s",NULL);
  sprintf(file->fname,"%s",NULL);

  sprintf(dir_path,"%s%s/%s.%s",raid_path,yr_str,mo_str,dy_str);
  //  fprintf(stderr,"\n%s%s/%s.%s\n",raid_path,yr_str,mo_str,dy_str);
  if((dp=opendir(dir_path))==NULL) {
    fprintf(stderr,"---- COULDN'T OPEN DATA DIRECTORY ----\n");
    return(file);
  }

  while( (ep=readdir(dp)) )
    {
      if( strstr(ep->d_name,".gz") != NULL) continue;
      if( strstr(ep->d_name,".bz2") != NULL) continue;
      if( strstr(ep->d_name,radar) != NULL)
	{
	  ftime=fname_to_time(ep->d_name);
	  diftime=time-ftime;
	  if( diftime>=0 && diftime<mindif )
	  {
	    sprintf(file->dir_path,"%s",dir_path);
	    sprintf(file->fname,"%s",ep->d_name);
	    mindif=diftime;
	  }
	}
    }
  return file;
}
