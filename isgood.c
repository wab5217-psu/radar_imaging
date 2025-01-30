#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image.h"

int isgood(double val){
  if( isnan(val) ) return 0;
  if( val == BAD_VALUE ) return 0;
  return 1;
    }
