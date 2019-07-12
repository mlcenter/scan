#include <math.h>
#include <stdio.h>

#ifndef  __IDM__
#define  __IDM__

#ifdef __cplusplus
extern "C" {
#endif 

extern int  HBblocks(char * fname);
extern int  LiLithF(char*fname);
extern int loopGamma(double * cs_gz, double *cs_gg);

#ifdef __cplusplus
}
#endif 


#endif
