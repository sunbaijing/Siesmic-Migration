#include <stdio.h>
//#include <assert.h>
#include <stdlib.h>
//#include <malloc.h>
//#include <sys/time.h>
//#include <time.h>
//#include <errno.h>
#include <math.h>

typedef unsigned int uint;

void fwd_step2D( float *, float *, float *, int *, float * );

void bc2Dx( float *, float *, float *, float *, 
            float *, int *, float *, float *, 
            float *, float * );

void bc2Dz( float *, float *, float *, float *, 
            float *, int *, float *, float *, 
            float *, float * );

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define TRANSLATE2D()\
    int stencil = prms[0];\
    int nz = prms[1];\
    int nx = prms[2];\
    int ny = prms[3];\
    int col = prms[9];\
    int col2 = prms[10];\
    int col3 = prms[11];\
    int col4 = prms[12];\
    int plane = prms[4];\
    int plane2 = prms[5];\
    int plane3 = prms[6];\
    int plane4 = prms[7];\
    int pln = prms[8];\
    int aux = prms[13];\
    float c0 = cfs[0];\
    float cx1 = cfs[1];\
    float cx2 = cfs[2];\
    float cx3 = cfs[3];\
    float cx4 = cfs[4];\
    float cy1 = cfs[9];\
    float cy2 = cfs[10];\
    float cy3 = cfs[11];\
    float cy4 = cfs[12];\
    float cz1 = cfs[5];\
    float cz2 = cfs[6];\
    float cz3 = cfs[7];\
    float cz4 = cfs[8];\
    float dt = cfs[15];

#define S1O()\
    tmp2 = tmp1 * ( acfs[0] * ( po[ indxe         ]                       ) +\
                    acfs[1] * ( po[ indxe - 1     ] + po[ indxe + 1     ] ) +\
                    acfs[2] * ( po[ indxe - col   ] + po[ indxe + col   ] ) +\
                    acfs[3] * ( po[ indxe - plane ] + po[ indxe + plane ] ) );
