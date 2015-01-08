#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* struct containing the parameters necessary for computations on the parts */
typedef struct {
  int linenumX, linenumY; /* geometries of the aperture part */
  int numofelems, numofresults; /* number of active elements and actuators respectively */
  int *lineX_start; /* start of the X chains */
  int *lineX_length; /* length of the X chains */
  int *lineX_starta; /* start of the common chain with the previous chain (X) */
  int *lineX_startb; /* start of the common chain with the next chain (X) */
  int *lineX_enda; /* end of the common chain with the previous chain (X) */
  int *lineX_endb; /* end of the common chain with the next chain (X) */
  int *actX_start; /* start of valid actuator positions, seen from X direction */
  int *actX_length; /* length of valid actuator positions, seen from X direction */
  int *lineY_start; /* start of the Y chains */
  int *lineY_length; /* length of the Y chains */
  int *lineY_starta; /* start of the common chain with the previous chain (Y) */
  int *lineY_startb; /* start of the common chain with the next chain (Y) */
  int *lineY_enda; /* end of the common chain with the previous chain (Y) */
  int *lineY_endb; /* end of the common chain with the next chain (Y) */
  int *I_sub; /* aperture mask (as integer field of zeros and ones) */
  int *empty; /* number of empty lines (x_start, x_end, y_start, y_end) */
  int *connect_length; /* length of connecting parts */
  int *connect_pos; /* positions of the actuators used for the connection */
  float *act_weight; /* weighting for the actuators */
  int *actconnect; /* gives the index in the aperture for the local actuator */
} p_par;

typedef struct {
  float *linesX; /* values of x chains */
  float *linesY; /* values of y chains */
  float **linestart; /* array of starts of y chains */
  float **extrastart; /* array of starts of x extra arrays (for alternative interpolation) */
  float *extraX; /* values of averaged gradients for alternative interpolation (X) */
  float *extraY; /* values of averaged gradients for alternative interpolation (Y) */
  float *meantX; /* mean values of X chains */
  float *meantY; /* mean values of Y chains */
  float *meanaX; /* mean values of X chains on common part with previous chain */
  float *meanaY; /* mean values of Y chains on common part with previous chain */
  float *meanbX; /* mean values of X chains on common part with next chain */
  float *meanbY; /* mean values of Y chains on common part with next chain */
  float *avgaX; /* average values of gradients on common part with previous chain (X) */
  float *avgaY; /* average values of gradients on common part with previous chain (Y) */
  float *avgbX; /* average values of gradients on common part with next chain (X) */
  float *avgbY; /* average values of gradients on common part with next chain (Y) */
  float *trendX; /* values of trend steps (X) */
  float *trendY; /* values of trend steps (Y) */
} p_arrays;

typedef struct {
  int numofelems; /* active subapertures */
  int numofresults; /* active actuators */
  int linenum; /* subapertures across the pupil */
  int ndivs; /* subdivision level in CuReD */
  int *I_sub; /* mask of active subapertures */
  int *I_fried; /* mask of actuators according to the Fried geometry */
  int *I_act; /* mask of active actuators ( = I_fried) */
  float *dataX; /* array holding measurements (X), needed to have a fixed adress for prepro step */
  float *dataY; /* array holding measurements (Y) */

} sysCure;

typedef struct {
  p_par *parts;
  p_arrays *arrays;
  int *S_connect;
  float *shift;
  float *actweight;
  float *result;
  float **S_p_iter;
  float **Sx_p;
  float **Sy_p;
  float **result_p;
  float **connect;
} parCure;

/* system setup - needs yao_params.txt for the system description */
/* creates the system struct sysCure */

sysCure* cureSystem(int linenum, int numofelems, int numofresults, int *I_sub, int ndivs);

/* creating structs needed too run CuReD */
parCure* cureInit(sysCure *sys);

/* call of CuReD */
int cured(sysCure* sys, parCure *par, float *data, float *result_vec, float *ttX=NULL, float* ttY=NULL);


void curefree(sysCure* sys, parCure *par);
