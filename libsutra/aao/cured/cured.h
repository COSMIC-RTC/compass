#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* struct containing the parameters necessary for computations on the parts */
typedef struct {
  int32_t linenumX, linenumY; /* geometries of the aperture part */
  int32_t numofelems,
      numofresults;  /* number of active elements and actuators respectively */
  int32_t *lineX_start;  /* start of the X chains */
  int32_t *lineX_length; /* length of the X chains */
  int32_t *lineX_starta; /* start of the common chain with the previous chain (X) */
  int32_t *lineX_startb; /* start of the common chain with the next chain (X) */
  int32_t *lineX_enda;   /* end of the common chain with the previous chain (X) */
  int32_t *lineX_endb;   /* end of the common chain with the next chain (X) */
  int32_t *
      actX_start; /* start of valid actuator positions, seen from X direction */
  int32_t *actX_length; /* length of valid actuator positions, seen from X direction
                     */
  int32_t *lineY_start; /* start of the Y chains */
  int32_t *lineY_length; /* length of the Y chains */
  int32_t *lineY_starta; /* start of the common chain with the previous chain (Y) */
  int32_t *lineY_startb; /* start of the common chain with the next chain (Y) */
  int32_t *lineY_enda;   /* end of the common chain with the previous chain (Y) */
  int32_t *lineY_endb;   /* end of the common chain with the next chain (Y) */
  int32_t *I_sub;        /* aperture mask (as integer field of zeros and ones) */
  int32_t *empty; /* number of empty lines (x_start, x_end, y_start, y_end) */
  int32_t *connect_length; /* length of connecting parts */
  int32_t *connect_pos;    /* positions of the actuators used for the connection */
  float *act_weight;   /* weighting for the actuators */
  int32_t *actconnect; /* gives the index in the aperture for the local actuator */
} p_par;

typedef struct {
  float *linesX;      /* values of x chains */
  float *linesY;      /* values of y chains */
  float **linestart;  /* array of starts of y chains */
  float **extrastart; /* array of starts of x extra arrays (for alternative
                         interpolation) */
  float *extraX; /* values of averaged gradients for alternative interpolation
                    (X) */
  float *extraY; /* values of averaged gradients for alternative interpolation
                    (Y) */
  float *meantX; /* mean values of X chains */
  float *meantY; /* mean values of Y chains */
  float
      *meanaX; /* mean values of X chains on common part with previous chain */
  float
      *meanaY; /* mean values of Y chains on common part with previous chain */
  float *meanbX; /* mean values of X chains on common part with next chain */
  float *meanbY; /* mean values of Y chains on common part with next chain */
  float *avgaX;  /* average values of gradients on common part with previous
                    chain (X) */
  float *avgaY;  /* average values of gradients on common part with previous
                    chain (Y) */
  float *avgbX;  /* average values of gradients on common part with next chain
                    (X) */
  float *avgbY;  /* average values of gradients on common part with next chain
                    (Y) */
  float *trendX; /* values of trend steps (X) */
  float *trendY; /* values of trend steps (Y) */
} p_arrays;

typedef struct {
  int32_t numofelems;   /* active subapertures */
  int32_t numofresults; /* active actuators */
  int32_t linenum;      /* subapertures across the pupil */
  int32_t ndivs;        /* subdivision level in CuReD */
  int32_t *I_sub;       /* mask of active subapertures */
  int32_t *I_fried;     /* mask of actuators according to the Fried geometry */
  int32_t *I_act;       /* mask of active actuators ( = I_fried) */
  float *dataX; /* array holding measurements (X), needed to have a fixed adress
                   for prepro step */
  float *dataY; /* array holding measurements (Y) */

} sysCure;

typedef struct {
  p_par *parts;
  p_arrays *arrays;
  int32_t *S_connect;
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

sysCure *cureSystem(int32_t linenum, int32_t numofelems, int32_t numofresults, int32_t *I_sub,
                    int32_t ndivs);

/* creating structs needed too run CuReD */
parCure *cureInit(sysCure *sys);

/* call of CuReD */
int32_t cured(sysCure *sys, parCure *par, float *data, float *result_vec,
          float *ttX = NULL, float *ttY = NULL);

void curefree(sysCure *sys, parCure *par);
