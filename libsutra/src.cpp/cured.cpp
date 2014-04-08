#include "cured.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* CuReD for general aperture
 * author: Matthias Rosensteiner */

#define max(a,b) (((a) > (b)) ? (a) : (b)) 
#define min(a,b) (((a) < (b)) ? (a) : (b)) 
#define round(a) ((a >= 0) ? (int)(a + 0.5) : (int)(a - 0.5))

sysCure*
cureSystem() {

  int i, j;
  int linenum;
  float readme;
  sysCure *sys;
  FILE *file;

  sys = (sysCure*) malloc(sizeof(sysCure));
  file = fopen("cured_params.txt", "r");

  fscanf(file, "%e", &readme);
  sys->linenum = (int) readme;
  linenum = sys->linenum;
  fscanf(file, "%e", &readme);
  sys->numofelems = (int) readme;
  fscanf(file, "%e", &readme);
  sys->numofresults = (int) readme;

  sys->I_sub = (int*) malloc(linenum * linenum * sizeof(int));
  sys->I_fried = (int*) malloc((linenum + 1) * (linenum + 1) * sizeof(int));
  sys->I_act = sys->I_fried;

  for (i = 0; i < linenum * linenum; i++) {
    fscanf(file, "%e", &readme);
    sys->I_sub[i] = (int) readme;
  }
  fscanf(file, "%e", &readme);
  sys->ndivs = (int) readme;

  for (i = 0; i < (linenum + 1) * (linenum + 1); i++)
    sys->I_fried[i] = 0;
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[i * (linenum + 1) + j] = sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[i * (linenum + 1) + j + 1] = sys->I_fried[i * (linenum + 1)
          + j + 1] | sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[(i + 1) * (linenum + 1) + j] = sys->I_fried[(i + 1)
          * (linenum + 1) + j] | sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[(i + 1) * (linenum + 1) + j + 1] = sys->I_fried[(i + 1)
          * (linenum + 1) + j + 1] | sys->I_sub[i * linenum + j];

  return sys;
}

sysCure*
cureSystem(int linenum, int numofelems, int numofresults, int *I_sub,
    int ndivs) {

  int i, j;
  //float readme;
  sysCure *sys;

  sys = (sysCure*) malloc(sizeof(sysCure));

  sys->linenum = linenum;
  sys->numofelems = numofelems;
  sys->numofresults = numofresults;

  sys->I_sub = (int*) malloc(linenum * linenum * sizeof(int));
  sys->I_fried = (int*) malloc((linenum + 1) * (linenum + 1) * sizeof(int));
  sys->I_act = sys->I_fried;

  memcpy(sys->I_sub, I_sub, linenum * linenum * sizeof(int));

  sys->ndivs = ndivs;

  for (i = 0; i < (linenum + 1) * (linenum + 1); i++)
    sys->I_fried[i] = 0;
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[i * (linenum + 1) + j] = sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[i * (linenum + 1) + j + 1] = sys->I_fried[i * (linenum + 1)
          + j + 1] | sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[(i + 1) * (linenum + 1) + j] = sys->I_fried[(i + 1)
          * (linenum + 1) + j] | sys->I_sub[i * linenum + j];
  for (i = 0; i < linenum; i++)
    for (j = 0; j < linenum; j++)
      sys->I_fried[(i + 1) * (linenum + 1) + j + 1] = sys->I_fried[(i + 1)
          * (linenum + 1) + j + 1] | sys->I_sub[i * linenum + j];

  return sys;
}

parCure*
cureInit(sysCure *sys) {

  int i, j, k, l, start, help, help2, row, sum;
  int ndivs = sys->ndivs; /* handover of system parameters */
  int numofelems = sys->numofelems;
  int numofresults = sys->numofresults;
//int numofacts = numofresults;
  int linenum = sys->linenum;
  int *I_sub = sys->I_sub;
  int parts = 1 << ndivs;
  int mainpartsize, rest;
  int *iter, *S_iter;
  int *partsize = (int*) malloc(parts * sizeof(int));
  int *S_connect = (int*) malloc(numofelems * sizeof(int));
  int **I_sub_p_iter = (int**) malloc(parts * parts * sizeof(int*));
  float *startweight;
  float **Sx_p = (float**) malloc(parts * parts * sizeof(float*));
  float **Sy_p = (float**) malloc(parts * parts * sizeof(float*));
  float **S_p_iter = (float**) malloc(parts * parts * sizeof(float*));
//float **S_p_iter2 = (float**)malloc(parts*parts*sizeof(float*));
  float **result_p = (float**) malloc(parts * parts * sizeof(float*));
  float **result_p_iter = (float**) malloc(parts * parts * sizeof(float*));
  float **connect = (float**) malloc(parts * parts * sizeof(float*));
  float *shift = (float*) malloc(parts * parts * sizeof(float));
  float *actweight = (float*) malloc(numofresults * sizeof(float));

  p_par *part_par = (p_par*) malloc(parts * parts * sizeof(p_par));
  p_arrays *part_arrays = (p_arrays*) malloc(parts * parts * sizeof(p_arrays));

  parCure *par = (parCure*) malloc(sizeof(parCure));

  /* not implemented: CuReD(0) */
  if (ndivs == 0) {
    return 0;
  }

  /* create part sizes */
  mainpartsize = linenum / parts;
  rest = linenum - parts * mainpartsize;
  for (i = 0; i < parts; i++)
    partsize[i] = mainpartsize;
  for (i = 0; i < rest; i++)
    partsize[(parts - rest) / 2 + i] += 1;
  for (i = 0; i < parts; i++)
    for (j = 0; j < parts; j++) {
      start = i * parts + j;
      part_par[start].I_sub = (int*) malloc(
          (partsize[i]) * (partsize[j]) * sizeof(int));
      I_sub_p_iter[start] = part_par[start].I_sub;
      part_par[start].numofelems = 0;
      part_par[start].linenumX = partsize[i];
      part_par[start].linenumY = partsize[j];
    }
  /* create I_sub parts and connection to I_sub */
  iter = I_sub;
  S_iter = S_connect;
  start = 0;
  for (i = 0; i < parts; i++) {
    for (j = 0; j < partsize[i]; j++)
      for (k = 0; k < parts; k++)
        for (l = 0; l < partsize[k]; l++) {
          *(I_sub_p_iter[start + k]) = *iter;
          part_par[start + k].numofelems += *iter;
          if (*iter) {
            *S_iter = start + k;
            *iter = start + k + 1;
            S_iter++;
          }
          I_sub_p_iter[start + k]++;
          iter++;
        }
    start += parts;
  }

  for (i = 0; i < parts; i++)
    for (j = 0; j < parts; j++) {
      start = i * parts + j;
      part_par[start].numofresults = part_par[start].numofelems
          + part_par[start].linenumX + part_par[start].linenumY + 1;
    }

  /* create the parameters and fields for the parts for faster computation */
  for (k = 0; k < parts * parts; k++) {
    part_par[k].lineX_start = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].lineX_length = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].lineX_starta = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].lineX_startb = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].lineX_enda = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].lineX_endb = (int*) malloc(
        (part_par[k].linenumX) * sizeof(int));
    part_par[k].actX_start = (int*) malloc(
        (part_par[k].linenumX + 1) * sizeof(int));
    part_par[k].actX_length = (int*) malloc(
        (part_par[k].linenumX + 1) * sizeof(int));
    part_par[k].lineY_start = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].lineY_length = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].lineY_starta = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].lineY_startb = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].lineY_enda = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].lineY_endb = (int*) malloc(
        (part_par[k].linenumY) * sizeof(int));
    part_par[k].connect_length = (int*) malloc(4 * sizeof(int));
    part_par[k].empty = (int*) malloc(4 * sizeof(int));
    part_par[k].connect_pos = (int*) malloc(
        (2 * part_par[k].linenumX + 2 * part_par[k].linenumY + 4)
            * sizeof(int));
    part_par[k].act_weight = (float*) malloc(
        (part_par[k].numofresults) * sizeof(float));
    part_par[k].actconnect = (int*) malloc(
        (part_par[k].numofresults) * sizeof(int));

    part_arrays[k].linestart = (float**) malloc(
        part_par[k].linenumY * sizeof(float*)); /* array of starts of y lines */
    part_arrays[k].extrastart = (float**) malloc(
        part_par[k].linenumY * sizeof(float*));
    part_arrays[k].linesX = (float*) malloc(
        (part_par[k].linenumX + part_par[k].numofelems) * sizeof(float)); /* x lines */
    part_arrays[k].linesY = (float*) malloc(
        (part_par[k].linenumY + part_par[k].numofelems) * sizeof(float)); /* y lines */
    part_arrays[k].extraX = (float*) malloc(
        (part_par[k].linenumY + part_par[k].numofelems) * sizeof(float));
    part_arrays[k].extraY = (float*) malloc(
        (part_par[k].linenumX + part_par[k].numofelems) * sizeof(float));
    part_arrays[k].meantX = (float*) malloc(
        part_par[k].linenumX * sizeof(float)); /* total mean */
    part_arrays[k].meantY = (float*) malloc(
        part_par[k].linenumY * sizeof(float));
    part_arrays[k].meanaX = (float*) malloc(
        part_par[k].linenumX * sizeof(float)); /* mean on common part with previous */
    part_arrays[k].meanaY = (float*) malloc(
        part_par[k].linenumY * sizeof(float));
    part_arrays[k].meanbX = (float*) malloc(
        part_par[k].linenumX * sizeof(float)); /* mean on common part with next */
    part_arrays[k].meanbY = (float*) malloc(
        part_par[k].linenumY * sizeof(float));
    part_arrays[k].avgaX = (float*) malloc(
        part_par[k].linenumY * sizeof(float)); /* average on common part with previous */
    part_arrays[k].avgaY = (float*) malloc(
        part_par[k].linenumX * sizeof(float));
    part_arrays[k].avgbX = (float*) malloc(
        part_par[k].linenumY * sizeof(float)); /* average on common part with next */
    part_arrays[k].avgbY = (float*) malloc(
        part_par[k].linenumX * sizeof(float));
    part_arrays[k].trendX = (float*) malloc(
        part_par[k].linenumX * sizeof(float));
    part_arrays[k].trendY = (float*) malloc(
        part_par[k].linenumY * sizeof(float));

    if (part_par[k].numofelems == 0) {
      part_par[k].numofresults = 0;
      continue;
    }

    /* setup of geometric parameters */
    for (i = 0; i < 4; i++)
      part_par[k].empty[i] = 0;
    /* x line setup */
    iter = part_par[k].I_sub;
    help2 = part_par[k].linenumX;
    for (i = 0; i < part_par[k].linenumX; i++) {
      part_par[k].lineX_start[i] = 0;
      part_par[k].lineX_length[i] = 0;
      for (j = 0; j < part_par[k].linenumY; j++) {
        part_par[k].lineX_length[i] += iter[j];
        if (part_par[k].lineX_length[i] == 0)
          part_par[k].lineX_start[i]++;
      }
      iter += part_par[k].linenumY;
      if (part_par[k].lineX_length[i] == 0) { /* remove empty lines */
        i--;
        part_par[k].linenumX--;
        if (i < 0)
          part_par[k].empty[0]++;
        else
          part_par[k].empty[1]++;
      }
    }

    /* y line setup */
    iter = part_par[k].I_sub;
    for (i = 0; i < part_par[k].linenumY; i++) {
      part_par[k].lineY_start[i] = 0;
      part_par[k].lineY_length[i] = 0;
    }
    for (i = 0; i < help2/* old part_par[k].linenumX*/; i++) {
      for (j = 0; j < part_par[k].linenumY; j++) {
        part_par[k].lineY_length[j] += iter[j];
        if (part_par[k].lineY_length[j] == 0)
          part_par[k].lineY_start[j]++;
      }
      iter += part_par[k].linenumY;
    }
    for (i = 0; i < part_par[k].linenumY; i++) { /* remove empty lines */
      if (part_par[k].lineY_length[i] == 0) {
        for (j = i; j < part_par[k].linenumY - 1; j++) {
          part_par[k].lineY_start[j] = part_par[k].lineY_start[j + 1];
          part_par[k].lineY_length[j] = part_par[k].lineY_length[j + 1];
        }
        i--;
        part_par[k].linenumY--;
        if (i < 0)
          part_par[k].empty[2]++;
        else
          part_par[k].empty[3]++;
      }
    }

    for (i = 0; i < part_par[k].linenumX; i++)
      part_par[k].lineX_start[i] -= part_par[k].empty[2];
    for (i = 0; i < part_par[k].linenumY; i++)
      part_par[k].lineY_start[i] -= part_par[k].empty[0];
    part_par[k].numofresults = part_par[k].numofelems + part_par[k].linenumX
        + part_par[k].linenumY + 1;

    part_par[k].lineX_starta[0] = part_par[k].lineX_start[0];
    part_par[k].lineX_enda[0] = part_par[k].lineX_start[0]
        + part_par[k].lineX_length[0];
    for (i = 1; i < part_par[k].linenumX; i++) {
      part_par[k].lineX_starta[i] =
          max(part_par[k].lineX_start[i-1],part_par[k].lineX_start[i]);
      part_par[k].lineX_enda[i] =
          min(part_par[k].lineX_start[i-1]+part_par[k].lineX_length[i-1],part_par[k].lineX_start[i]+part_par[k].lineX_length[i]);
    }
    part_par[k].lineX_startb[part_par[k].linenumX - 1] =
        part_par[k].lineX_start[part_par[k].linenumX - 1];
    part_par[k].lineX_endb[part_par[k].linenumX - 1] =
        part_par[k].lineX_start[part_par[k].linenumX - 1]
            + part_par[k].lineX_length[part_par[k].linenumX - 1];
    for (i = 0; i < part_par[k].linenumX - 1; i++) {
      part_par[k].lineX_startb[i] =
          max(part_par[k].lineX_start[i+1],part_par[k].lineX_start[i]);
      part_par[k].lineX_endb[i] =
          min(part_par[k].lineX_start[i+1]+part_par[k].lineX_length[i+1],part_par[k].lineX_start[i]+part_par[k].lineX_length[i]);
    }

    part_par[k].actX_start[0] = part_par[k].lineX_start[0];
    part_par[k].actX_length[0] = part_par[k].lineX_length[0] + 1;
    for (i = 1; i < part_par[k].linenumX; i++) {
      part_par[k].actX_start[i] =
          min(part_par[k].lineX_start[i-1],part_par[k].lineX_start[i]);
      help =
          max(part_par[k].lineX_start[i-1]+part_par[k].lineX_length[i-1],part_par[k].lineX_start[i]+part_par[k].lineX_length[i]);
      part_par[k].actX_length[i] = help + 1 - part_par[k].actX_start[i];
    }
    part_par[k].actX_start[part_par[k].linenumX] =
        part_par[k].lineX_start[part_par[k].linenumX - 1];
    part_par[k].actX_length[part_par[k].linenumX] =
        part_par[k].lineX_length[part_par[k].linenumX - 1] + 1;

    startweight = part_par[k].act_weight;
    for (i = 0; i < part_par[k].actX_length[0]; i++)
      startweight[i] = 2;
    startweight += part_par[k].actX_length[0];
    for (i = 1; i < part_par[k].linenumX; i++) {
      for (j = 0; j < part_par[k].actX_length[i]; j++) {
        row = j + part_par[k].actX_start[i];
        sum = 0;
        if (part_par[k].actX_start[i - 1] <= row
            && part_par[k].actX_start[i - 1] + part_par[k].actX_length[i - 1]
                > row)
          sum++;
        if (part_par[k].actX_start[i + 1] <= row
            && part_par[k].actX_start[i + 1] + part_par[k].actX_length[i + 1]
                > row)
          sum++;
        startweight[j] = 2 * sum;
      }
      startweight += part_par[k].actX_length[i];
    }
    for (i = 0; i < part_par[k].actX_length[part_par[k].linenumX]; i++)
      startweight[i] = 2;

    part_par[k].lineY_starta[0] = part_par[k].lineY_start[0];
    part_par[k].lineY_enda[0] = part_par[k].lineY_start[0]
        + part_par[k].lineY_length[0];
    for (i = 1; i < part_par[k].linenumY; i++) {
      part_par[k].lineY_starta[i] =
          max(part_par[k].lineY_start[i-1],part_par[k].lineY_start[i]);
      part_par[k].lineY_enda[i] =
          min(part_par[k].lineY_start[i-1]+part_par[k].lineY_length[i-1],part_par[k].lineY_start[i]+part_par[k].lineY_length[i]);
    }
    part_par[k].lineY_startb[part_par[k].linenumY - 1] =
        part_par[k].lineY_start[part_par[k].linenumY - 1];
    part_par[k].lineY_endb[part_par[k].linenumY - 1] =
        part_par[k].lineY_start[part_par[k].linenumY - 1]
            + part_par[k].lineY_length[part_par[k].linenumY - 1];
    for (i = 0; i < part_par[k].linenumY - 1; i++) {
      part_par[k].lineY_startb[i] =
          max(part_par[k].lineY_start[i+1],part_par[k].lineY_start[i]);
      part_par[k].lineY_endb[i] =
          min(part_par[k].lineY_start[i+1]+part_par[k].lineY_length[i+1],part_par[k].lineY_start[i]+part_par[k].lineY_length[i]);
    }

  }

  /* setup for connection of the parts */
  for (k = 0; k < parts * parts; k++) {
    part_par[k].connect_length[0] = 0; /* set to -1 if empty lines are cut at that side */
    if (part_par[k].empty[0] || part_par[k].numofelems == 0)
      part_par[k].connect_length[0] = -1;
    part_par[k].connect_length[1] = 0;
    if (part_par[k].empty[1] || part_par[k].numofelems == 0)
      part_par[k].connect_length[1] = -1;
    part_par[k].connect_length[2] = 0;
    if (part_par[k].empty[2] || part_par[k].numofelems == 0)
      part_par[k].connect_length[2] = -1;
    part_par[k].connect_length[3] = 0;
    if (part_par[k].empty[3] || part_par[k].numofelems == 0)
      part_par[k].connect_length[3] = -1;
  }
  for (i = 1; i < parts; i++) {
    for (j = 0; j < parts; j++) {
      start = i * parts + j;
      if (part_par[start].connect_length[0] != -1
          && part_par[start - parts].connect_length[1] != -1) {
        help = part_par[start - parts].linenumX;
        sum =
            min(part_par[start].actX_start[0]+part_par[start].actX_length[0]+part_par[start].empty[2],part_par[start-parts].actX_start[help]+part_par[start-parts].actX_length[help]+part_par[start-parts].empty[2]);
        row =
            max(part_par[start].actX_start[0]+part_par[start].empty[2],part_par[start-parts].actX_start[help]+part_par[start-parts].empty[2]);
        sum -= row;
        if (sum < 0)
          sum = 0;
        part_par[start].connect_length[0] = sum;
        part_par[start - parts].connect_length[1] = sum;
        for (k = 0; k < sum; k++) {
          part_par[start].connect_pos[k] = row - part_par[start].actX_start[0]
              - part_par[start].empty[2] + k;
          part_par[start - parts].connect_pos[part_par[start - parts].linenumY
              + 1 + k] = part_par[start - parts].numofresults;
          part_par[start - parts].connect_pos[part_par[start - parts].linenumY
              + 1 + k] -= part_par[start - parts].actX_length[help];
          part_par[start - parts].connect_pos[part_par[start - parts].linenumY
              + 1 + k] += row - part_par[start - parts].actX_start[help]
              - part_par[start - parts].empty[2];
          part_par[start - parts].connect_pos[part_par[start - parts].linenumY
              + 1 + k] += k;
        }
      }
    }
  }
  for (i = 0; i < parts; i++) {
    for (j = 1; j < parts; j++) {
      start = i * parts + j;
      if (part_par[start].connect_length[2] != -1
          && part_par[start - 1].connect_length[3] != -1) {
        help = 0;
        help2 = -1;
        row = max(part_par[start].empty[0],part_par[start-1].empty[0]);
        sum = partsize[i] - row
            - max(part_par[start].empty[1],part_par[start-1].empty[1]);
        for (k = 0; k < row - part_par[start].empty[0]; k++)
          help += part_par[start].actX_length[k];
        for (k = 0; k < row - part_par[start - 1].empty[0]; k++)
          help2 += part_par[start - 1].actX_length[k];
        for (k = 0; k < sum + 1; k++) {
          help2 += part_par[start - 1].actX_length[k + row
              - part_par[start - 1].empty[0]];
          if (part_par[start].actX_start[k + row - part_par[start].empty[0]]
              == 0
              && part_par[start - 1].actX_start[k + row
                  - part_par[start - 1].empty[0]]
                  + part_par[start - 1].actX_length[k + row
                      - part_par[start - 1].empty[0]]
                  == part_par[start - 1].linenumY + 1) {
            part_par[start].connect_pos[2 * part_par[start].linenumY + 2
                + part_par[start].connect_length[2]] = help;
            part_par[start].connect_length[2]++;
            part_par[start - 1].connect_pos[2 * part_par[start - 1].linenumY
                + part_par[start - 1].linenumX + 3
                + part_par[start - 1].connect_length[3]] = help2;
            part_par[start - 1].connect_length[3]++;
          }
          help +=
              part_par[start].actX_length[k + row - part_par[start].empty[0]];
        }
      }
    }
  }

  for (i = 0; i < parts * parts; i++) {
    if (part_par[i].connect_length[0] < 1)
      part_par[i].connect_length[0] = 1;
    if (part_par[i].connect_length[1] < 1)
      part_par[i].connect_length[1] = 1;
    if (part_par[i].connect_length[2] < 1)
      part_par[i].connect_length[2] = 1;
    if (part_par[i].connect_length[3] < 1)
      part_par[i].connect_length[3] = 1;
  }

  /* create connection from results on parts to result */
  for (i = 0; i < parts * parts; i++)
    I_sub_p_iter[i] = part_par[i].actconnect;
  for (i = 0; i < numofresults; i++)
    actweight[i] = 0;
  startweight = actweight;
  start = 0;
  if (I_sub[0]) {
    (*startweight)++;
    startweight++;
    *(I_sub_p_iter[I_sub[0] - 1]++) = start;
    start++;
  }
  for (j = 0; j < linenum - 1; j++) {
    if (I_sub[j]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[j] - 1]++) = start;
    }
    if (I_sub[j + 1] != 0 && I_sub[j + 1] != I_sub[j]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[j + 1] - 1]++) = start;
    }
    if (*startweight) {
      startweight++;
      start++;
    }
  }
  if (I_sub[linenum - 1]) {
    (*startweight)++;
    *(I_sub_p_iter[I_sub[linenum - 1] - 1]++) = start;
    startweight++;
    start++;
  }
  for (i = 0; i < linenum - 1; i++) {
    if (I_sub[i * linenum]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[i * linenum] - 1]++) = start;
    }
    if (I_sub[(i + 1) * linenum] != 0
        && I_sub[(i + 1) * linenum] != I_sub[i * linenum]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[(i + 1) * linenum] - 1]++) = start;
    }
    if (*startweight) {
      startweight++;
      start++;
    }
    for (j = 0; j < linenum - 1; j++) {
      if (I_sub[i * linenum + j]) {
        (*startweight)++;
        *(I_sub_p_iter[I_sub[i * linenum + j] - 1]++) = start;
      }
      if (I_sub[i * linenum + j + 1] != 0
          && I_sub[i * linenum + j + 1] != I_sub[i * linenum + j]) {
        (*startweight)++;
        *(I_sub_p_iter[I_sub[i * linenum + j + 1] - 1]++) = start;
      }
      if (I_sub[(i + 1) * linenum + j] != 0
          && I_sub[(i + 1) * linenum + j] != I_sub[i * linenum + j]
          && I_sub[(i + 1) * linenum + j] != I_sub[i * linenum + j + 1]) {
        (*startweight)++;
        *(I_sub_p_iter[I_sub[(i + 1) * linenum + j] - 1]++) = start;
      }
      if (I_sub[(i + 1) * linenum + j + 1] != 0
          && I_sub[(i + 1) * linenum + j + 1] != I_sub[(i + 1) * linenum + j]
          && I_sub[(i + 1) * linenum + j + 1] != I_sub[i * linenum + j + 1]) {
        (*startweight)++;
        *(I_sub_p_iter[I_sub[(i + 1) * linenum + j + 1] - 1]++) = start;
      }
      if (*startweight) {
        startweight++;
        start++;
      }
    }
    if (I_sub[(i + 1) * linenum - 1]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[(i + 1) * linenum - 1] - 1]++) = start;
    }
    if (I_sub[(i + 2) * linenum - 1] != 0
        && I_sub[(i + 2) * linenum - 1] != I_sub[(i + 1) * linenum - 1]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[(i + 2) * linenum - 1] - 1]++) = start;
    }
    if (*startweight) {
      startweight++;
      start++;
    }
  }
  if (I_sub[(linenum - 1) * linenum]) {
    (*startweight)++;
    startweight++;
    *(I_sub_p_iter[I_sub[(linenum - 1) * linenum] - 1]++) = start;
    start++;
  }
  for (j = 0; j < linenum - 1; j++) {
    if (I_sub[(linenum - 1) * linenum + j]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[(linenum - 1) * linenum + j] - 1]++) = start;
    }
    if (I_sub[(linenum - 1) * linenum + j + 1] != 0
        && I_sub[(linenum - 1) * linenum + j + 1]
            != I_sub[(linenum - 1) * linenum + j]) {
      (*startweight)++;
      *(I_sub_p_iter[I_sub[(linenum - 1) * linenum + j + 1] - 1]++) = start;
    }
    if (*startweight) {
      startweight++;
      start++;
    }
  }
  if (I_sub[linenum * linenum - 1]) {
    (*startweight)++;
    startweight++;
    *(I_sub_p_iter[I_sub[linenum * linenum - 1] - 1]++) = start;
    start++;
  }

  /* create the data field for the parts */
  for (i = 0; i < parts * parts; i++) {
    Sx_p[i] = (float*) malloc(part_par[i].numofelems * sizeof(float));
    Sy_p[i] = (float*) malloc(part_par[i].numofelems * sizeof(float));
    result_p[i] = (float*) malloc(part_par[i].numofresults * sizeof(float));
    result_p_iter[i] = result_p[i];
    connect[i] = (float*) malloc(4 * sizeof(float));
  }
  for (i = 0; i < parts * parts; i++)
    for (j = 0; j < 4; j++)
      connect[i][j] = 0;

  /* put everything in the big struct to handle it over to the reconstruction function */
  par->parts = part_par;
  par->arrays = part_arrays;
  par->S_p_iter = S_p_iter;
  par->Sx_p = Sx_p;
  par->Sy_p = Sy_p;
  par->S_connect = S_connect;
  par->result_p = result_p;
  par->connect = connect;
  par->shift = shift;
  par->actweight = actweight;
  par->result = (float *) malloc(numofresults * sizeof(float));

  return par;
}

/* CuRe implementation for computing the parts */
int cure(p_par part_par, p_arrays part_arrays, float *dataX, float *dataY,
    float *result, float *connect) {

  int i, j, row, sum;
  int linenumX = part_par.linenumX;
  int linenumY = part_par.linenumY;
//int numofelems = part_par.numofelems;
  int numofresults = part_par.numofresults;
//int *I_sub = part_par.I_sub;
  int *par_lineX_start = part_par.lineX_start;
  int *par_lineX_length = part_par.lineX_length;
  int *par_lineX_starta = part_par.lineX_starta;
  int *par_lineX_startb = part_par.lineX_startb;
  int *par_lineX_enda = part_par.lineX_enda;
  int *par_lineX_endb = part_par.lineX_endb;
  int *par_actX_start = part_par.actX_start;
  int *par_actX_length = part_par.actX_length;
  int *par_lineY_start = part_par.lineY_start;
  int *par_lineY_length = part_par.lineY_length;
  int *par_lineY_starta = part_par.lineY_starta;
  int *par_lineY_startb = part_par.lineY_startb;
  int *par_lineY_enda = part_par.lineY_enda;
  int *par_lineY_endb = part_par.lineY_endb;
  float *par_act_weight = part_par.act_weight;

  float *startdataX = dataX, *startdataY = dataY;
  float *startline, *startextra, *startresult;
  float **linestart = part_arrays.linestart; /* array of starts of y lines */
  float **extrastart = part_arrays.extrastart;
  float *linesX = part_arrays.linesX; /* x lines */
  float *linesY = part_arrays.linesY; /* y lines */
  float *extraX = part_arrays.extraX;
  float *extraY = part_arrays.extraY;
  float *meantX = part_arrays.meantX; /* total mean */
  float *meantY = part_arrays.meantY;
  float *meanaX = part_arrays.meanaX; /* mean on common part with previous */
  float *meanaY = part_arrays.meanaY;
  float *meanbX = part_arrays.meanbX; /* mean on common part with next */
  float *meanbY = part_arrays.meanbY;
  float *avgaX = part_arrays.avgaX; /* average on common part with previous */
  float *avgaY = part_arrays.avgaY;
  float *avgbX = part_arrays.avgbX; /* average on common part with next */
  float *avgbY = part_arrays.avgbY;
  float *trendX = part_arrays.trendX;
  float *trendY = part_arrays.trendY;
  float mX, mY, help, help2, help3;

  /* computation of x lines */
  startdataX = dataX;
  startline = linesX;

  for (i = 0; i < linenumY; i++) {
    avgaX[i] = 0.0;
    avgbX[i] = 0.0;
  }
  for (i = 0; i < linenumX; i++) {
    meantX[i] = 0.0;
    meanaX[i] = 0.0;
    meanbX[i] = 0.0;
  }
  for (i = 0; i < linenumX; i++) {
    startline[0] = 0.0;
    for (j = 0; j < par_lineX_length[i]; j++) {
      startline[j + 1] = startline[j] + startdataX[j];
      meantX[i] += startline[j + 1];
      row = par_lineX_start[i] + j;
      if (par_lineX_starta[i] <= row && row < par_lineX_enda[i]) {
        meanaX[i] += startline[j];
        meanaX[i] += startline[j + 1];
      }
      if (par_lineX_startb[i] <= row && row < par_lineX_endb[i]) {
        meanbX[i] += startline[j];
        meanbX[i] += startline[j + 1];
      }
      if (par_lineY_starta[row] <= i && i < par_lineY_enda[row])
        avgaX[row] += startdataX[j];
      if (par_lineY_startb[row] <= i && i < par_lineY_endb[row])
        avgbX[row] += startdataX[j];
    }
    meantX[i] -= startline[par_lineX_length[i]] / 2;
    startline += par_lineX_length[i] + 1;
    startdataX += par_lineX_length[i];
  }
  for (i = 0; i < linenumX; i++) {
    meantX[i] /= par_lineX_length[i];
    meanaX[i] /= 2 * (par_lineX_enda[i] - par_lineX_starta[i]);
    meanbX[i] /= 2 * (par_lineX_endb[i] - par_lineX_startb[i]);
  }
  for (i = 0; i < linenumY; i++) {
    avgaX[i] /= 2 * (par_lineY_enda[i] - par_lineY_starta[i]);
    avgbX[i] /= 2 * (par_lineY_endb[i] - par_lineY_startb[i]);
  }

  /* computation of y lines */
  startdataY = dataY;
  for (i = 0; i < linenumY; i++) {
    meantY[i] = 0.0;
    meanaY[i] = 0.0;
    meanbY[i] = 0.0;
  }
  for (i = 0; i < linenumX; i++) {
    avgaY[i] = 0.0;
    avgbY[i] = 0.0;
  }
  linestart[0] = linesY;
  linestart[0][0] = 0.0;
  for (i = 1; i < linenumY; i++) {
    linestart[i] = linestart[i - 1] + par_lineY_length[i - 1] + 1;
    *(linestart[i]) = 0.0;
  }
  for (i = 0; i < linenumX; i++) {
    for (j = 0; j < par_lineX_length[i]; j++) {
      row = par_lineX_start[i] + j;
      linestart[row][1] = linestart[row][0] + startdataY[j];
      meantY[row] += *(linestart[row]);
      if (par_lineY_starta[row] <= i && i < par_lineY_enda[row])
        meanaY[row] += *(linestart[row]);
      if (par_lineY_startb[row] <= i && i < par_lineY_endb[row])
        meanbY[row] += *(linestart[row]);
      if (par_lineX_starta[i] <= row && row < par_lineX_enda[i])
        avgaY[i] += startdataY[j];
      if (par_lineX_startb[i] <= row && row < par_lineX_endb[i])
        avgbY[i] += startdataY[j];
      linestart[row]++;
    }

    avgaY[i] /= 2 * (par_lineX_enda[i] - par_lineX_starta[i]);
    avgbY[i] /= 2 * (par_lineX_endb[i] - par_lineX_startb[i]);
    startdataY += par_lineX_length[i];
  }
  startline = linesY;
  for (i = 0; i < linenumY; i++) {
    meantY[i] += (*(startline + par_lineY_length[i])) / 2;
    meantY[i] /= par_lineY_length[i];
    meanaY[i] -= (*(startline + par_lineY_starta[i] - par_lineY_start[i])) / 2;
    meanaY[i] += (*(startline + par_lineY_enda[i] - par_lineY_start[i])) / 2;
    meanaY[i] /= par_lineY_enda[i] - par_lineY_starta[i];
    meanbY[i] -= (*(startline + par_lineY_startb[i] - par_lineY_start[i])) / 2;
    meanbY[i] += (*(startline + par_lineY_endb[i] - par_lineY_start[i])) / 2;
    meanbY[i] /= par_lineY_endb[i] - par_lineY_startb[i];
    startline += par_lineY_length[i] + 1;
  }

  /* trend computations */
  trendX[0] = 0.0;
  trendY[0] = 0.0;
  for (i = 1; i < linenumX; i++)
    trendX[i] = trendX[i - 1] + meanbX[i - 1] - meanaX[i] + avgbY[i - 1]
        + avgaY[i];
  for (i = 1; i < linenumY; i++)
    trendY[i] = trendY[i - 1] + meanbY[i - 1] - meanaY[i] + avgbX[i - 1]
        + avgaX[i];

  /* global mean value computation */
  mX = 0.0;
  for (i = 0; i < linenumX; i++)
    mX += par_lineX_length[i] * (trendX[i] + meantX[i]);
  sum = 0;
  for (i = 0; i < linenumX; i++)
    sum += par_lineX_length[i];
  mX /= sum;
  mY = 0.0;
  for (i = 0; i < linenumY; i++)
    mY += par_lineY_length[i] * (trendY[i] + meantY[i]);
  sum = 0;
  for (i = 0; i < linenumY; i++)
    sum += par_lineY_length[i];
  mY /= sum;
  /* update computations */
  for (i = 0; i < linenumX; i++)
    trendX[i] -= mX;
  for (i = 0; i < linenumY; i++)
    trendY[i] -= mY;

  /* shift lines */
  startline = linesX;
  for (i = 0; i < linenumX; i++) {
    for (j = 0; j < par_lineX_length[i] + 1; j++)
      startline[j] += trendX[i];
    startline += par_lineX_length[i] + 1;
  }
  startline = linesY;
  for (i = 0; i < linenumY; i++) {
    for (j = 0; j < par_lineY_length[i] + 1; j++)
      startline[j] += trendY[i];
    startline += par_lineY_length[i] + 1;
  }

  /* interpolation */
  for (i = 0; i < numofresults; i++)
    result[i] = 0;
  /* alternative interpolation extras */
  startdataX = dataY;
  startextra = extraY;
  for (i = 0; i < linenumX; i++) {
    startextra[0] = startdataX[0] / 2;
    for (j = 0; j < par_lineX_length[i] - 1; j++)
      startextra[j + 1] = (startdataX[j] + startdataX[j + 1]) / 4;
    startextra[par_lineX_length[i]] = startdataX[par_lineX_length[i] - 1] / 2;
    startextra += par_lineX_length[i] + 1;
    startdataX += par_lineX_length[i];
  }
  startdataY = dataX;
  linestart[0] = extraX;
  linestart[0][0] = 0.0;
  for (i = 1; i < linenumY; i++) {
    linestart[i] = linestart[i - 1] + par_lineY_length[i - 1] + 1;
    *(linestart[i]) = 0.0;
  }
  for (i = 0; i < linenumX; i++) {
    for (j = 0; j < par_lineX_length[i]; j++) {
      row = par_lineX_start[i] + j;
      linestart[row][0] += startdataY[j];
      linestart[row][0] /= 4;
      linestart[row][1] = startdataY[j];
      linestart[row]++;
    }
    startdataY += par_lineX_length[i];
  }
  startextra = extraX;
  for (i = 0; i < linenumY; i++) {
    *startextra *= 2;
    startextra[par_lineY_length[i]] /= 2;
    startextra += par_lineY_length[i] + 1;
  }

  /* add x lines and y extras */
  startresult = result;
  startline = linesX;
  startextra = extraY;
  for (i = 0; i < linenumX; i++) {
    for (j = 0; j < par_lineX_length[i] + 1; j++) {
      *startresult += *startline;
      *startresult -= *startextra;
      startresult++;
      startline++;
      startextra++;
    }
    if (i > 0
        && par_lineX_start[i - 1] + par_lineX_length[i - 1]
            > par_lineX_start[i] + par_lineX_length[i])
      startresult += par_lineX_start[i - 1] + par_lineX_length[i - 1]
          - par_lineX_start[i] - par_lineX_length[i];
    if (i < linenumX - 1 && par_lineX_start[i] < par_lineX_start[i + 1])
      startresult += par_lineX_start[i + 1] - par_lineX_start[i];
  }
  startresult = result + par_lineX_length[0] + 1;
  startline = linesX;
  startextra = extraY;
  for (i = 0; i < linenumX; i++) {
    if (i < linenumX - 1 && par_lineX_start[i + 1] < par_lineX_start[i])
      startresult += par_lineX_start[i] - par_lineX_start[i + 1];
    for (j = 0; j < par_lineX_length[i] + 1; j++) {
      *startresult += *startline;
      *startresult += *startextra;
      startresult++;
      startline++;
      startextra++;
    }
    if (i < linenumX - 1
        && par_lineX_start[i] + par_lineX_length[i]
            < par_lineX_start[i + 1] + par_lineX_length[i + 1])
      startresult += par_lineX_start[i + 1] + par_lineX_length[i + 1]
          - par_lineX_start[i] - par_lineX_length[i];
  }
  for (i = 0; i < numofresults; i++) {
    result[i] /= par_act_weight[i];
  }
  /* add y lines and x extras */
  linestart[0] = linesY;
  extrastart[0] = extraX;
  for (i = 1; i < linenumY; i++) {
    linestart[i] = linestart[i - 1] + par_lineY_length[i - 1] + 1;
    extrastart[i] = extrastart[i - 1] + par_lineY_length[i - 1] + 1;
  }
  startresult = result;
  for (i = 0; i < linenumX + 1; i++) {
    row = par_actX_start[i];
    help2 = *(linestart[row]);
    help3 = *(extrastart[row]);
    *startresult += (help2 - help3) / 2;
    help = help2 + help3;
    linestart[row]++;
    extrastart[row]++;
    startresult++;
    for (j = 1; j < par_actX_length[i] - 1; j++) {
      row = par_actX_start[i] + j;
      help2 = *(linestart[row]);
      help3 = *(extrastart[row]);
      *startresult += (help + help2 - help3) / 4;
      help = help2 + help3;
      linestart[row]++;
      extrastart[row]++;
      startresult++;
    }
    *startresult += help / 2;
    startresult++;
  }

  /* compute connection to neighbouring parts */
  connect[0] = 0;
  if (part_par.connect_length[0] > 1) {
    connect[0] += result[part_par.connect_pos[0]] / 2;
    for (i = 1; i < part_par.connect_length[0] - 1; i++)
      connect[0] += result[part_par.connect_pos[i]];
    connect[0] += result[part_par.connect_pos[part_par.connect_length[0] - 1]]
        / 2;
  }
  connect[1] = 0;
  if (part_par.connect_length[1] > 1) {
    connect[1] += result[part_par.connect_pos[linenumY + 1]] / 2;
    for (i = 1; i < part_par.connect_length[1] - 1; i++)
      connect[1] += result[part_par.connect_pos[linenumY + 1 + i]];
    connect[1] += result[part_par.connect_pos[linenumY
        + part_par.connect_length[1]]] / 2;
  }
  connect[2] = 0;
  if (part_par.connect_length[2] > 1) {
    connect[2] += result[part_par.connect_pos[2 * linenumY + 2]] / 2;
    for (i = 1; i < part_par.connect_length[2] - 1; i++)
      connect[2] += result[part_par.connect_pos[2 * linenumY + 2 + i]];
    connect[2] += result[part_par.connect_pos[2 * linenumY + 1
        + part_par.connect_length[2]]] / 2;
  }
  connect[3] = 0;
  if (part_par.connect_length[3] > 1) {
    connect[3] += result[part_par.connect_pos[2 * linenumY + linenumX + 3]] / 2;
    for (i = 1; i < part_par.connect_length[3] - 1; i++)
      connect[3] +=
          result[part_par.connect_pos[2 * linenumY + linenumX + 3 + i]];
    connect[3] += result[part_par.connect_pos[2 * linenumY + linenumX + 2
        + part_par.connect_length[3]]] / 2;
  }

  return 0;
}

int cured(sysCure* sys, parCure *par, float *data, float *result_vec,
    float gain) {

  /* variable definitions */
  int i, j, k, l, m, start, help, len1, len2, len3, len4;
  int ndivs = sys->ndivs;
  int parts = 1 << ndivs;
  int numofelems = sys->numofelems;
  int numofresults = sys->numofresults;
  int *S_connect = par->S_connect;
  float error, v1, v2, diff1, diff2, diff3, diff4, shift1, shift2, shift3;
//float tiptiltX, tiptiltY;
  float *dataX;
  float *dataY;
  float *shift = par->shift;
  float *actweight = par->actweight;
  float *result = par->result;
  float **S_p_iter = par->S_p_iter;
  float **Sx_p = par->Sx_p;
  float **Sy_p = par->Sy_p;
  float **result_p = par->result_p;
  float **connect = par->connect;

  p_par *part_par = par->parts;
  p_arrays *part_arrays = par->arrays;

  /* connection from data array to Sx, Sy */
  dataX = &(data[0]);
  dataY = &(data[numofelems]);

  /* distribute Sx, Sy */
  for (i = 0; i < parts * parts; i++)
    S_p_iter[i] = Sx_p[i];
  for (i = 0; i < numofelems; i++) {
    *(S_p_iter[S_connect[i]]++) = dataX[i];
  }
  for (i = 0; i < parts * parts; i++)
    S_p_iter[i] = Sy_p[i];
  for (i = 0; i < numofelems; i++) {
    *(S_p_iter[S_connect[i]]++) = dataY[i];
  }

  /* compute results on parts */
  for (i = 0; i < parts * parts; i++)
    if (part_par[i].numofelems)
      cure(part_par[i], part_arrays[i], Sx_p[i], Sy_p[i], result_p[i],
          connect[i]);

  /* connect parts */
  for (i = 0; i < parts * parts; i++)
    shift[i] = 0;

  for (i = 0; i < ndivs; i++) {
    for (j = 0; j < parts; j += (2 << i))
      for (k = 0; k < parts; k += (2 << i)) {
        start = j * parts + k;
        v1 = 0;
        v2 = 0;
        len1 = 0;
        for (l = 0; l < (1 << i); l++) {
          help = start + ((1 << i) - 1) * parts;
          v1 += shift[help + l] * (part_par[help + l].connect_length[1] - 1)
              + connect[help + l][1];
          v2 += shift[help + l + parts]
              * (part_par[help + l + parts].connect_length[0] - 1)
              + connect[help + l + parts][0];
          len1 += (part_par[help + l].connect_length[1] - 1);
        }
        if (len1 < 1)
          len1 = -1;
        diff1 = (v1 - v2) / len1;
        /*			diff[0] = connect[0][1]/(part_par[0].connect_length[1]-1) - connect[2][0]/(part_par[2].connect_length[0]-1); */
        v1 = 0;
        v2 = 0;
        len2 = 0;
        for (l = 0; l < (1 << i); l++) {
          help = start + ((1 << i) - 1) * parts + (1 << i);
          v1 += shift[help + l] * (part_par[help + l].connect_length[1] - 1)
              + connect[help + l][1];
          v2 += shift[help + l + parts]
              * (part_par[help + l + parts].connect_length[0] - 1)
              + connect[help + l + parts][0];
          len2 += (part_par[help + l].connect_length[1] - 1);
        }
        if (len2 < 1)
          len2 = -1;
        diff2 = (v1 - v2) / len2;
        /*			diff[1] = connect[1][1]/(part_par[1].connect_length[1]-1) - connect[3][0]/(part_par[3].connect_length[0]-1); */
        v1 = 0;
        v2 = 0;
        len3 = 0;
        for (l = 0; l < (1 << i); l++) {
          help = start + ((1 << i) - 1);
          v1 += shift[help + l * parts]
              * (part_par[help + l * parts].connect_length[3] - 1)
              + connect[help + l * parts][3];
          v2 += shift[help + 1 + l * parts]
              * (part_par[help + 1 + l * parts].connect_length[2] - 1)
              + connect[help + 1 + l * parts][2];
          len3 += (part_par[help + l * parts].connect_length[3] - 1);
        }
        if (len3 < 1)
          len3 = -1;
        diff3 = (v1 - v2) / len3;
        /*			diff[2] = connect[0][3]/(part_par[0].connect_length[3]-1) - connect[1][2]/(part_par[1].connect_length[2]-1); */
        v1 = 0;
        v2 = 0;
        len4 = 0;
        for (l = 0; l < (1 << i); l++) {
          help = start + (1 << i) * (parts + 1) - 1;
          v1 += shift[help + l * parts]
              * (part_par[help + l * parts].connect_length[3] - 1)
              + connect[help + l * parts][3];
          v2 += shift[help + 1 + l * parts]
              * (part_par[help + 1 + l * parts].connect_length[2] - 1)
              + connect[help + 1 + l * parts][2];
          len4 += (part_par[help + l * parts].connect_length[3] - 1);
        }
        if (len4 < 1)
          len4 = -1;
        diff4 = (v1 - v2) / len4;
        /*			diff[3] = connect[2][3]/(part_par[2].connect_length[3]-1) - connect[3][2]/(part_par[3].connect_length[2]-1); */

        if (len1 == -1) {
          shift2 = diff3;
          shift3 = shift2 + diff2;
          shift1 = shift3 - diff4;
        } else if (len2 == -1) {
          shift2 = diff3;
          shift1 = diff1;
          shift3 = shift1 + diff4;
        } else if (len3 == -1) {
          shift1 = diff1;
          shift3 = shift1 + diff4;
          shift2 = shift3 - diff2;
        } else if (len4 == -1) {
          shift1 = diff1;
          shift2 = diff3;
          shift3 = shift2 + diff2;
        } else {
          error = (diff1 + diff4) - (diff3 + diff2);
          diff1 = diff1 - error / 4;
          diff2 = diff2 + error / 4;
          diff3 = diff3 + error / 4;
          shift1 = diff1;
          shift2 = diff3;
          shift3 = shift2 + diff2;
        }
        for (l = 0; l < (1 << i); l++)
          for (m = 0; m < (1 << i); m++)
            shift[start + (1 << i) + m + l * parts] += shift2;
        for (l = 0; l < (1 << i); l++)
          for (m = 0; m < (1 << i); m++)
            shift[start + (1 << i) * parts + m + l * parts] += shift1;
        for (l = 0; l < (1 << i); l++)
          for (m = 0; m < (1 << i); m++)
            shift[start + (1 << i) * (parts + 1) + m + l * parts] += shift3;
      }
  }

  for (i = 0; i < parts * parts; i++)
    for (j = 0; j < part_par[i].numofresults; j++)
      result_p[i][j] += shift[i];

  /* write results */
  for (i = 0; i < numofresults; i++)
    result[i] = 0;

  for (i = 0; i < parts * parts; i++)
    for (j = 0; j < part_par[i].numofresults; j++) {
      result[part_par[i].actconnect[j]] += result_p[i][j];
    }

  error = 0;
  for (i = 0; i < numofresults; i++) {
    result[i] /= actweight[i];
    error += result[i];
  }
  error /= numofresults;
  for (i = 0; i < numofresults; i++)
    result[i] -= error;

  /* integral control */
  for (i = 0; i < numofresults; i++)
    result_vec[i] += gain * result[i];

  return 0;
}
