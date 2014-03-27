//kp_init_from_matio.h
//function for make initializations of vector, matrises and sparce matrises
//using matio library (which read MAT matlab files)

#ifndef __KP_INIT_FROM_MATIO_H__
#define __KP_INIT_FROM_MATIO_H__


#include "kp_smatrix.h"
#include "kp_matrix.h"
#include <string>
#include <matio.h>
using namespace std;

void   kp_init4matio_smatrix(kp_smatrix& M, mat_t* mat, string name);
void   kp_init4matio_matrix(kp_matrix&M, mat_t* mat, string name);
void   kp_init4matio_vector(kp_vector&M, mat_t* mat, string name);
void   kp_init4matio_vector(vector<double>&v, mat_t* mat, string name);
void   kp_init4matio_vector(vector<int>&vi,   mat_t* mat, string name);
double kp_init4matio_double(mat_t* mat, string name);
int    kp_init4matio_int(mat_t* mat, string name);

matvar_t * kp_matio_readvar(mat_t* mat, string name);


//in matlab numeration starts from 1
//so we just supprime 1;
void   kp_init4matio_index_vector(vector<int>&vi,   mat_t* mat, string name);


#endif

