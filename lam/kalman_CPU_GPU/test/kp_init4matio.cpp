// kp_init_from_matio.cpp

#include "kp_init4matio.h"
#include <math.h>
#include <string.h>
#include <iostream>
#include "kp_vector.h"

//
void kp_init4matio_vector(vector<double>& v, mat_t* mat, string name) {
  kp_vector<float> kpv;
  kp_init4matio_vector(kpv, mat, name);
  v.resize(kpv.size());
  for (size_t i = 0; i < v.size(); i++) v[i] = kpv[i];
}
//
void kp_init4matio_vector(vector<int>& v, mat_t* mat, string name) {
  kp_vector<float> kpv;
  kp_init4matio_vector(kpv, mat, name);
  v.resize(kpv.size());
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = kpv[i];
    if (fabs((double)v[i] - kpv[i]) > 1e-10) {
      cerr << "Error | kp_init4matio_vector{int} | " << kpv[i]
           << "not look like integer" << endl;
      exit(EXIT_FAILURE);
    }
  }
}
//
double kp_init4matio_double(mat_t* mat, string name) {
  double rez;
  matvar_t* var = kp_matio_readvar(mat, name);
  if (var->class_type != MAT_C_DOUBLE || var->data_type != MAT_T_DOUBLE ||
      var->rank != 2 || var->isComplex || var->data_size != sizeof(double) ||
      var->dims[1] != 1 || var->dims[1] != 1) {
    cerr << "Error | kp_init4matio_vector | variable " << name
         << " is not double value" << endl;
    exit(EXIT_FAILURE);
  }

  rez = ((double*)var->data)[0];
  Mat_VarFree(var);
  return rez;
}
//
int kp_init4matio_int(mat_t* mat, string name) {
  double rezd = kp_init4matio_double(mat, name);
  int rezi = (int)rint(rezd);
  if (fabs((double)rezi - rezd) > 1e-10) {
    cerr << "Error | kp_init4matio_int | " << rezd << "not look like integer"
         << endl;
    exit(EXIT_FAILURE);
  }
  return rezi;
}
//
matvar_t* kp_matio_readvar(mat_t* mat, string name) {
  matvar_t* var = Mat_VarRead(mat, name.c_str());
  if (var == NULL) {
    Mat_Rewind(mat);
    var = Mat_VarRead(mat, name.c_str());
  }
  if (var == NULL) {
    cerr << "error | kp_matio_readvar | cannot read " << name << endl;
    exit(EXIT_FAILURE);
  }
  return var;
}
//
void kp_init4matio_index_vector(vector<int>& v, mat_t* mat, string name) {
  kp_init4matio_vector(v, mat, name);
  for (size_t i = 0; i < v.size(); i++) {
    v[i] -= 1;
    if (v[i] < 0) {
      cerr << "Error | kp_init4matio_index_vector | not a index " << v[i]
           << endl;
      exit(EXIT_FAILURE);
    }
  }
}
