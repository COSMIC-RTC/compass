/*
 * carmasparsehostobj.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: ???
 */
#include "carma_sparse_host_obj.h"
#include <algorithm>

template<class T_data>
carma_sparse_host_obj<T_data>::carma_sparse_host_obj() {
  _create(0, 0, 0);
}

template<class T_data>
carma_sparse_host_obj<T_data>::~carma_sparse_host_obj() {
  _clear();
}

template<class T_data>
carma_sparse_host_obj<T_data>::carma_sparse_host_obj(
    carma_sparse_host_obj<T_data>& M) {
  _create(0, 0, 0);
  this->operator=(M);
}

template<class T_data>
void carma_sparse_host_obj<T_data>::resize(int nnz_, int dim1_, int dim2_) {
  if (nz_elem != nnz_) {
    _clear();
    _create(nnz_, dim1_, dim2_);
  } else {
    dims_data[0] = 2;
    dims_data[1] = dim1_;
    dims_data[2] = dim2_;
    majorDim = 'U';
  }
}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::operator=(
    carma_sparse_host_obj<T_data>& M) {
  resize(M.nz_elem, M.getDims(1), M.getDims(2));
  memcpy(h_data, M.h_data, sizeof(T_data) * nz_elem);
  memcpy(rowind, M.rowind, sizeof(int) * nz_elem);
  memcpy(colind, M.colind, sizeof(int) * nz_elem);

  majorDim = M.majorDim;

}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::init_from_rowidx(
    carma_sparse_host_obj<T_data>* M, vector<int>* idx) {
  if (this == M) {
    cerr
        << "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | the same matrix"
        << endl;
    throw "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | the same matrix";
    //exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < (*idx).size(); i++)
    if ((*idx)[i] < 0 || (*idx)[i] >= M->getDims(1)) {
      cerr
          << "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | index error"
          << endl;
      cout << (*idx)[i] << " " << M->getDims(1) << endl;
      throw "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | index error";
      //exit(EXIT_FAILURE);
    }

  for (size_t i = 1; i < (*idx).size(); i++)
    if ((*idx)[i] <= (*idx)[i - 1]) {
      cerr
          << "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | idx haven't been sorted"
          << endl;
      throw "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | idx haven't been sorted";
      //exit(EXIT_FAILURE);
    }

  //sort rows
  vector<pair<int, pair<int, T_data> > > sM(M->nz_elem);
  //sM[i].first         -> rowind
  //sM[i].second.first  -> colind
  //sM[i].second.second -> value
  for (int i = 0; i < M->nz_elem; i++) {
    sM[i].first = M->rowind[i];
    sM[i].second.first = M->colind[i];
    sM[i].second.second = M->h_data[i];
  }
  sort(sM.begin(), sM.end());

  //make first loop to calculate number of non zero elements
  size_t i_idx = 0;
  int i_this = 0;
  //idx in zero-based indexing but matrix in one-based indexing
  for (int i = 0; i < M->nz_elem; i++) {
    while (i_idx < (*idx).size() && ((*idx)[i_idx] + 1) < sM[i].first)
      i_idx++;
    if (i_idx < (*idx).size() && sM[i].first == (*idx)[i_idx] + 1) {
      i_this++;
    }
  }
  resize(i_this, (*idx).size(), M->getDims(2));

  //make second loop to initialize non zero elements

  i_idx = 0;
  i_this = 0;
  for (int i = 0; i < M->nz_elem; i++) {
    while (i_idx < (*idx).size() && (*idx)[i_idx] + 1 < sM[i].first)
      i_idx++;
    if (i_idx < (*idx).size() && sM[i].first == (*idx)[i_idx] + 1) {
      rowind[i_this] = i_idx + 1;
      colind[i_this] = sM[i].second.first;
      h_data[i_this] = sM[i].second.second;
      i_this++;
    }
  }

  majorDim = 'R'; // *this est une matrice row-major a la fin de cette methode

}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::init_from_transpose(
    carma_sparse_host_obj<T_data>*M) {
  resize(M->nz_elem, M->getDims(1), M->getDims(2));
  memcpy(h_data, M->h_data, sizeof(T_data) * nz_elem);
  memcpy(colind, M->rowind, sizeof(int) * nz_elem);
  memcpy(rowind, M->colind, sizeof(int) * nz_elem);

  if (M->majorDim == 'C')
    majorDim = 'R';
  else if (M->majorDim == 'R')
    majorDim = 'C';
  else
    majorDim = 'U';

}
//

template<class T_data>
void carma_sparse_host_obj<T_data>::init_from_matrix(
    carma_host_obj<T_data>* B, char majorDim='R') {

  long nb_el = B->getDims(1) * B->getDims(2);
  T_data* valB = B->getData();
  long cmpt = 0;
  for (long i = 0; i < nb_el; i++) {
    if (*valB != 0) {
      cmpt++;
    }
    valB++;
  }

  resize(cmpt, B->getDims(1), B->getDims(2));
  int *dim1, *dim2, ld;
  if(majorDim=='R'){
    dim1=this->rowind;
    dim2=this->colind;
    ld=dims_data[1];
  } else {
    dim1=this->colind;
    dim2=this->rowind;
    ld=dims_data[2];
  }
  valB = B->getData();
  long index=0;
  for(long i = 0;(index<cmpt)||(i < nb_el); i++) {
    if (*valB != 0) {
      dim1[index] = i / ld;
      dim2[index] = i - dim1[index] * ld;
      h_data[index++]=*valB;
    }
    valB++;
  }

  if(majorDim=='C'){
    //sort rows
    vector<pair<int, pair<int, T_data> > > sM(nz_elem);
    //sM[i].first         -> rowind
    //sM[i].second.first  -> colind
    //sM[i].second.second -> value
    for (int i = 0; i < nz_elem; i++) {
      sM[i].first = rowind[i];
      sM[i].second.first = colind[i];
      sM[i].second.second = h_data[i];
    }
    sort(sM.begin(), sM.end());
    for (int i = 0; i < nz_elem; i++) {
      rowind[i] = sM[i].first;
      colind[i] = sM[i].second.first;
      h_data[i] = sM[i].second.second;
    }
  }

  this->majorDim = majorDim;
}

template<class T_data>
void carma_sparse_host_obj<T_data>::check() {
  cerr << "CHECK smatrix (DEBUG)" << endl;
  for (int i = 0; i < nz_elem; i++) {
    if (rowind[i] < 0 || rowind[i] > dims_data[1] || colind[i] < 0
        || colind[i] > dims_data[2]) {
      cerr << "Error | carma_sparse_host_obj<T_data>::check | bad index"
          << endl;
      cerr << dims_data[1] << "x" << dims_data[2] << " " << rowind[i] << "x"
          << colind[i] << endl;
      throw "Error | carma_sparse_host_obj<T_data>::check | bad index";
      //exit(EXIT_FAILURE);
    }
  }
}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::_create(int nnz_, int dim1_, int dim2_) {
  nz_elem = nnz_;
  dims_data[0] = 2;
  dims_data[1] = dim1_;
  dims_data[2] = dim2_;
  h_data = new T_data[nz_elem];
  rowind = new int[nz_elem];
  colind = new int[nz_elem];
  majorDim = 'U';
}

//
template<class T_data>
void carma_sparse_host_obj<T_data>::_clear() {
  if (h_data == NULL || rowind == NULL || colind == NULL) {
    cerr << "Error | carma_sparse_host_obj<T_data>::_clear | double clear"
        << endl;
    throw "Error | carma_sparse_host_obj<T_data>::_clear | double clear";
    //exit(EXIT_FAILURE);
  }
  delete[] h_data;
  delete[] rowind;
  delete[] colind;
  h_data = NULL;
  rowind = NULL;
  colind = NULL;
  nz_elem = 0;
  dims_data[0] = 2;
  dims_data[1] = 0;
  dims_data[2] = 0;
  majorDim = 'U';
}
//
template<class T_data>
void carma_gemv(T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* x, T_data betta, carma_host_obj<T_data>* y,
    void (*ptr_coomv)(char *transa, long *m, long *k, T_data *alpha,
        char *matdescra, T_data *val, int *rowind, int *colind, int *nnz,
        T_data *x, T_data *beta, T_data *y)) {
  if (A->getDims(2) != x->getDims(1) || A->getDims(1) != y->getDims(1)) {
    cerr << "Error | kp_cscmv | diminsion problem" << endl;
    throw "Error | kp_cscmv | diminsion problem";
    //exit(EXIT_FAILURE);
  }
//   A.check();
  long dim1 = A->getDims(1);
  long dim2 = A->getDims(2);
  int nz_elem = A->nz_elem;

  char transa[]="N";
  char matdescra[]="GFFFFFFFF";
  ptr_coomv(transa, &dim1, &dim2, &alpha, matdescra, A->getData(), A->rowind,
      A->colind, &nz_elem, x->getData(), &betta, y->getData());
//   mkl_dcscmv("N", &A.dim1, &A.dim2, &alpha, "GFFFFFFF",
//        A.values, A.rows, A.pointerB, A.pointerB + 1,
//        x.d,  &betta, y.d);

}

template<class T_data>
void carma_gemm(char op_A, T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C) {
//   A.check();
  /* FIXME: à coder sans MKL ou pas...
   long dimA1=A->getDims(1);
   long dimA2=A->getDims(2);
   long dimB1=B->getDims(1);
   long dimB2=B->getDims(2);
   long dimC1=C->getDims(1);
   long dimC2=C->getDims(2);
   long nz_elem=A->nz_elem;

   mkl_dcoomm(&op_A, &dimA1, &dimC2, &dimA2, &alpha, "GFFFFFFF", A->getData(),
   A->rowind, A->colind, &nz_elem, B->getData(), &dimB1, &betta, C->getData(), &dimC1);
   */

  //   mkl_dcscmm("N", &A.dim1, &C.dim2,  &A.dim2, &alpha, "GCCCCCCC",
//        A.values, A.rows, A.pointerB, A.pointerB+1,
//        B.d, &B.tda,  &betta, C.d, & C.tda);
}

template<class T_data>
void carma_sparse_host_obj<T_data>::resize2rowMajor() {
  int i;
  int indice;
  vector<vector<pair<int, int> > > position(dims_data[1],
      vector<pair<int, int> >());
  int* rowind2 = new int[nz_elem];
  int* colind2 = new int[nz_elem];
  T_data* values2 = new T_data[nz_elem];

  //cout<<"dim1="<<this->dim1<<" dim2="<<this->dim2<<" nnz="<<this->nnz<<endl;
  //cout<<"i="<<172248<<" rowind[i]="<<rowind[172248]<<" colind[i]="<<colind[172248]<<endl;

  for (i = 0; i < nz_elem; i++) {
    //cout<<"i="<<i<<" rowind[i]="<<rowind[i]-1<<" colind[i]="<<colind[i]-1<<endl;
    position[rowind[i]].push_back(make_pair(colind[i], i));

  }

  indice = 0;
  for (i = 0; i < dims_data[1]; i++) {
    if (position[i].size() > 0) {
      sort(position[i].begin(), position[i].end());
      for ( size_t j = 0; j < position[i].size(); j++) {
        rowind2[indice] = rowind[(position[i][j]).second];
        colind2[indice] = colind[(position[i][j]).second];
        values2[indice] = h_data[(position[i][j]).second];
        indice++;
      }
    }
  }
  if (indice != nz_elem) {
    cerr
        << "Erreur | carma_sparse_host_obj<T_data>::resize2rowMajor | erreur lors de la conversion."
        << endl;
    throw "Erreur | carma_sparse_host_obj<T_data>::resize2rowMajor | erreur lors de la conversion.";
    //exit(EXIT_FAILURE);
  }

  delete[] rowind;
  rowind = rowind2;
  delete[] colind;
  colind = colind2;
  delete[] h_data;
  h_data = values2;

  majorDim = 'R';

}
template<class T_data>
void carma_sparse_host_obj<T_data>::resize2colMajor() {
  int i;
  int indice;
  vector<vector<pair<int, int> > > position(dims_data[2],
      vector<pair<int, int> >());
  int* rowind2 = new int[nz_elem];
  int* colind2 = new int[nz_elem];
  T_data* values2 = new T_data[nz_elem];

  for (i = 0; i < nz_elem; i++) {
    position[colind[i]].push_back(make_pair(rowind[i], i));
  }

  indice = 0;
  for (i = 0; i < dims_data[2]; i++) {
    if (position[i].size() > 0) {
      sort(position[i].begin(), position[i].end());
      for (size_t j = 0; j < position[i].size(); j++) {
        colind2[indice] = colind[(position[i][j]).second];
        rowind2[indice] = rowind[(position[i][j]).second];
        values2[indice] = h_data[(position[i][j]).second];
        indice++;
      }
    }
  }
  if (indice != nz_elem) {
    cerr
        << "Erreur | carma_sparse_host_obj<T_data>::resize2colMajor | erreur lors de la conversion."
        << endl;
    throw "Erreur | carma_sparse_host_obj<T_data>::resize2colMajor | erreur lors de la conversion.";
    //exit(EXIT_FAILURE);
  }

  delete[] rowind;
  rowind = rowind2;
  delete[] colind;
  colind = colind2;
  delete[] h_data;
  h_data = values2;

  majorDim = 'C';

}

#define EXPLICITE_TEMPLATE(T_data) template carma_sparse_host_obj<T_data>::carma_sparse_host_obj(); \
    template carma_sparse_host_obj<T_data>::~carma_sparse_host_obj(); \
    template carma_sparse_host_obj<T_data>::carma_sparse_host_obj( \
        carma_sparse_host_obj<T_data>& M); \
    template void carma_sparse_host_obj<T_data>::resize(int nnz_, int dim1_, int dim2_); \
    template void carma_sparse_host_obj<T_data>::operator=( carma_sparse_host_obj<T_data>& M); \
    template void carma_sparse_host_obj<T_data>::init_from_rowidx( \
        carma_sparse_host_obj<T_data>* M, vector<int>* idx); \
    template void carma_sparse_host_obj<T_data>::init_from_transpose( \
        carma_sparse_host_obj<T_data>*M); \
    template void carma_sparse_host_obj<T_data>::init_from_matrix(carma_host_obj<T_data>* B, char majorDim='R'); \
    template void carma_sparse_host_obj<T_data>::check(); \
    template void carma_sparse_host_obj<T_data>::_create(int nnz_, int dim1_, int dim2_); \
    template void carma_sparse_host_obj<T_data>::_clear(); \
    template void carma_gemv(T_data alpha, carma_sparse_host_obj<T_data>* A, \
        carma_host_obj<T_data>* x, T_data betta, carma_host_obj<T_data>* y, \
        void (*ptr_coomv)(char *transa, long *m, long *k, T_data *alpha, \
            char *matdescra, T_data *val, int *rowind, int *colind, \
            int *nnz, T_data *x, T_data *beta, T_data *y)); \
    template void carma_gemm(char op_A, T_data alpha, carma_sparse_host_obj<T_data>* A, \
        carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C); \
    template void carma_sparse_host_obj<T_data>::resize2rowMajor(); \
    template void carma_sparse_host_obj<T_data>::resize2colMajor();

EXPLICITE_TEMPLATE(double)


#undef EXPLICITE_TEMPLATE
