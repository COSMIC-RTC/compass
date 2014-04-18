/*
 * carmasparsehostobj.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: ???
 */
#include "carma_sparse_host_obj.h"

template<class T_data>
carma_sparse_host_obj<T_data>::carma_sparse_host_obj() {
  _create(0, 0, 0);
}

template<class T_data>
carma_sparse_host_obj<T_data>::~carma_sparse_host_obj() {
  _clear();
}

//
template<class T_data>
carma_sparse_host_obj<T_data>::carma_sparse_host_obj(
   const carma_sparse_host_obj<T_data>& M) {
  _create(0, 0, 0);
  this->operator=(*M);
}
//
//
template<class T_data>
void carma_sparse_host_obj<T_data>::resize(int nnz_, int dim1_, int dim2_) {
  if (nz_elem != nnz_) {
    _clear();
    _create(nnz_, dim1_, dim2_);
  } else {
    dims_data[0]=2;
    dims_data[1]=dim1_;
    dims_data[2]=dim2_;
    majorDim = 'U';
  }
}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::operator=(
    const carma_sparse_host_obj<T_data>& M) {
  resize(M->nz_elem, M->getData(1), M->getData(2));
  memcpy(h_data, M->h_data, sizeof(T_data) * nz_elem);
  memcpy(rowind, M->rowind, sizeof(int) * nz_elem);
  memcpy(colind, M->colind, sizeof(int) * nz_elem);

  majorDim = M->majorDim;

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
    if ((*idx)[i] < 0 || (*idx)[i] >= M->dim1) {
      cerr
          << "Error | carma_sparse_host_obj<T_data>::init_from_rowidx | index error"
          << endl;
      cout << (*idx)[i] << " " << M->dim1 << endl;
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
  vector<pair<int, pair<int, T_data> > > sM(M->nnz);
  //sM[i].first         -> rowind
  //sM[i].second.first  -> colind
  //sM[i].second.second -> value
  for (int i = 0; i < M->nnz; i++) {
    sM[i].first = M->rowind[i];
    sM[i].second.first = M->colind[i];
    sM[i].second.second = M->values[i];
  }
  sort(sM.begin(), sM.end());

  //make first loop to calculate number of non zero elements
  size_t i_idx = 0;
  int i_this = 0;
  //idx in zero-based indexing but matrix in one-based indexing
  for (int i = 0; i < M->nnz; i++) {
    while (i_idx < (*idx).size() && ((*idx)[i_idx] + 1) < sM[i].first)
      i_idx++;
    if (i_idx < (*idx).size() && sM[i].first == (*idx)[i_idx] + 1) {
      i_this++;
    }
  }
  resize(i_this, (*idx).size(), M->dim2);

  //make second loop to initialize non zero elements

  i_idx = 0;
  i_this = 0;
  for (int i = 0; i < M->nnz; i++) {
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
  resize(M->nnz, M->dim1, M->dim2);
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
void carma_sparse_host_obj<T_data>::init_from_matrix(carma_host_obj<T_data>* B) {
  int nb_el = B.dim1 * B.dim2;
  int* row = new int[nb_el];
  int* col = new int[nb_el];
  T_data* val = new T_data[nb_el];
  int cmpt = 0;
  for (int j = 0; j < B.dim2; j++) {
    for (int i = 0; i < B.dim1; i++) {
      if (B(i, j) != 0) {
        row[cmpt] = i + 1;
        col[cmpt] = j + 1;
        val[cmpt] = B(i, j);
        cmpt++;
      }
    }
  }
  resize(cmpt, B.dim1, B.dim2);
  for (int i = 0; i < nz_elem; i++) {
    rowind[i] = row[i];
    colind[i] = col[i];
    h_data[i] = val[i];
  }

  delete[] row;
  delete[] col;
  delete[] val;
}

template<class T_data>
void carma_sparse_host_obj<T_data>::check() {
  cerr << "CHECK smatrix (DEBUG)" << endl;
  for (int i = 0; i < nz_elem; i++) {
    if (rowind[i] < 1 || rowind[i] > dims_data[1] || colind[i] < 1
        || colind[i] > dims_data[2]) {
      cerr << "Error | carma_sparse_host_obj<T_data>::check | bad index"
          << endl;
      cerr << dims_data[1] << "x" << dims_data[2] << " " << rowind[i] << "x" << colind[i]
          << endl;
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
        char *matdescra, T_data *val, long *rowind, long *colind,
        long *nnz, T_data *x, T_data *beta, T_data *y)) {
  if (A->dim2 != x->getDims(1) || A->dim1 != y->getDims(1)) {
    cerr << "Error | kp_cscmv | diminsion problem" << endl;
    throw "Error | kp_cscmv | diminsion problem";
    //exit(EXIT_FAILURE);
  }
//   A.check();
  ptr_coomv("N", &(A->dim1), &(A->dim2), &alpha, "GFFFFFFFF", A->values,
      A->rowind, A->colind, &(A->nnz), x, &betta, y.d);
//   mkl_dcscmv("N", &A.dim1, &A.dim2, &alpha, "GFFFFFFF",
//        A.values, A.rows, A.pointerB, A.pointerB + 1,
//        x.d,  &betta, y.d);

}
//
/*void kp_gemm(T_data alpha, carma_sparse_host_obj<T_data>* A, carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C)
 {
 if (C.dim1 != A.dim1 || A.dim2 != B.dim1 || C.dim2 != B.dim2)
 {
 cerr<<"Error | kp_cscmv | diminsion problem"<<endl;
 exit(EXIT_FAILURE);
 }
 //   A.check();

 //we make matrix multiplication using vector multiplication because for some reason is
 //better parallised
 #pragma omp parallel for
 for (int i = 0 ; i < B.dim2 ; i++) //for each column of matrix B
 {
 #ifdef ST_SINGLE
 mkl_scoomv("N", &A.dim1, &A.dim2, &alpha,"GFFFFFFF",
 A.values, A.rowind, A.colind, &A.nnz,
 B.d + i * B.dim1, &betta, C.d + i * C.dim1);
 #else
 mkl_dcoomv("N", &A.dim1, &A.dim2, &alpha,"GFFFFFFF",
 A.values, A.rowind, A.colind, &A.nnz,
 B.d + i * B.dim1, &betta, C.d + i * C.dim1);
 #endif
 }
 }*/
//
template<class T_data>
void kp_gemm(char op_A, T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C) {
  int opA_dim1, opA_dim2;
  kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);

  if (C.dim1 != opA_dim1 || opA_dim2 != B.dim1 || C.dim2 != B.dim2) {
    cerr << "Error | kp_gemm (sparce) | diminsion problem" << endl;
    throw "Error | kp_gemm (sparce) | diminsion problem";
    //exit(EXIT_FAILURE);
  }
//   A.check();
#ifdef KP_SINGLE
  mkl_scoomm(&op_A, &A.dim1, &C.dim2, &A.dim2,&alpha,"GFFFFFFF",
      A.values, A.rowind, A.colind, &A.nnz,
      B.d, &B.dim1, &betta, C.d, & C.dim1);
#else
  mkl_dcoomm(&op_A, &A.dim1, &C.dim2, &A.dim2, &alpha, "GFFFFFFF", A.values,
      A.rowind, A.colind, &A.nnz, B.d, &B.dim1, &betta, C.d, &C.dim1);
//   mkl_dcscmm("N", &A.dim1, &C.dim2,  &A.dim2, &alpha, "GCCCCCCC",
//        A.values, A.rows, A.pointerB, A.pointerB+1,
//        B.d, &B.tda,  &betta, C.d, & C.tda);
#endif

#ifdef __WITH_FLOPS_CALCULATOR__
  kp_fpc_global.add_sparce(2 * (long)A.nnz * (long)B.dim2);
#endif
}
//
template<class T_data>
void kp_check_op_set_dim(int op, const carma_sparse_host_obj<T_data>*M,
    int* dim1, int* dim2) {
  if (op == 'N') {
    dim1 = M->dim1;
    dim2 = M->dim2;
  } else if (op == 'T') {
    dim1 = M->dim2;
    dim2 = M->dim1;
  } else {
    cerr << "Error | kp_check_op_set_dim | op should ge either N either T"
        << endl;
    throw "Error | kp_check_op_set_dim | op should ge either N either T";
    //exit(EXIT_FAILURE);
  }
}
//
template<class T_data>
void carma_sparse_host_obj<T_data>::resize2rowMajor() {
  int i, j;
  int indice;
  vector< vector< pair<int, int> > > position(dims_data[1], vector< pair<int, int> >());
  int* rowind2 = new int[nz_elem];
  int* colind2 = new int[nz_elem];
  T_data* values2 = new T_data[nz_elem];

  //cout<<"dim1="<<this->dim1<<" dim2="<<this->dim2<<" nnz="<<this->nnz<<endl;
  //cout<<"i="<<172248<<" rowind[i]="<<rowind[172248]<<" colind[i]="<<colind[172248]<<endl;

  for (i = 0; i < nz_elem; i++) {
    //cout<<"i="<<i<<" rowind[i]="<<rowind[i]-1<<" colind[i]="<<colind[i]-1<<endl;
    position[rowind[i] - 1].push_back(make_pair(colind[i], i));

  }

  indice = 0;
  // FIXME: problem with sort
//  for (i = 0; i < dim1; i++) {
//    if (position[i].size() > 0) {
//      std::sort<vector< pair<int, int> >::const_iterator>(position[i].begin(), position[i].end());
//      for (int j = 0; j < position[i].size(); j++) {
//        rowind2[indice] = rowind[(position[i][j]).second];
//        colind2[indice] = colind[(position[i][j]).second];
//        values2[indice] = values[(position[i][j]).second];
//        indice++;
//      }
//    }
//  }
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
  int i, j;
  int indice;
  vector< vector< pair<int, int> > > position(dims_data[2], vector<pair<int, int> >());
  int* rowind2 = new int[nz_elem];
  int* colind2 = new int[nz_elem];
  T_data* values2 = new T_data[nz_elem];

  for (i = 0; i < nz_elem; i++) {
    position[colind[i] - 1].push_back(make_pair(rowind[i], i));
  }

  indice = 0;
  // FIXME: problem with sort
//  for (i = 0; i < dim2; i++) {
//    if (position[i].size() > 0) {
//      sort(position[i].begin(), position[i].end());
//      for (int j = 0; j < position[i].size(); j++) {
//        colind2[indice] = colind[(position[i][j]).second];
//        rowind2[indice] = rowind[(position[i][j]).second];
//        values2[indice] = values[(position[i][j]).second];
//        indice++;
//      }
//    }
//  }
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
