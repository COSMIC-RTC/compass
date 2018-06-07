// kp_kalman_core_sparse_GPU.h
// Kalman controller (core)

#ifndef __LAM__KP_KALMAN_CORE_SPARSE_GPU__H__
#define __LAM__KP_KALMAN_CORE_SPARSE_GPU__H__

#include "kp_kalman_core_sparse.h"

#include "kp_cu_matrix.h"
#include "kp_cu_smatrix.h"
#include "kp_cu_vector.h"

class kp_kalman_core_sparse_GPU : public kp_kalman_core_sparse {
 public:
  kp_kalman_core_sparse_GPU(const kp_smatrix<KFPP>& D_Mo_,
                            const kp_smatrix<KFPP>& N_Act_,
                            const kp_smatrix<KFPP>& PROJ_, bool isZonal_,
                            cublasHandle_t cublasHandle_,
                            cusparseHandle_t cusparseHandle_);

  virtual void calculate_gain(float bruit_pix, float k_W,
                              const kp_matrix<double>& SigmaV,
                              const kp_vector<double>& atur_,
                              const kp_vector<double>& btur_);

  virtual void next_step(const kp_vector<KFPP>& Y_k, kp_vector<KFPP>& U_k);

  virtual ~kp_kalman_core_sparse_GPU();

 private:
  kp_cu_vector<KFPP> cu_U_km2;
  kp_cu_vector<KFPP> cu_U_km1;
  kp_cu_vector<KFPP> cu_X_kskm1;

  kp_cu_smatrix<KFPP> cu_D_Mo;
  kp_cu_smatrix<KFPP> cu_N_Act;
  kp_cu_smatrix<KFPP> cu_PROJ;

  kp_cu_vector<KFPP> cu_atur;
  kp_cu_vector<KFPP> cu_btur;

  kp_cu_matrix<KFPP> cu_H_inf;

  cublasHandle_t cublasHandle;
  cusparseHandle_t cusparseHandle;

  // variables de next_step
  kp_cu_vector<KFPP> cu_Y_k;
  kp_cu_vector<KFPP> cu_U_k;
  kp_cu_vector<KFPP> cu_Nact_Ukm2;
  kp_cu_vector<KFPP> cu_tmp_vec1;
  kp_cu_vector<KFPP> cu_innovation;
  kp_cu_vector<KFPP> cu_X_kskm1_tmp;
  kp_cu_vector<KFPP> cu_X_kp1sk;
  kp_cu_vector<KFPP> cu_X_kp1sk_tmp;
  kp_cu_vector<KFPP> cu_Y_kskm1;
  kp_cu_vector<KFPP> cu_A1_00_Xkdebut;
  kp_cu_vector<KFPP> cu_A1_01_Xkfin;
  kp_cu_vector<KFPP> cu_X_kp1sk_debut;
};

#endif
