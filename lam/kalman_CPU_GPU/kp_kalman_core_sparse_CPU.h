// kp_kalman_core_sparse_CPU.h
// Kalman controller (core)

#ifndef __LAM__KP_KALMAN_CORE_SPARSE_CPU__H__
#define __LAM__KP_KALMAN_CORE_SPARSE_CPU__H__

#include "kp_kalman_core_sparse.h"

#include "kp_matrix.h"
#include "kp_smatrix.h"
#include "kp_vector.h"

#include "kp_randn.h"

class kp_kalman_core_sparse_CPU : public kp_kalman_core_sparse {
 public:
  kp_kalman_core_sparse_CPU(const kp_smatrix<KFPP>& D_Mo_,
                            const kp_smatrix<KFPP>& N_Act_,
                            const kp_smatrix<KFPP>& PROJ_, bool isZonal_);

  virtual void calculate_gain(float bruit_pix, float k_W,
                              const kp_matrix<double>& SigmaV,
                              const kp_vector<double>& atur_,
                              const kp_vector<double>& btur_);

  virtual void next_step(const kp_vector<KFPP>& Y_k, kp_vector<KFPP>& U_k);

  virtual ~kp_kalman_core_sparse_CPU() {}

 private:
  kp_vector<KFPP> U_km2;
  kp_vector<KFPP> U_km1;
  kp_vector<KFPP> X_kskm1;
  kp_smatrix<KFPP> D_Mo, N_Act, PROJ;
  kp_vector<KFPP> atur;
  kp_vector<KFPP> btur;
  kp_matrix<KFPP> H_inf;

  // variables de next_step
  kp_vector<KFPP> Nact_Ukm2;
  kp_vector<KFPP> tmp_vec1;
  kp_vector<KFPP> Y_kskm1;
  kp_vector<KFPP> innovation;
  kp_vector<KFPP> X_kskm1_tmp;
  kp_vector<KFPP> A1_00_Xkdebut;
  kp_vector<KFPP> A1_01_Xkfin;
  kp_vector<KFPP> X_kp1sk_debut;
  kp_vector<KFPP> X_kp1sk_tmp;
  kp_vector<KFPP> X_kp1sk;
};

#endif
