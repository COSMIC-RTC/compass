// kp_kalman_core.h
// Kalman controller (core)

#ifndef __LAM__KP_KALMAN_CORE_FULL__H__
#define __LAM__KP_KALMAN_CORE_FULL__H__

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "kp_KFPP.h"
#include "kp_cu_timer.h"

#include "kp_matrix.h"
#include "kp_vector.h"

using namespace std;

class kp_kalman_core_full {
 public:
  // make initialization and calculate gain
  kp_kalman_core_full(const kp_matrix<KFPP>& D_Mo_,
                      const kp_matrix<KFPP>& N_Act_,
                      const kp_matrix<KFPP>& PROJ_, bool isZonal_);

  virtual void calculate_gain(float bruit, float k_W,
                              const kp_matrix<double>& SigmaV,
                              const kp_vector<double>& atur_,
                              const kp_vector<double>& btur_) = 0;

  virtual void next_step(const kp_vector<KFPP>& Y_k, kp_vector<KFPP>& U_k) = 0;

  virtual ~kp_kalman_core_full() {}

 protected:
  int ordreAR;
  int nb_p, nb_act, nb_z, nb_az, nb_n;
  bool gainComputed;
  bool isZonal;

 public:
  kp_cu_timer temps_boucle;
  kp_cu_timer temps_boucle_op1;
  kp_cu_timer temps_boucle_op2;
  kp_cu_timer temps_boucle_op3;
};

#endif
