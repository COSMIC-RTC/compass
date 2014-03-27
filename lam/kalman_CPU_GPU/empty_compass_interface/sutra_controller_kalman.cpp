//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

//nslope = dim1(D_Mo)
//nactu  = dim1(N_Act)
sutra_controller_kalman::sutra_controller_kalman(carma_context* context,
						 carma_obj<float>& cD_Mo,
						 carma_obj<float>& cN_Act,
						 carma_obj<float>& cPROJ,
						 bool is_zonal)
:sutra_controller(context, 0, 0)
{
}
//                                                                                     
void sutra_controller_kalman::calculate_gain(double bruit,
					     double k_W,
					     carma_obj<float>& cSigmaV,
					     carma_obj<float>& catur,
					     carma_obj<float>& cbtur)
{
}
//                                                                                     
int sutra_controller_kalman::comp_com()
{
   return -378;
}
//                                                                                     
