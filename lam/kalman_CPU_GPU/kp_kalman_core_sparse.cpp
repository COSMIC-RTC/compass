//kp_kalman_core_sparse.cpp

#include "kp_kalman_core_sparse.h"



kp_kalman_core_sparse::kp_kalman_core_sparse(const kp_smatrix& D_Mo_,
		const kp_smatrix& N_Act_,
		const kp_smatrix& PROJ_,
                bool isZonal_)
:isZonal(isZonal_)
{
   nb_p   = D_Mo_.dim1;
   nb_act = PROJ_.dim1;
   if (!isZonal) 
   {
	   nb_z = N_Act_.dim1;
	   nb_az = nb_z;
   }
   else
   {
	   nb_z = -1;
	   nb_az = nb_act;
   }

   nb_n = 2*nb_az;
   gainComputed = false;

   if ((N_Act_.dim1 != D_Mo_.dim2) || N_Act_.dim1 != PROJ_.dim2 || N_Act_.dim1 != nb_az)
     {
	cerr<<"Error | kp_kalman_core_sparse::kp_kalman_core_sparse | dimension problem N_Act.dim1 or D_Mo.dim2 or PROJ.dim2"<<endl;
	exit(EXIT_FAILURE);
     }

   if ((N_Act_.dim2 != PROJ_.dim1) || N_Act_.dim2 != nb_act)
     {
	cerr<<"Error | kp_kalman_core_sparse::kp_kalman_core_sparse | dimension problem N_Act.dim2 or PROJ.dim1"<<endl;
	exit(EXIT_FAILURE);
     }
} 
