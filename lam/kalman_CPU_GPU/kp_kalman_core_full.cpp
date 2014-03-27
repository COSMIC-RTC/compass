//kp_kalman_core_full.cpp

#include "kp_kalman_core_full.h"


kp_kalman_core_full::kp_kalman_core_full(const kp_matrix& D_Mo_,
			       const kp_matrix& N_Act_,
			       const kp_matrix& PROJ_,
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
	cerr<<"Error | kp_kalman_core_full::kp_kalman_core_full | dimension problem N_Act.dim1 or D_Mo.dim2 or PROJ.dim2"<<endl;
	exit(EXIT_FAILURE);
     }

   if ((N_Act_.dim2 != PROJ_.dim1) || N_Act_.dim2 != nb_act)
     {
	cerr<<"Error | kp_kalman_core_full::kp_kalman_core_full | dimension problem N_Act.dim2 or PROJ.dim1"<<endl;
	exit(EXIT_FAILURE);
     }
} 
