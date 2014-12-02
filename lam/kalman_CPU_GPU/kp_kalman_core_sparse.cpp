//kp_kalman_core_sparse.cpp

#include "kp_kalman_core_sparse.h"

#include "kp_vector.h"
#include "kp_matrix.h"
#include "kp_smatrix.h"


kp_kalman_core_sparse::kp_kalman_core_sparse(const kp_smatrix<KFPP>& D_Mo_,
		const kp_smatrix<KFPP>& N_Act_,
		const kp_smatrix<KFPP>& PROJ_,
                bool isZonal_)
:isZonal(isZonal_)
{
   nb_p   = D_Mo_.getDim1();
   nb_act = PROJ_.getDim1();
   if (!isZonal) 
   {
	   nb_z = N_Act_.getDim1();
	   nb_az = nb_z;
   }
   else
   {
	   nb_z = -1;
	   nb_az = nb_act;
   }

   nb_n = 2*nb_az;
   gainComputed = false;

   if ((N_Act_.getDim1() != D_Mo_.getDim2()) || N_Act_.getDim1() != PROJ_.getDim2() || N_Act_.getDim1() != nb_az)
     {
	cerr<<"Error | kp_kalman_core_sparse::kp_kalman_core_sparse | dimension problem N_Act.getDim1() or D_Mo.getDim2() or PROJ.getDim2()"<<endl;
	exit(EXIT_FAILURE);
     }

   if ((N_Act_.getDim2() != PROJ_.getDim1()) || N_Act_.getDim2() != nb_act)
     {
	cerr<<"Error | kp_kalman_core_sparse::kp_kalman_core_sparse | dimension problem N_Act.getDim2() or PROJ.getDim1()"<<endl;
	exit(EXIT_FAILURE);
     }
} 
