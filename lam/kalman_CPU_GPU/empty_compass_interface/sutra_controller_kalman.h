//sutra_controller_kalman.h

#ifndef __LAM__SUTRA_CONTROLLER_KALMAN_H__
#define __LAM__SUTRA_CONTROLLER_KALMAN_H__

#include "sutra_controller.h"


class sutra_controller_kalman:public sutra_controller
{
 public:
   sutra_controller_kalman(carma_context* context, 
			   carma_obj<float>& D_Mo,
			   carma_obj<float>& N_Act,
			   carma_obj<float>& PROJ,
			   bool is_zonal);
   
   void calculate_gain(double bruit,
		       double k_W,
		       carma_obj<float>&  SigmaV,
		       carma_obj<float>& atur,
		       carma_obj<float>& btur);
   
   virtual string get_type(){return "kalman";};
   
   virtual int comp_com();   
};


#endif
