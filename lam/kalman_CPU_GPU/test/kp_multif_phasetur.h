//kp_multif_phasetur.h
//class for read multiple turbulent phase files

#ifndef __SEGER__KP_MULTIF_PHASETUR_H__
#define __SEGER__KP_MULTIF_PHASETUR_H__

#include <vector>
#include <string>
#include "kp_matrix.h"
using namespace std;

class kp_multif_phasetur
{
 public:
   void init(const vector<string>& files);   
   bool set_next_phase(kp_vector& CPT_km1);
   
   int get_boucle() {return boucle;};
private:
   bool read_next();
 private:
   vector<string> files;
   int fidx_next;
   kp_matrix PHAS_TUR;
   int boucle;
   int boucle_local;
};


#endif
