//kp_multif_phasetur.cpp

#include "kp_multif_phasetur.h"
#include "kp_init4matio.h"
#include <matio.h>
#include <iostream>

void kp_multif_phasetur::init(const vector<string>& files_)
{
   files     = files_;
   if (files.size() == 0)
     {
	cerr<<"Error | kp_multif_phasetur::init | empty list of phas_tur files"<<endl;
	exit(EXIT_FAILURE);
     }
   fidx_next = 0;
   boucle    = 0;
   boucle_local = 0;
   read_next();
}
//                                                                                            
bool kp_multif_phasetur::set_next_phase(kp_vector& CPT_km1)
{
   if (boucle_local >= PHAS_TUR.dim2)
     {
	if (!read_next())
	  return false;
     }
   CPT_km1.init_from_matrix_column(PHAS_TUR, boucle_local);
   boucle_local++;
   boucle++;
   return true;
}
//                                                                                            
bool kp_multif_phasetur::read_next()
{
   if (fidx_next >= files.size())
     return false;
   string f = files[fidx_next++];
   mat_t* mat = Mat_Open(f.c_str(), MAT_ACC_RDONLY);
   if (!mat)
     {
	cerr<<"Error | kp_multif_phasetur::read_next | cannot open "<<f<<endl;
	exit(EXIT_FAILURE);
     }
   kp_init4matio_matrix(PHAS_TUR, mat, "PHAS_TUR");
   if (PHAS_TUR.dim2 <= 0)
     {
	cerr<<"Error | kp_multif_phasetur::read_next | in "<<f<<" we have empty PHAS_TUR"<<endl;
	exit(EXIT_FAILURE);
     }
   boucle_local = 0;
   Mat_Close(mat);
   return true;
}
//                                                                                            
