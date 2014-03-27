//kp_matio_writer.cpp

#include "kp_matio_writer.h"

kp_matio_writer::kp_matio_writer(string fn)
{
   matfp = Mat_CreateVer(fn.c_str(), NULL, MAT_FT_MAT73);
   if ( NULL == matfp )
     {
	cerr<<"Error | kp_matio_writer::kp_matio_writer | Error creating MAT file"<<endl;
	exit(EXIT_FAILURE);
     }
}
//                                                                        
kp_matio_writer::~kp_matio_writer()
{
   Mat_Close(matfp);
}
//                                                                        
void kp_matio_writer::add(string name, const kp_vector& v)
{
   size_t    dims[2];
   dims[0] = v.size();
   dims[1] = 1;
   matvar_t* matvar = Mat_VarCreate(name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, v.d,0);
   if ( NULL == matvar )
     {
	cerr<<"Error | kp_matio_writer::add | Error creating variable"<<endl;
	exit(EXIT_FAILURE);
     }
   Mat_VarWrite(matfp,matvar, MAT_COMPRESSION_NONE);
   Mat_VarFree(matvar);
}
//                                                                        
void kp_matio_writer::add(string name, const kp_matrix& m)
{
   size_t    dims[2];
   dims[0] = m.dim1;
   dims[1] = m.dim2;
   matvar_t* matvar = Mat_VarCreate(name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, m.d, 0);
   if ( NULL == matvar )
     {
	cerr<<"Error | kp_matio_writer::add | Error creating variable"<<endl;
	exit(EXIT_FAILURE);
     }
   Mat_VarWrite(matfp,matvar, MAT_COMPRESSION_NONE);
   Mat_VarFree(matvar);   
}
//                                                                        
