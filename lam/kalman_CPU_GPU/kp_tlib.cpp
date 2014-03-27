//kp_smatrix2.cpp

#include "kp_tlib.h"
#include <iostream>


// indr1, indc1, indr2, indc2 are ZERO-BASED
// ne peut etre utilise que lorsque les indices varient par colonne PUIS par ligne (column major)
void init_smatrix_from_larger_smatrix(kp_smatrix& T, const kp_smatrix&M, const int indr1, const int indc1, const int indr2, const int indc2)
{
   if (&T == &M)
     {
	cerr<<"Error | init_smatrix_from_larger_smatrix | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }

   if( (indr1>indr2) || (indc1>indc2) || (indr2>M.dim1) || (indc2>M.dim2) )
     {	
	cerr<<"Error | kp_smatrix::init_from_larger_matrix | index are not correct"<<endl;
	exit(EXIT_FAILURE);
     }
   
   int ind_nnz=0;
   int* rowind_tmp = (int*)malloc(M.nnz*sizeof(int));
   int* colind_tmp = (int*)malloc(M.nnz*sizeof(int));
   real* values_tmp = (real*)malloc(M.nnz*sizeof(real));
   for (int i= 0 ; i<M.nnz ; i++)
     {
	if ( M.colind[i]>(indc2+1) ) break;
	else if( (M.rowind[i]>indr1) && (M.rowind[i]<=(indr2+1)) && (M.colind[i]>indc1) )
	  {
	     rowind_tmp[ind_nnz]=M.rowind[i];
	     colind_tmp[ind_nnz]=M.colind[i];
	     values_tmp[ind_nnz]=M.values[i];
	     ind_nnz++;
	  }
	
     }
   T.resize(ind_nnz, indr2-indr1+1, indc2-indc1+1);
   int rowind_im1 = -1;
   int colind_im1 = -1;
   for (int i=0 ; i<ind_nnz ; i++)
     {
	T.rowind[i] = rowind_tmp[i];
	T.colind[i] = colind_tmp[i];
	T.values[i] = values_tmp[i];

	if( (colind_im1>colind_tmp[i]) || (colind_im1==colind_tmp[i] && rowind_im1>rowind_tmp[i]) )
	{
		cerr<<"Error | init_smatrix_from_larger_smatrix | the larger sparse matrix must be column major"<<endl;
		exit(EXIT_FAILURE);
		break;
	}

	rowind_im1 = rowind_tmp[i];
	colind_im1 = colind_tmp[i];
     }
   free(rowind_tmp);
   free(colind_tmp);
   free(values_tmp);
}
