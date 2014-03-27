//kp_cu_vector.h
//Kalman project vector

#ifndef __SEGER__KP_CU_VECTOR_H__
#define __SEGER__KP_CU_VECTOR_H__

#include "kp_vector.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "kernels.h"
#include "kp_cu_cuda.h"


using namespace std;

class kp_cu_matrix;
class kp_vector;
class kp_matrix;


class kp_cu_vector
{
 public:
	kp_cu_vector(int s_)  {_create(s_);};
	kp_cu_vector()        {_create(0);};
	kp_cu_vector(const kp_cu_vector&);
	kp_cu_vector(const kp_vector&);
	kp_cu_vector(const int vec_size, real val);
	void zeros();

	void operator=( const kp_vector& v);
	void operator=( const kp_cu_vector& v);

	void operator+=(const kp_cu_vector& v);
	void operator-=(const kp_cu_vector& v);
	void operator*=(const kp_cu_vector& v);
	void operator/=(const kp_cu_vector& v);

	void operator+=(real val);
	void operator-=(real val);
	void operator*=(real val);
	void operator/=(real val);

	

      void resize(int s_);
      int size() const {return s;};
      void init_from_matrix_column(const kp_cu_matrix& M, int col);
      void init_from_matrix_column(const kp_matrix& M, int col);

   //MATLAB code this(ind_begin + 1, int_end)
   void init_from_vector(const kp_cu_vector&, int ind_begin, int int_end);
   void init_from_vector(const kp_vector&, int ind_begin, int int_end);

   //in idx we have indexes of v which we should use
   //MATLAB code: this = v(idx)
   void init_from_idx(const kp_vector& v, const vector<int>& idx);
   void init_from_idx(const kp_cu_vector& v, const vector<int>& idx);

   //subv is vector which was created by init_from_rowidx
   //MATLAB code: this(idx) = subv
   void set_from_subvector(const kp_vector& subv, const vector<int>& idx);
   void set_from_subvector(const kp_cu_vector& subv, const vector<int>& idx);


 private:
	void _create(int s_){s = s_ ; if(s>0) kp_cu_cudaMalloc((void**)&d_cu, sizeof(real)*s);else d_cu=NULL;};
       	void _clear();



	
 public:
	 real* d_cu;
	 int s;

};

void kp_cu_inverse(kp_cu_vector& v);
void sqrt(kp_cu_vector& v);


#endif
