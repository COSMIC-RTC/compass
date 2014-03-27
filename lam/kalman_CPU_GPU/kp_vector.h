//kp_vector.h
//Kalman project vector

#ifndef __SEGER__KP_VECTOR_H__
#define __SEGER__KP_VECTOR_H__

#include <stdlib.h>
#include <vector>
#include <iostream>
#include "kp_real.h"
using namespace std;

class kp_matrix;


class kp_vector
{
 public:
   kp_vector(int s_)      {_create(s_);};
   kp_vector()            {_create(0);};
   kp_vector(const kp_vector&);
   
   //MATLAB code this(i1 + 1, i2)
   kp_vector(const kp_vector& v, int i1, int i2) {_create(0); init_from_vector(v, i1, i2);};
   kp_vector(const vector<double>&);
   kp_vector(int vec_size, double val); //create vector and fill it with val
   ~kp_vector()           {_clear();};
   void resize(int s_);
   
   int size() const {return s;};
   
   
   real&       operator[](int i)       {return d[i];};
   const real& operator[](int i) const {return d[i];};
   real&               el(int i)       {return d[i];};
   
   void operator=( const kp_vector&);
   void zeros();

   //apply operation to each element
   void operator /= (real val); 
   void operator *= (real val); 
   void operator += (real val); 
   void operator -= (real val); 
   
   //apply operation element by element
   void operator += (const kp_vector& v); 
   void operator -= (const kp_vector& v); 
   void operator *= (const kp_vector& v);  

   
   void init_from_matrix_column(const kp_matrix&, int col);
   
   //MATLAB code this(ind_begin + 1, int_end)
   void init_from_vector(const kp_vector&, int ind_begin, int int_end);
   
   //in idx we have indexes of v which we should use
   //MATLAB code: this = v(idx)
   void init_from_idx(const kp_vector& v, const vector<int>& idx);
   
   //subv is vector which was created by init_from_rowidx
   //MATLAB code: this(idx) = subv
   void set_from_subvector(const kp_vector& subv, const vector<int>& idx);
   
   //put subv after position
   void set_from_subvector(const kp_vector& subv, int start);
   
   //calulate variance (square of standart deviation)
   //analuges of MATLAB var
   //sum( (xi - mean(xi)) ) / (n - 1)
   real var();
   
   //calculate mean of vector
   real mean();
   
   //sum of square of elements
   //this' * this  (or this * this')
   real sum_sqr();
 private:
   void _create(int s_) {s = s_; d = new real[s];};
   void _clear();
 public:
   real* d;
   int s;
};

//v[i] = 1 / v[i]
void kp_inverse(kp_vector& v);

//v[i] = sqrt(v[i])
void sqrt(kp_vector& v);

double kp_calc_diff(const kp_vector& v1, const kp_vector& v2);


void kp_vertcat(kp_vector& rez, vector<kp_vector*> vs);
ostream& operator<<(ostream& out, kp_vector& v);

#endif
