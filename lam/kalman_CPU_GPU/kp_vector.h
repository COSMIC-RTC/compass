//kp_vector.h
//Kalman project vector

#ifndef __SEGER__KP_VECTOR_H__
#define __SEGER__KP_VECTOR_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include <carma_host_obj.h>

#include "kp_matrix.h"
#include <typeinfo>

template<typename real> class kp_matrix;                                                                                  

using namespace std;


template<typename real>                                                                                   
class kp_vector
{
	template<typename T> friend class kp_vector;

#ifndef KP_WITH_CARMA
 public:
   int size() const {return s;};
   real* getData()const{return d;};
   real* getData(int i)const{return &(d[i]);};
   //MATLAB code this(i1 + 1, i2)
   kp_vector(const kp_vector<real>& v, int i1, int i2) {_create(0); init_from_vector(v, i1, i2);};
   kp_vector(const vector<double>&);
   kp_vector(int vec_size, double val); //create vector and fill it with val
     
   real& operator[](int i){return d[i];};
   const real& operator[](int i) const {return d[i];};
   

 private:
   real& el(int i) {return d[i];};

   


#else
 public:
   int size() const {return carma_host_object->getDims(1)};
   real* getData()const{return carma_host_object->getData();};
   real* getData(int i)const{return carma_host_object->getData(i);};
   real& operator[](int i){return *(carma_host_object->getData(i));};
   const real& operator[](int i) const {return *(carma_host_object->getData(i));};


 private:
   carma_host_obj<real>* carma_host_object;
   real& el(int i){return *(carma_host_object->getData(i));};

#endif

   
// same methods (do not depend on whether KP_WITH_CARMA is defined or not)	
 public:
   kp_vector();
   kp_vector(int s_);
   ~kp_vector();
   //apply operation to each element
   void operator /= (real val); 
   void operator *= (real val); 
   void operator += (real val); 
   void operator -= (real val); 
   //apply operation element by element
   void operator += (const kp_vector<real>& v); 
   void operator -= (const kp_vector<real>& v); 
   void operator *= (const kp_vector<real>& v);  
   void init_from_matrix_column(const kp_matrix<real>&, int col);
   //MATLAB code this(ind_begin + 1, int_end)
   void init_from_vector(const kp_vector<real>&, int ind_begin, int int_end);
   //in idx we have indexes of v which we should use
   //MATLAB code: this = v(idx)
   void init_from_idx(const kp_vector<real>& v, const vector<int>& idx);
   //subv is vector which was created by init_from_rowidx
   //MATLAB code: this(idx) = subv
   void set_from_subvector(const kp_vector<real>& subv, const vector<int>& idx);
   //put subv after position
   void set_from_subvector(const kp_vector<real>& subv, int start);
   //calulate variance (square of standard deviation)
   //analuges of MATLAB var
   //sum( (xi - mean(xi)) ) / (n - 1)
   real var();
   //calculate mean of vector
   real mean();
   //sum of square of elements
   //this' * this  (or this * this')
   real sum_sqr();


// each method exists twice : one if KP_WITH_CARMA is not defined ; another one if KP_WITH_CARMA is defined
   template<typename T>  kp_vector(const kp_vector<T>&);
   kp_vector(const kp_vector<real>&);
   void resize(int s_);
   template<typename T> void operator=( const kp_vector<T>&);
   void operator=( const kp_vector<real>&);
   void zeros();



 private:
// each method exists twice : one if KP_WITH_CARMA is not defined ; another one if KP_WITH_CARMA is defined
   void _create(int s_);
   void _clear();

   real* d;
   int s;




};

//v[i] = 1 / v[i]
template<typename real> 
void kp_inverse(kp_vector<real>& v);

//v[i] = sqrt(v[i])
template<typename real> 
void sqrt(kp_vector<real>& v);

template<typename real>
double kp_calc_diff(const kp_vector<real>& v1, const kp_vector<real>& v2);


template<typename real> 
void kp_vertcat(kp_vector<real>& rez, vector<kp_vector<real>*> vs);

template<typename real>
ostream& operator<<(ostream& out, kp_vector<real>& v);

#include "kp_vector.hpp"

#endif
