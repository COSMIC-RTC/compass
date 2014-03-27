//kp_tlib.h
//some usefull functions (written by Thomas)


#ifndef __KP_TLIB_H__
#define __KP_TLIB_H__
#include "kp_smatrix.h"


// indr1, indc1, indr2, indc2 are ZERO-BASED
// ne peut etre utilise que lorsque les indices varient par colonne PUIS par ligne (column major)
void init_smatrix_from_larger_smatrix(kp_smatrix&M_to, const kp_smatrix&M_from, const int indr1, const int indc1, const int indr2, const int indc2);

#endif
