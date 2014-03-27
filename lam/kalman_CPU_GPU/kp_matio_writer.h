//kp_matio_writer.h
//class for write in MATLAB file


#ifndef __SEGER__KP_MATIO_WRITER_H__
#define __SEGER__KP_MATIO_WRITER_H__

#include <iostream>
#include "kp_vector.h"
#include "kp_matrix.h"
#include <string>
#include <matio.h>

class kp_matio_writer
{
 public:
   kp_matio_writer(string fn);  //open file
   ~kp_matio_writer();          //close file
   void add(string name, const kp_vector&);
   void add(string name, const kp_matrix&);   
 private:
   mat_t* matfp;
};

#endif
