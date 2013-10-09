#include <yoga_phase.h>

yoga_phase::yoga_phase(yoga_context *current_context, long size)
{
  this->screen_size = size;
  
  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  dims_data2[1] = this->screen_size; 
  dims_data2[2] = this->screen_size; 
  
  this->d_screen = new yoga_obj<float>(current_context, dims_data2);
  delete[] dims_data2;
}

yoga_phase::~yoga_phase()
{
  delete this->d_screen;
}
