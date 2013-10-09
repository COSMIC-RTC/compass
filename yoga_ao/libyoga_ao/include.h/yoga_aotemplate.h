#ifndef _YOGA_AOTEMPLATE_H_
#define _YOGA_AOTEMPLATE_H_

#include <yoga_wfs.h>

using namespace std;

class yoga_aotemplate {
 public:
  int                      device; // # device
  string                   type;   // a name for your data
  long                     dim;    // # of elements

  yoga_obj<float>          *d_data;// the data
  yoga_obj<float>          *d_res; // the result


  yoga_context *current_context;   // the context in which it has been created 

 public:
  yoga_aotemplate(yoga_context *context, const char* type, long dim, int device);
  yoga_aotemplate(const yoga_aotemplate& aotemplate);
  ~yoga_aotemplate();

  int fill_data(float *idata);
  int fill_data();
  int do_compute();
};
template <class T> void comp_aotemplate(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif // _YOGA_AOTEMPLATE_H_
