/**
 * \file sutra_aotemplate.h
 *
 * \class sutra_aotemplate
 *
 * \ingroup libsutra
 *
 * \brief this class provides a class template to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_AOTEMPLATE_H_
#define _SUTRA_AOTEMPLATE_H_

#include <sutra_wfs.h>

class sutra_aotemplate {
 public:
  int device;        // # device
  std::string type;  // a name for your data
  long dim;          // # of elements

  carma_obj<float> *d_data;  // the data
  carma_obj<float> *d_res;   // the result

  carma_context *current_context;  // the context in which it has been created

 public:
  sutra_aotemplate(carma_context *context, const char *type, long dim,
                   int device);
  sutra_aotemplate(const sutra_aotemplate &aotemplate);
  ~sutra_aotemplate();

  int fill_data(float *idata);
  int fill_data();
  int do_compute();
};
template <class T>
void comp_aotemplate(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif  // _SUTRA_AOTEMPLATE_H_
