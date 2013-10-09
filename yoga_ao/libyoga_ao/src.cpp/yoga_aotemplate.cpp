#include <yoga_aotemplate.h>

yoga_aotemplate::yoga_aotemplate(yoga_context *context, const char* type, long dim, int device)
{
  // some inits
  this->current_context=context;
  this->dim        = dim;
  this->type       = type;
  this->device     = device;

  // allocate data and result
  long *dims_data1 = new long[2];
  dims_data1[0] = 1; dims_data1[1] = dim;     
  this->d_data     = new yoga_obj<float>(context, dims_data1);
  this->d_res      = new yoga_obj<float>(context, dims_data1);
  delete[] dims_data1;

}

yoga_aotemplate::~yoga_aotemplate()
{
  // delete data and result
  delete this->d_data;
  delete this->d_res;
}


int yoga_aotemplate::fill_data(float *idata)
{
  // fill data with an external array
  this->d_data->host2device(idata);

  return EXIT_SUCCESS;
}

int yoga_aotemplate::fill_data()
{
  // fill data with random numbers
  this->d_data->init_prng(this->device);
  this->d_data->prng('N');

  return EXIT_SUCCESS;
}

int yoga_aotemplate::do_compute()
{
  // do computation on data and store in result
  int nthreads = 0,nblocks = 0;
  getNumBlocksAndThreads(this->device,this->dim, nblocks, nthreads);

  comp_aotemplate(nthreads,nblocks, this->d_data->getData(),this->d_res->getData(),this->d_data->getNbElem());

  return EXIT_SUCCESS;
}

