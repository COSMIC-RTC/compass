#include <yoga_dm.h>

yoga_dm::yoga_dm(yoga_context *context, const char* type, long dim, long ninflu, long influsize, long ninflupos, long n_npoints, float push4imat, int device)
{

  this->d_influ = NULL; 
  this->d_influpos = NULL;
  this->d_npoints = NULL;
  this->d_istart = NULL;
  this->d_xoff = NULL;
  this->d_yoff = NULL;
  this->d_coeffs = NULL;
  this->d_mask = NULL;
  this->d_zr = NULL;
  this->d_ztheta = NULL;

  this->current_context = context;
  this->ninflu          = ninflu;
  this->dim             = dim;
  this->influsize       = influsize;
  this->d_shape         = new yoga_phase(context, dim);
  this->type            = type;
  this->push4imat       = push4imat;
  this->device          = device;
  
  long *dims_data1 = new long[2];
  dims_data1[0] = 1; dims_data1[1] = ninflu;     
  this->d_comm = new yoga_obj<float>(context, dims_data1);

  if (strcmp(type, "kl") != 0) {
    long *dims_data3 = new long[4];
    dims_data3[0] = 3; dims_data3[1] = influsize; dims_data3[2] = influsize;  dims_data3[3] = ninflu;
    this->d_influ = new yoga_obj<float>(context, dims_data3);
    delete[] dims_data3;

    dims_data1[1] = ninflu;     
    this->d_xoff = new yoga_obj<int>(context, dims_data1);
    this->d_yoff = new yoga_obj<int>(context, dims_data1);
  }

  if (strcmp(type, "pzt")==0) {
    dims_data1[1] = ninflupos;     
    this->d_influpos = new yoga_obj<int>(context, dims_data1);

    dims_data1[1] = n_npoints;     
    this->d_npoints = new yoga_obj<int>(context, dims_data1);
    this->d_istart = new yoga_obj<int>(context, dims_data1);
  }

  if (strcmp(type, "kl")==0) {
    this->d_kl    = new yoga_kl(context, influsize, ninflupos, n_npoints, ninflu, device);
  }

  delete[] dims_data1;

}

yoga_dm::~yoga_dm()
{
  //delete this->current_context;

  delete this->d_shape;
  delete this->d_comm;

  if (this->d_influ != NULL) delete this->d_influ; 
  if (this->d_influpos != NULL) delete this->d_influpos;
  if (this->d_npoints != NULL) delete this->d_npoints;
  if (this->d_istart != NULL) delete this->d_istart;
  if (this->d_xoff != NULL) delete this->d_xoff;
  if (this->d_yoff != NULL) delete this->d_yoff;
  if (this->d_coeffs != NULL) delete this->d_coeffs;
  if (this->d_mask != NULL) delete this->d_mask;
  if (this->d_zr != NULL) delete this->d_zr;
  if (this->d_ztheta != NULL) delete this->d_ztheta;
}


int yoga_dm::pzt_loadarrays(float *influ, int *influpos,int *npoints, int *istart, int *xoff, int *yoff)
{
  this->d_influ->host2device(influ);
  this->d_xoff->host2device(xoff);
  this->d_yoff->host2device(yoff);
  this->d_istart->host2device(istart);     // size = dim x dim
  this->d_influpos->host2device(influpos); // size ???
  this->d_npoints->host2device(npoints);   // size : dim x dim

  return EXIT_SUCCESS;
}

int yoga_dm::kl_loadarrays(float *rabas, float *azbas,int *ord, float *cr, float *cp)
{
  this->d_kl->d_rabas->host2device(rabas);
  this->d_kl->d_azbas->host2device(azbas);
  this->d_kl->h_ord->fill_from(ord);
  this->d_kl->d_ord->host2device(ord);
  this->d_kl->d_cr->host2device(cr);
  this->d_kl->d_cp->host2device(cp);

  return EXIT_SUCCESS;
}

int yoga_dm::reset_shape()
{
  current_context->set_activeDevice(device);

  cutilSafeCall(cudaMemset(this->d_shape->d_screen->getData(), 0,
			   sizeof(float)*this->d_shape->d_screen->getNbElem()));

  return EXIT_SUCCESS;
}

int yoga_dm::comp_shape(float *comvec)
{
  current_context->set_activeDevice(device);
  this->reset_shape();

  int nthreads = 0,nblocks = 0;
  getNumBlocksAndThreads(this->device,this->d_shape->d_screen->getNbElem(),nblocks, nthreads);

  if (this->type == "pzt")
    comp_dmshape(nthreads,nblocks, this->d_influ->getData(), this->d_shape->d_screen->getData(),
		 this->d_influpos->getData(), this->d_istart->getData(), this->d_npoints->getData(),
		 comvec, this->influsize*this->influsize,this->d_shape->d_screen->getNbElem());

  if (this->type == "tt")
    comp_fulldmshape(nthreads,nblocks, this->d_influ->getData(), this->d_shape->d_screen->getData(),
		     this->ninflu, this->influsize*this->influsize,comvec,this->d_shape->d_screen->getNbElem());

  if (this->type == "kl") {
    int xoff = (int)((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)/2.0f);
    int yoff = xoff;
    this->d_kl->do_combi(comvec,this->d_shape->d_screen->getData(),this->d_shape->d_screen->getDims()[1],xoff,yoff);
  }
  return EXIT_SUCCESS;
}

int yoga_dm::comp_shape()
{
  return this->comp_shape(this->d_comm->getData());
}

int yoga_dm::comp_oneactu(int nactu, float ampli)
{
  this->reset_shape();
  int nthreads = 0,nblocks = 0;
  //getNumBlocksAndThreads(this->device,this->dim * this->dim, nblocks, nthreads);
  getNumBlocksAndThreads(this->device,this->influsize * this->influsize, nblocks, nthreads);
  if (this->type == "pzt")
    oneactu(nthreads,nblocks, this->d_influ->getData(), this->d_shape->d_screen->getData(),
	    nactu, ampli, this->d_xoff->getData(), this->d_yoff->getData(), this->dim,this->influsize,
	    this->d_shape->d_screen->getNbElem());
  if (this->type == "tt")
    oneactu(nthreads,nblocks, this->d_influ->getData(), this->d_shape->d_screen->getData(),
	    nactu, ampli, this->dim,this->influsize,this->d_shape->d_screen->getNbElem());
  if (this->type == "kl") {
    int xoff = (int)((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)/2.0f);
    int yoff = xoff;
    this->d_kl->do_compute(ampli,this->d_shape->d_screen->getData(),nactu,
			   this->d_shape->d_screen->getDims()[1],xoff,yoff);
  }

  return EXIT_SUCCESS;
}

yoga_dms::yoga_dms(int ndm)
{
  this->ndm     = ndm;
}

yoga_dms::~yoga_dms()
{
  for( std::map< type_screen,yoga_dm * >::iterator it = this->d_dms.begin();
       this->d_dms.end()!= it ; it++) {
    delete it->second;
  }
  this->d_dms.clear();
}

int yoga_dms::add_dm(yoga_context *context, const char* type,float alt, long dim, long ninflu, long influsize, long ninflupos, long n_npoints, float push4imat, int device)
{
  this->d_dms[make_pair(type,alt)] = new yoga_dm(context, type,dim,ninflu,influsize,ninflupos,n_npoints,push4imat,device);

  return EXIT_SUCCESS;
}

int yoga_dms::remove_dm(const char* type,float alt)
{
  delete this->d_dms[make_pair(type,alt)];
  this->d_dms.erase(make_pair(type,alt));
  
  return EXIT_SUCCESS;
}

