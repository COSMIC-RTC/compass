#include <yoga_wfs.h>
#include <yoga_ao_utils.h>
#include <yoga_utils.h>

yoga_thread_barrier thread_barrier;

yoga_wfs::yoga_wfs(yoga_context *context, const char* type, long nxsub, long nvalid, long npix, long nphase, 
		   long nrebin, long nfft,long ntot, long npup,float pdiam,float nphotons, int lgs, int device)
{
  this->d_camplipup     = 0L;
  this->d_camplifoc     = 0L;
  this->d_fttotim       = 0L;
  this->d_ftkernel      = 0L;
  this->d_pupil         = 0L;
  this->d_hrimg         = 0L;
  this->d_bincube       = 0L;
  this->d_binimg        = 0L;
  this->d_subsum        = 0L;
  this->d_offsets       = 0L;
  this->d_fluxPerSub    = 0L;
  this->d_sincar        = 0L;
  this->d_submask       = 0L;
  this->d_hrmap         = 0L;
  this->d_slopes        = 0L;
  this->image_telemetry = 0L;
  this->d_phasemap      = 0L;
  this->d_binmap        = 0L;
  this->d_validsubsx    = 0L;
  this->d_validsubsy    = 0L;
  this->d_istart        = 0L;
  this->d_jstart        = 0L;
  this->d_psum          = 0L;
  this->d_phalfxy       = 0L;
  this->d_poffsets      = 0L;
  this->pyr_cx          = 0L;
  this->pyr_cy          = 0L;
  
  this->type    = type;
  this->nxsub   = nxsub;
  this->nvalid  = nvalid;
  this->npix    = npix;
  this->nphase  = nphase;
  this->nrebin  = nrebin;
  this->nfft    = nfft;
  this->ntot    = ntot;
  this->npup    = npup;
  this->subapd  = pdiam;
  this->nphot   = nphotons;
  this->lgs     = (lgs == 1 ? true : false);
  this->device  = device; 
  this->nmaxhr  = nvalid;
  this->nffthr  = 1;

  this->kernconv       = false;

  // according to context determine number of hosts
  // and device per host
  // //////////////////////////////////////////////////////////////////
  // !!!!! just working for single or pair of devices per hosts !!!!!!!
  // //////////////////////////////////////////////////////////////////

  if ((context->get_ndevices() > 1) && (context->is_p2pactive(0,0)) && (context->is_p2pactive(0,1))) {
    this->ndevices = context->get_ndevices();
  } else this->ndevices = 1;

  //this->ndevices = 1;

  this->nhosts = context->get_nhosts();
  //int nslaves=(nhosts > 1?nhosts -1:0);

  cout << "yoga_wfs will use " <<  this->ndevices << 
    " devices on " << this->nhosts << " hosts"<<endl;

  // who am i ?
  if (device == 0) {
    this->is_master = true;
    this->child_id = 0;
  } else {
    this->is_master = false;
    this->child_id = device-1;
  }

  this->current_context = context;
  context->set_activeDevice(device);

  if (this->type == "sh") this->init_sh();
  if (this->type == "pyr") this->init_pyr();

  long *dims_data1 = new long[2]; dims_data1[0] = 1;
  long *dims_data2 = new long[3]; dims_data2[0] = 2; 

  this->nstreams = 1;//nvalid/10;
  while(nvalid % this->nstreams!=0) nstreams--;
  cerr << "wfs uses " << nstreams << " streams" << endl;
  streams = new yoga_streams(nstreams);

  if (this->is_master) {
    dims_data1[1] = 2 * nvalid;
    this->d_slopes = new yoga_obj<float>(context,dims_data1);

    dims_data2[1] = nphase*nphase; dims_data2[2] = nvalid; 
    this->d_phasemap = new yoga_obj<int>(context,dims_data2);

    if (this->ndevices > 1) {
      for (int i=1;i<this->ndevices;i++) {
	child_wfs.push_back(new yoga_wfs(context,type,nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,npup,
					 pdiam,nphot,lgs,i));
      }
    }
  }
}

int yoga_wfs::init_sh()
{
  long *dims_data1 = new long[2]; dims_data1[0] = 1;
  long *dims_data2 = new long[3]; dims_data2[0] = 2; 
  long *dims_data3 = new long[4]; dims_data3[0] = 3;

  if (this->is_master) {
    dims_data2[1] = npix*nxsub; dims_data2[2] =npix*nxsub ; 
    
    this->d_binimg = new yoga_obj<float>(current_context,dims_data2);
    // using 1 stream for telemetry
    this->image_telemetry = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK,1);
  }

  dims_data3[1] = npix; dims_data3[2] = npix;
  dims_data3[3] = nvalid/(this->ndevices*this->nhosts);

  this->d_bincube = new yoga_obj<float>(current_context,dims_data3);

  dims_data3[1] = nfft; dims_data3[2] = nfft;  
  dims_data3[3] = nvalid/(this->ndevices*this->nhosts);

  this->d_hrimg     = new yoga_obj<float>(current_context,dims_data3);
  this->d_submask   = new yoga_obj<float>(current_context,dims_data3);
  this->d_camplipup = new yoga_obj<cuFloatComplex>(current_context,dims_data3);
  this->d_camplifoc = new yoga_obj<cuFloatComplex>(current_context,dims_data3);

  int mdims[2]; 
  mdims[0] = (int)dims_data3[1]; mdims[1] = (int)dims_data3[2];
  cufftSafeCall(cufftPlanMany(&(this->pup_plan), 2 ,mdims,NULL,1,0,NULL,1,0,
			       CUFFT_C2C ,(int)dims_data3[3]));

  if (this->ntot != this->nfft) {
    // this is the big array => we use nmaxhr and treat it sequentially
    int mnmax = 500;
    if (nvalid > 2*mnmax) {
      this->nmaxhr = mnmax;
      int tmp0 = nvalid % mnmax;
      int tmp = 0;
      for (int cc = 1;cc<mnmax/5;cc++) {
	tmp = nvalid % (mnmax + cc);
	if ((tmp > tmp0) || (tmp == 0)) {
	  if (tmp == 0) tmp0 = 2*mnmax;
	  else tmp = tmp0;
	  this->nmaxhr  = mnmax + cc;
	}
      }
      this->nffthr  = (nvalid % this->nmaxhr == 0 ? 
		       nvalid/this->nmaxhr : nvalid/this->nmaxhr+1);
    } 
      
    dims_data3[1] = ntot; dims_data3[2] =ntot; dims_data3[3] = nmaxhr;  
    this->d_fttotim   = new yoga_obj<cuFloatComplex>(current_context,dims_data3);
    mdims[0] = (int)dims_data3[1];
    mdims[1] = (int)dims_data3[2];
    cufftSafeCall( cufftPlanMany(&(this->im_plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
				 (int)dims_data3[3]));
    
    dims_data1[1] = nfft * nfft;  
    this->d_hrmap = new yoga_obj<int>(current_context,dims_data1);
    
  } else {
    if (this->lgs) {
      dims_data3[1] = ntot; dims_data3[2] =ntot; dims_data3[3] = nvalid/(ndevices*nhosts);  
      this->d_fttotim   = new yoga_obj<cuFloatComplex>(current_context,dims_data3);
      mdims[0] = (int)dims_data3[1];
      mdims[1] = (int)dims_data3[2];
      cufftSafeCall( cufftPlanMany(&(this->im_plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
				   (int)dims_data3[3]));
    }
  }
  
  dims_data2[1] = ntot; dims_data2[2] = ntot; 
  this->d_ftkernel = new yoga_obj<cuFloatComplex>(current_context,dims_data2);
  
  dims_data2[1] = npup; dims_data2[2] = npup; 
  this->d_pupil = new yoga_obj<float>(current_context,dims_data2);
  
  dims_data2[1] = nphase; dims_data2[2] = nphase; 
  this->d_offsets = new yoga_obj<float>(current_context,dims_data2);
  
  dims_data1[1] = nxsub;
  this->d_istart = new yoga_obj<int>(current_context,dims_data1);
  this->d_jstart = new yoga_obj<int>(current_context,dims_data1);
  
  dims_data2[1] = nrebin*nrebin; dims_data2[2] = npix * npix; 
  this->d_binmap = new yoga_obj<int>(current_context,dims_data2);

  dims_data1[1] = nvalid/(this->ndevices*this->nhosts);
  this->d_validsubsx = new yoga_obj<int>(current_context,dims_data1);
  this->d_validsubsy = new yoga_obj<int>(current_context,dims_data1);
  this->d_subsum = new yoga_obj<float>(current_context,dims_data1);
  this->d_fluxPerSub = new yoga_obj<float>(current_context,dims_data1);


  return EXIT_SUCCESS;
}

int yoga_wfs::init_pyr()
{
  long *dims_data1 = new long[2]; dims_data1[0] = 1;
  long *dims_data2 = new long[3]; dims_data2[0] = 2; 
  long *dims_data3 = new long[4]; dims_data3[0] = 3;

  if (this->is_master) {
    dims_data2[1] = npix*nxsub + 3; // padding
    dims_data2[2] = npix*nxsub + 3; 
    this->d_binimg = new yoga_obj<float>(current_context,dims_data2);
    // using 1 stream for telemetry
    this->image_telemetry = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK,1);
  }

  dims_data3[1] = nfft; dims_data3[2] = nfft;  
  dims_data3[3] = 4;
  this->d_hrimg = new yoga_obj<float>(current_context,dims_data3);

  dims_data2[1] = ntot; dims_data2[2] = ntot; 
  this->d_submask = new yoga_obj<float>(current_context,dims_data2);
  this->d_camplipup = new yoga_obj<cuFloatComplex>(current_context,dims_data2);
  this->d_camplifoc = new yoga_obj<cuFloatComplex>(current_context,dims_data2);

  cufftSafeCall( cufftPlan2d(&(this->pup_plan),dims_data2[1],dims_data2[2],CUFFT_C2C ));

  dims_data3[1] = nfft/nrebin; dims_data3[2] = nfft/nrebin;

  this->d_bincube = new yoga_obj<float>(current_context,dims_data3);

  dims_data3[1] = nfft; dims_data3[2] =nfft; dims_data3[3] = 4; 
  this->d_fttotim   = new yoga_obj<cuFloatComplex>(current_context,dims_data3);
  int mdims[2]; 
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];
  //cufftHandle *plan=this->d_fttotim->getPlan();///< FFT plan
  cufftSafeCall( cufftPlanMany(&(this->im_plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
			       (int)dims_data3[3]));
  
  dims_data2[1] = ntot; dims_data2[2] = ntot; 
  this->d_pupil = new yoga_obj<float>(current_context,dims_data2);
  
  dims_data1[1] = npup;
  this->pyr_cx = new yoga_host_obj<int>(dims_data1,MA_WRICOMB);
  this->pyr_cy = new yoga_host_obj<int>(dims_data1,MA_WRICOMB);
  
  dims_data2[1] = nfft; dims_data2[2] = nfft; 
  this->d_poffsets = new yoga_obj<cuFloatComplex>(current_context,dims_data2);
  dims_data2[1] = ntot; dims_data2[2] = ntot; 
  this->d_phalfxy = new yoga_obj<cuFloatComplex>(current_context,dims_data2);
  dims_data2[1] = nfft; dims_data2[2] = nfft; 
  this->d_sincar = new yoga_obj<float>(current_context,dims_data2);

  dims_data1[1] = nvalid;
  this->d_psum = new yoga_obj<float>(current_context,dims_data1);

  dims_data1[1] = nvalid;
  this->d_validsubsx = new yoga_obj<int>(current_context,dims_data1);
  this->d_validsubsy = new yoga_obj<int>(current_context,dims_data1);
  this->d_subsum = new yoga_obj<float>(current_context,dims_data1);
  this->d_fluxPerSub = new yoga_obj<float>(current_context,dims_data1);

  return EXIT_SUCCESS;
}


yoga_wfs::~yoga_wfs()
{
  current_context->set_activeDevice(device);
  if (this->d_camplipup != 0L) delete this->d_camplipup;
  if (this->d_camplifoc != 0L) delete this->d_camplifoc;

  if (this->d_fttotim != 0L) delete this->d_fttotim;

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_pupil != 0L) delete this->d_pupil;
  if (this->d_hrimg != 0L) delete this->d_hrimg;
  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_subsum != 0L) delete this->d_subsum;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_submask != 0L) delete this->d_submask;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  cufftSafeCall( cufftDestroy(this->im_plan) );
  cufftSafeCall( cufftDestroy(this->pup_plan) );

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_binmap != 0L) delete this->d_binmap;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;
  if (this->d_istart != 0L) delete this->d_istart;
  if (this->d_jstart != 0L) delete this->d_jstart;

  if (this->d_psum != 0L) delete this->d_psum;
  if (this->d_phalfxy != 0L) delete this->d_phalfxy;
  if (this->d_poffsets != 0L) delete this->d_poffsets;
  if (this->pyr_cx != 0L) delete this->pyr_cx;
  if (this->pyr_cy != 0L) delete this->pyr_cy;
  
  if (this->lgs) delete this->d_gs->d_lgs;

  delete this->d_gs;

  delete this->streams;

  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    delete this->child_wfs[(this->child_wfs).size()-1];
    child_wfs.pop_back();
  } 
  //delete this->current_context;
}

int yoga_wfs::wfs_initgs(float xpos,float ypos,float lambda, float mag, long size,float noise,long seed)
{
  this->d_gs = new yoga_source(current_context,xpos,ypos,lambda,mag,size,"wfs",this->device);
  this->noise = noise;
  if (noise > -1) {
    this->d_bincube->init_prng(this->device);
    this->d_bincube->prng('N',noise,0.0f);
  }

  if (this->lgs) {
    this->d_gs->d_lgs  = new yoga_lgs(current_context,this->nvalid,this->ntot,this->nmaxhr);
    this->d_gs->lgs = this->lgs;
  }

  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    this->child_wfs[idx]->d_gs = new yoga_source(current_context,xpos,ypos,lambda,mag,size,"wfs",this->device);
    this->child_wfs[idx]->noise = noise;
    if (noise > -1) {
      this->child_wfs[idx]->d_bincube->init_prng(this->child_wfs[idx]->device);
      this->child_wfs[idx]->d_bincube->prng('N',noise,0.0f);
    }
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::wfs_initarrays(int *phasemap,int *hrmap, int *binmap,float *offsets, 
			     float *pupil, float *fluxPerSub, int *validsubsx, int *validsubsy, 
			     int *istart, int *jstart, cuFloatComplex *kernel)
{
  // copy only to master
  this->d_phasemap->host2device(phasemap);

  // need to copy to child
  this->d_fluxPerSub->host2device(fluxPerSub);
  this->d_offsets->host2device(offsets);
  this->d_pupil->host2device(pupil);
  this->d_binmap->host2device(binmap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);
  this->d_istart->host2device(istart);
  this->d_jstart->host2device(jstart);
  if (this->ntot != this->nfft) 
    this->d_hrmap->host2device(hrmap);
  this->d_ftkernel->host2device(kernel);

  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    int idx_valid = (idx + 1) * this->nvalid/(this->ndevices*this->nhosts);
    this->child_wfs[idx]->d_fluxPerSub->host2device(&(fluxPerSub[idx_valid]));
    this->child_wfs[idx]->d_offsets->host2device(offsets);
    this->child_wfs[idx]->d_pupil->host2device(pupil);
    this->child_wfs[idx]->d_binmap->host2device(binmap);
    this->child_wfs[idx]->d_validsubsx->host2device(&(validsubsx[idx_valid]));
    this->child_wfs[idx]->d_validsubsy->host2device(&(validsubsy[idx_valid]));
    this->child_wfs[idx]->d_istart->host2device(istart);
    this->child_wfs[idx]->d_jstart->host2device(jstart);
    if (this->ntot != this->nfft) 
      this->child_wfs[idx]->d_hrmap->host2device(hrmap);
    this->child_wfs[idx]->d_ftkernel->host2device(kernel);
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::wfs_initarrays(cuFloatComplex *halfxy,cuFloatComplex *offsets, float *focmask, 
			     float *pupil, int *cx, int *cy, float *sincar, int *phasemap, 
			     int *validsubsx, int *validsubsy)
{
  this->d_phalfxy->host2device(halfxy);
  this->d_poffsets->host2device(offsets);
  this->d_submask->host2device(focmask);
  this->d_pupil->host2device(pupil);
  this->pyr_cx->fill_from(cx);
  this->pyr_cy->fill_from(cy);
  this->d_sincar->host2device(sincar);
  this->d_phasemap->host2device(phasemap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);

  return EXIT_SUCCESS;
}

int yoga_wfs::load_kernels(float *lgskern)
{
  if (this->lgs) this->d_gs->d_lgs->load_kernels(lgskern,this->device);

  return EXIT_SUCCESS;
}

int yoga_wfs::sensor_trace(yoga_atmos *yatmos)
{
  //do raytracing to get the phase
  this->d_gs->raytrace(yatmos);
  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    cutilSafeCall(cudaMemcpy(this->child_wfs[idx]->d_gs->d_phase->d_screen->getData(), 
			     this->d_gs->d_phase->d_screen->getData(), 
			     sizeof(float)*this->d_gs->d_phase->d_screen->getNbElem(),
			     cudaMemcpyDeviceToDevice));
  }

  // same with textures
  // apparently not working for large sizes...
  //this->d_gs->raytrace_shm(yatmos);

  return EXIT_SUCCESS;
}

int yoga_wfs::sensor_trace(yoga_dms *ydm, int rst)
{
  //do raytracing to get the phase
  this->d_gs->raytrace(ydm,rst);
  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    cutilSafeCall(cudaMemcpy(this->child_wfs[idx]->d_gs->d_phase->d_screen->getData(), 
			     this->d_gs->d_phase->d_screen->getData(), 
			     sizeof(float)*this->d_gs->d_phase->d_screen->getNbElem(),
			     cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::comp_sh_psf(cuFloatComplex *data_out,cuFloatComplex *data_inter,float *phase,
			  float *offsets,float *pupil, float scale,
			  int *istart,int *jstart,int *validx,int *validy,
			  int nphase, int screen_size, int nshim, int nssp,
			  cufftHandle plan, int device) 
{
  // segment phase and fill cube of complex ampli with exp(i*phase_seg)
  fillcamplipup(data_inter, 
		phase, 
		offsets, 
		pupil, 
		scale, 
		istart, 
		jstart, 
		validx, 
		validy, 
		nphase, 
		screen_size, 
		nshim, 
		nphase * nphase * nssp, 
		device);

  // do fft of the cube  
  yoga_fft(data_inter,data_out, 1,plan);

  // get the hrimage by taking the | |^2
  // keep it in amplifoc to save mem space
  abs2c(data_out,data_out, nshim * nshim * nssp, device);

  return EXIT_SUCCESS;
}

int yoga_wfs::kernelconv_sh_psf(cuFloatComplex *data_out, cuFloatComplex *kernel, int nout, cufftHandle plan, 
				int device) 
{
  yoga_fft(data_out, data_out, 1, plan);
  
  convolve(data_out,kernel,nout,device);
  
  yoga_fft(data_out, data_out, -1, plan);

  return EXIT_SUCCESS;
}
 
int yoga_wfs::kernelconv_sh_psf(cuFloatComplex *data_out, cuFloatComplex *kernel, int nout, int nkernel, cufftHandle plan, 
				int device) 
{
  yoga_fft(data_out, data_out, 1, plan);
  
  convolve_cube(data_out, kernel, nout, nkernel, device);
  
  yoga_fft(data_out, data_out, -1,plan);

  return EXIT_SUCCESS;
}

int yoga_wfs::comp_sh_generic() 
{
  //current_context->set_activeDevice(device);
  int nssp;
  if (this->ndevices > 1)
    nssp  = this->nvalid / (this->ndevices * this->nhosts);
  else
    nssp  = this->nvalid;

  comp_sh_psf(this->d_camplifoc->getData(),
	      this->d_camplipup->getData(),
	      this->d_gs->d_phase->d_screen->getData(),
	      this->d_offsets->getData(),
	      this->d_pupil->getData(), 
	      this->d_gs->scale,
	      this->d_istart->getData(),
	      this->d_jstart->getData(),
	      this->d_validsubsx->getData(),
	      this->d_validsubsy->getData(),
	      this->nphase, 
	      this->d_gs->d_phase->d_screen->getDims(1),
	      this->nfft, 
	      nssp,
	      this->pup_plan,
	      this->device); 

  //set bincube to 0 or noise
  cutilSafeCall(cudaMemset(this->d_bincube->getData(), 0,
			   sizeof(float) * this->d_bincube->getNbElem()));

  // increase fov if required
  // and fill bincube with data from hrimg
  // we need to do this sequentially if nvalid > nmaxhr to
  // keep raesonable mem occupancy
  if (this->ntot != this->nfft) {
    for (int cc = 0; cc < this->nffthr; cc++) {
      cutilSafeCall(cudaMemset(this->d_fttotim->getData(), 0,
			       sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));
      
      int indxstart1, indxstart2, indxstart3;
      
      if ((cc == this->nffthr - 1)
	  && (this->nvalid % this->nmaxhr != 0)) {
	indxstart1 = this->d_camplifoc->getNbElem()
	  - this->nfft * this->nfft * this->nmaxhr;
	if (this->lgs)
	  indxstart2 = this->ntot * this->nvalid
	    - this->ntot * this->nmaxhr;
	indxstart3 = this->d_bincube->getNbElem()
	  - this->npix * this->npix * this->nmaxhr;
      } else {
	indxstart1 = this->nfft * this->nfft * this->nmaxhr * cc;
	if (this->lgs)
	  indxstart2 = this->ntot * this->nmaxhr * cc;
	indxstart3 = this->npix * this->npix * this->nmaxhr * cc;
      }
      
      cuFloatComplex *data = this->d_camplifoc->getData();
      indexfill(this->d_fttotim->getData(),
		&(data[indxstart1]),
		this->d_hrmap->getData(), this->nfft, this->ntot,
		this->nfft * this->nfft * this->nmaxhr, this->device);
      
      if (this->lgs) {
	// compute lgs spot on the fly from binned profile image
	this->d_gs->d_lgs->lgs_makespot(this->device, indxstart2);
	// convolve with psf
	kernelconv_sh_psf(this->d_fttotim->getData(), this->d_gs->d_lgs->d_ftlgskern->getData(), 
			  this->d_fttotim->getNbElem(),this->im_plan, this->device); 
      }
      
      if (this->kernconv) {
	kernelconv_sh_psf(this->d_fttotim->getData(), this->d_ftkernel->getData(), this->d_fttotim->getNbElem(),
			  this->d_ftkernel->getNbElem(),this->im_plan, this->device); 
      }
      
      float *data2=this->d_bincube->getData();
      if (this->nstreams > 1) 
	fillbincube_async(this->streams, &(data2[indxstart3]),
			  this->d_fttotim->getData(), this->d_binmap->getData(),
			  this->ntot * this->ntot, this->npix * this->npix,
			  this->nrebin * this->nrebin, this->nmaxhr, this->device);
      else
	fillbincube(&(data2[indxstart3]),this->d_fttotim->getData(), this->d_binmap->getData(),
		    this->ntot * this->ntot, this->npix * this->npix,
		    this->nrebin * this->nrebin, this->nmaxhr, this->device);
    }
  } else {
    if (this->lgs) {
      this->d_gs->d_lgs->lgs_makespot(this->device, 0);

      kernelconv_sh_psf(this->d_camplifoc->getData(), 
			this->d_gs->d_lgs->d_ftlgskern->getData(), 
			this->d_fttotim->getNbElem(),
			this->im_plan, 
			this->device); 

      if (this->nstreams > 1) 
	fillbincube_async(this->streams, 
			  this->d_bincube->getData(), 
			  this->d_camplifoc->getData(),
			  this->d_binmap->getData(), 
			  this->nfft * this->nfft,
			  this->npix * this->npix, 
			  this->nrebin * this->nrebin,
			  this->nvalid, 
			  this->device);
      else
	fillbincube(this->d_bincube->getData(), 
		    this->d_camplifoc->getData(),
		    this->d_binmap->getData(), 
		    this->nfft * this->nfft,
		    this->npix * this->npix, 
		    this->nrebin * this->nrebin,
		    this->nvalid, 
		    this->device);
    } else {
      if (this->kernconv) {

	kernelconv_sh_psf(this->d_camplifoc->getData(), 
			  this->d_ftkernel->getData(), 
			  this->d_camplifoc->getNbElem(),
			  this->d_ftkernel->getNbElem(),
			  this->pup_plan, 
			  this->device); 

      }

      if (this->nstreams > 1) 
    	fillbincube_async(this->streams, 
			  this->d_bincube->getData(), 
			  this->d_camplifoc->getData(),
			  this->d_binmap->getData(), 
			  this->nfft * this->nfft,
			  this->npix * this->npix, 
			  this->nrebin * this->nrebin,
			  nssp, 
			  this->device);
      else
    	fillbincube(this->d_bincube->getData(), 
		    this->d_camplifoc->getData(),
		    this->d_binmap->getData(), 
		    this->nfft * this->nfft,
		    this->npix * this->npix, 
		    this->nrebin * this->nrebin,
		    nssp, 
		    this->device);
    }
    
  }
  // normalize images :
  // get the sum value per subap

  if (this->nstreams > 1) 
    subap_reduce_async(this->npix * this->npix, 
		       nssp, 
		       this->streams,
		       this->d_bincube->getData(), 
		       this->d_subsum->getData());
  else {
    subap_reduce(this->d_bincube->getNbElem(), 
		 this->npix * this->npix,
		 nssp, 
		 this->d_bincube->getData(), 
		 this->d_subsum->getData());
  }

  if (this->nstreams > 1) 
    subap_norm_async(this->d_bincube->getData(), 
		     this->d_bincube->getData(),
		     this->d_fluxPerSub->getData(), 
		     this->d_subsum->getData(), 
		     this->nphot,
		     this->npix * this->npix, 
		     this->npix * this->npix * nssp,
		     this->streams,
		     this->device);
  else {
    // multiply each subap by nphot*fluxPersub/sumPerSub
    subap_norm(this->d_bincube->getData(), 
	       this->d_bincube->getData(),
	       this->d_fluxPerSub->getData(), 
	       this->d_subsum->getData(), 
	       this->nphot,
	       this->npix * this->npix, 
	       this->npix * this->npix * nssp, 
	       this->device);
  }

  // add noise
  if (this->noise > -1) {
    //cout << "adding poisson noise" << endl;
    this->d_bincube->prng('P');
  }
  if (this->noise > 0) {
    //cout << "adding detector noise" << endl;
    this->d_bincube->prng('N',this->noise,1.0f);
  }
  return EXIT_SUCCESS;
}

int yoga_wfs::comp_pyr_generic() 
{
  pyr_getpup(this->d_camplipup->getData(),this->d_gs->d_phase->d_screen->getData(), this->d_phalfxy->getData(),
	     this->d_pupil->getData(),this->ntot , this->device);
  
  yoga_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), -1,
	   this->pup_plan);
  
  pyr_submask(this->d_camplifoc->getData(), this->d_submask->getData(),this->ntot, this->device);

  cutilSafeCall(cudaMemset(this->d_hrimg->getData(), 0,
			   sizeof(float) * this->d_hrimg->getNbElem()));

  //this->npup = 1;
  for (int cpt=0;cpt<this->npup;cpt++) {
    // modulation loop
    // computes the high resolution images
    cutilSafeCall(cudaMemset(this->d_fttotim->getData(), 0,
			     sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));
    
    pyr_rollmod(this->d_fttotim->getData(),this->d_camplifoc->getData(), this->d_poffsets->getData()
		,(this->pyr_cx->getData())[cpt],(this->pyr_cy->getData())[cpt],this->ntot , this->nfft,
		this->device);
    /*
    pyr_rollmod(this->d_fttotim->getData(),this->d_camplifoc->getData(), this->d_poffsets->getData(),0,
		0,this->ntot , this->nfft, this->device);
    */
    
    yoga_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
	     this->im_plan);

    float fact = 1.0f/this->nfft/this->nfft/this->nfft/2.0;
    //if (cpt == this->npup-1) fact = fact / this->npup;
    
    pyr_abs2(this->d_hrimg->getData(), this->d_fttotim->getData(), fact, this->nfft, 4, this->device);
  }
  /*
  // spatial filtering by the pixel extent:
  yoga_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
	   *this->d_fttotim->getPlan());

  pyr_submask3d(this->d_fttotim->getData(), this->d_sincar->getData(),this->nfft, 4, this->device);

  yoga_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
	   *this->d_fttotim->getPlan());
  
  pyr_abs(this->d_hrimg->getData(), this->d_fttotim->getData(),this->nfft, 4, this->device);

  pyr_fact(this->d_hrimg->getData(),1.0f/this->nfft/this->nfft,this->nfft,4,this->device);
	   */	   
    
  if (this->noise > 0) {
    this->d_bincube->prng('N', this->noise);
  } else
    cutilSafeCall(cudaMemset(this->d_bincube->getData(), 0,
			     sizeof(float) * this->d_bincube->getNbElem()));

  pyr_fillbin(this->d_bincube->getData(),this->d_hrimg->getData(),this->nrebin,
	      this->nfft,this->nfft/this->nrebin,4, this->device);
  
  pyr_subsum(this->d_subsum->getData(),this->d_bincube->getData(),this->d_validsubsx->getData(),
	     this->d_validsubsy->getData(),this->nfft/this->nrebin,this->nvalid, 4, this->device);
  
  
  int blocks,threads;
  getNumBlocksAndThreads(this->device,this->nvalid, blocks,threads);
  reduce(this->nvalid,threads,blocks,this->d_subsum->getData(),this->d_subsum->getData());

  pyr_fact(this->d_bincube->getData(),this->nphot,this->d_subsum->getData(),this->nfft/this->nrebin,4,this->device);

  pyr_subsum(this->d_subsum->getData(),this->d_bincube->getData(),this->d_validsubsx->getData(),
	     this->d_validsubsy->getData(),this->nfft/this->nrebin,this->nvalid, 4, this->device);
  
  return EXIT_SUCCESS;
}

int yoga_wfs::comp_image()
{
  
  int nthreads = this->ndevices;

  thread_barrier = yoga_create_barrier(nthreads);
  yoga_thread thread[nthreads];

  thread[0] = yoga_start_thread(comp_image_thread,(void *)(this));
  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    thread[idx+1] = yoga_start_thread(comp_image_thread,(void *)(this->child_wfs[idx]));
  }

  yoga_wait4barrier(&thread_barrier);

  for (size_t idx = 0; idx < (this->child_wfs).size()+1; idx++) {
    yoga_destroy_thread(thread[idx]);
  }

  yoga_destroy_barrier(&thread_barrier);

  return EXIT_SUCCESS;
}

void* comp_image_thread(void *thread_data) {
  yoga_wfs *wfs =(yoga_wfs *)thread_data;
  
  wfs->current_context->set_activeDevice(wfs->device);

  //int result;
  if (wfs->type == "sh")
    wfs->comp_sh_generic();
  

  if (wfs->type == "pyr") 
    wfs->comp_pyr_generic();

  yoga_increment_barrier(&thread_barrier);

  return (void*)0L;
}

int yoga_wfs::comp_image_tele()
{
  int nssp;
  if (this->ndevices > 1)
    nssp  = this->nvalid / (this->ndevices * this->nhosts);
  else
    nssp  = this->nvalid;
  
  fillbinimg_async(this->image_telemetry,
		   this->d_binimg->getData(),this->d_bincube->getData(),this->npix,nssp,
		   this->npix*this->nxsub,this->d_validsubsx->getData(),this->d_validsubsy->getData(),
		   this->d_binimg->getNbElem(),false,this->device);
  
  for (size_t idx = 0; idx < (this->child_wfs).size(); idx++) {
    fillbinimg_async(this->image_telemetry,
		     this->d_binimg->getData(),this->child_wfs[idx]->d_bincube->getData(),this->npix,nssp,
		     this->npix*this->nxsub,this->child_wfs[idx]->d_validsubsx->getData(),
		     this->child_wfs[idx]->d_validsubsy->getData(),
		     this->d_binimg->getNbElem(),false,this->device);
  } 

  int result = comp_image();
  
  return result;
}

int yoga_wfs::slopes_geom(int type, float *slopes)
{
/* 
normalization notes :
 σ² = 0.17 (λ/D)^2 (D/r_0)^(5/3) , σ² en radians d'angle
 σ = sqrt(0.17 (λ/D)^2 (D/r_0)^(5/3)) * 206265 , σ en secondes

    // computing subaperture phase difference at edges

 todo : integrale( x * phase ) / integrale (x^2);
 with x = span(-0.5,0.5,npixels)(,-:1:npixels) * subap_diam * 2 * pi / lambda / 0.206265
*/
  if (type == 0) {
    // this is to convert in arcsec
    //> 206265* 0.000001/ 2 / 3.14159265 = 0.0328281

    float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    phase_reduce(this->nphase,this->nvalid,this->d_gs->d_phase->d_screen->getData(),
		 slopes,this->d_phasemap->getData(),alpha);
  }

  if (type == 1) {
    float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    phase_derive(this->nphase * this->nphase * this->nvalid,this->nphase * this->nphase,
		 this->nvalid,this->nphase,this->d_gs->d_phase->d_screen->getData(),
		 slopes,this->d_phasemap->getData(),this->d_pupil->getData(),
		 alpha,this->d_fluxPerSub->getData());
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::slopes_geom(int type)
{
  this->slopes_geom(type,this->d_slopes->getData());

  return EXIT_SUCCESS;
}



yoga_sensors::yoga_sensors(yoga_context *context, const char* type, int nwfs,long *nxsub,long *nvalid,long *npix,long *nphase, long *nrebin,
			   long *nfft, long *ntot,long npup,float *pdiam, float *nphot,  int *lgs, int device)
{
  this->nsensors = nwfs;

  for (int i=0;i<nwfs;i++) {
    d_wfs.push_back(new yoga_wfs(context,type,nxsub[i],nvalid[i],npix[i],nphase[i],nrebin[i],nfft[i],ntot[i],npup,
				 pdiam[i],nphot[i],lgs[i],device));
  }
}

yoga_sensors::~yoga_sensors()
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    delete this->d_wfs[(this->d_wfs).size()-1];
    d_wfs.pop_back();
  } 
}

int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, 
				 float *noise, long *seed)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],noise[idx],
				   seed[idx]);
  } 
  return EXIT_SUCCESS;
}
int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, 
				 float *noise)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],noise[idx],
				   1234*idx);
  } 
  return EXIT_SUCCESS;
}
int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],-1,1234);
  } 
  return EXIT_SUCCESS;
}



