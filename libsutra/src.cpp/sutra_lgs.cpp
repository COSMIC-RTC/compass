#include <sutra_ao_utils.h>
#include <sutra_lgs.h>
#include<sutra_wfs.h>

sutra_lgs::sutra_lgs(carma_context *context, sutra_sensors *sensors, long nvalid, long npix,
    long nmaxhr) {
  this->current_context = context;
  this->device = current_context->get_activeDevice();

  this->nvalid = nvalid;
  this->npix = npix;
  this->nprof = nprof;
  this->nmaxhr = nmaxhr;
  this->d_doffaxis = this->d_prof1d = this->d_profcum = 0;
  this->hg = this->pixsize = this->h0 = this->deltah = 0;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data1[1] = npix;
  this->d_beam = new carma_obj<float>(context, dims_data1);
  this->d_ftbeam = new carma_obj<cuFloatComplex>(context, dims_data1);

  dims_data1[1] = nvalid;
  this->d_azimuth = new carma_obj<float>(context, dims_data1);

  dims_data2[1] = npix;
  dims_data2[2] = nvalid;
  this->d_prof2d = new carma_obj<cuFloatComplex>(context, dims_data2);

  int mdims[2];
  mdims[0] = (int) dims_data2[1];

  cufftHandle *plan = this->d_prof2d->getPlan(); ///< FFT plan
  carmafftSafeCall(
      cufftPlanMany(plan, 1 ,mdims,NULL,1,0,NULL,1,0, CUFFT_C2C ,(int)dims_data2[2]));

  dims_data3[1] = npix;
  dims_data3[2] = npix;
  dims_data3[3] = nmaxhr;
  //this->d_lgskern = new carma_obj<float>(context, dims_data3);
  //this->d_ftlgskern = new carma_obj<cuFloatComplex>(context, dims_data3);
  this->d_lgskern = sensors->d_lgskern;
  this->d_ftlgskern = sensors->d_ftlgskern;

  mdims[0] = (int) dims_data3[1];
  mdims[1] = (int) dims_data3[2];

  vector<int> vdims (dims_data3+1,dims_data3 + 4);

  if (sensors->ftlgskern_plans.find(vdims) == sensors->ftlgskern_plans.end()) {
	  //DEBUG_TRACE("Creating FFT plan : %d %d %d",mdims[0],mdims[1],dims_data3[3]);printMemInfo();
	  cufftHandle *plan = (cufftHandle*)malloc(sizeof(cufftHandle));// = this->d_camplipup->getPlan(); ///< FFT plan
	  carmafftSafeCall(
			cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, CUFFT_C2C ,(int)dims_data3[3]));

	  sensors->ftlgskern_plans.insert(pair<vector<int>, cufftHandle*>(vdims, plan));

	  this->ftlgskern_plan = plan;
	  //DEBUG_TRACE("FFT plan created");printMemInfo();
  }
  else{
	  this->ftlgskern_plan = sensors->ftlgskern_plans.at(vdims);
	  //DEBUG_TRACE("FFT plan already exists : %d %d %d",mdims[0],mdims[1],dims_data3[3]);
  }

 // plan = this->d_ftlgskern->getPlan();
 // carmafftSafeCall(
 //     cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C , (int)dims_data3[3]));

  /*
   cudaExtent volumeSize = make_cudaExtent(npix,npix,nvalid);

   this->channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
   carmaSafeCall( cudaMalloc3DArray(&(this->d_spotarray), &(this->channelDesc), volumeSize) );

   // prepare 3d cppy
   this->copyParams.srcPtr = make_cudaPitchedPtr((void*)(this->d_lgskern->getData()), volumeSize.width*sizeof(float),
   volumeSize.width, volumeSize.height);
   copyParams.dstArray = this->d_spotarray;
   copyParams.extent   = volumeSize;
   copyParams.kind     = cudaMemcpyDeviceToDevice;
   */
  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;
}

sutra_lgs::~sutra_lgs() {
  current_context->set_activeDevice(device,1);
  delete this->d_doffaxis;
  delete this->d_azimuth;
  delete this->d_prof1d;
  delete this->d_profcum;
  delete this->d_prof2d;
  delete this->d_beam;
  delete this->d_ftbeam;
  //delete this->d_lgskern;
  //delete this->d_ftlgskern;

  //delete this->current_context;

  //carmaSafeCall(cudaFreeArray(this->d_spotarray));
}

int sutra_lgs::lgs_init(int nprof, float hg, float h0, float deltah,
    float pixsize, float *doffaxis, float *prof1d, float *profcum, float *beam,
    cuFloatComplex *ftbeam, float *azimuth) {
  current_context->set_activeDevice(device,1);
  this->nprof = nprof;
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;
  this->pixsize = pixsize;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = this->nvalid;
  this->d_doffaxis = new carma_obj<float>(current_context, dims_data1);
  this->d_doffaxis->host2device(doffaxis);

  dims_data1[1] = nprof + 1;
  this->d_prof1d = new carma_obj<float>(current_context, dims_data1);
  this->d_profcum = new carma_obj<float>(current_context, dims_data1);

  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->d_beam->host2device(beam);
  this->d_ftbeam->host2device(ftbeam);
  this->d_azimuth->host2device(azimuth);

  delete[] dims_data1;

  return EXIT_SUCCESS;
}

int sutra_lgs::load_prof(float *prof1d, float *profcum, float hg, float h0,
    float deltah) {
  current_context->set_activeDevice(device,1);
  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;

  return EXIT_SUCCESS;
}

int sutra_lgs::lgs_update(carma_device *device) {
  interp_prof(this->d_prof2d->getData(), this->d_prof1d->getData(),
      this->d_profcum->getData(), this->npix, this->d_doffaxis->getData(),
      this->hg, this->pixsize, this->h0, this->deltah, this->nprof,
      this->d_prof2d->getNbElem(), device);

  // convolution by beam
  // do fft on prof2d
  carma_fft(this->d_prof2d->getData(), this->d_prof2d->getData(), 1,
      *this->d_prof2d->getPlan());

  // mult by beamft
  times_ftbeam(this->d_prof2d->getData(), this->d_ftbeam->getData(), this->npix,
      this->d_prof2d->getNbElem(), device);

  // fft back
  carma_fft(this->d_prof2d->getData(), this->d_prof2d->getData(), -1,
      *this->d_prof2d->getPlan());

  return EXIT_SUCCESS;
}

int sutra_lgs::lgs_makespot(carma_device *device, int nin) {
  carmaSafeCall(
      cudaMemset(this->d_lgskern->getData(), 0,
          sizeof(float) * this->d_lgskern->getNbElem()));
  // build final image
  // get abs of real and roll
  cuFloatComplex *data = this->d_prof2d->getData();
  rollbeamexp(this->d_lgskern->getData(), &(data[nin]), this->d_beam->getData(),
      this->npix, this->npix*this->npix*this->nmaxhr/*this->d_lgskern->getNbElem()*/, device);

  // rotate image and fill kernels ft
  float *data2 = this->d_azimuth->getData();
  lgs_rotate(this->d_ftlgskern->getData(), this->d_lgskern->getData(),
      this->npix, this->npix, &(data2[nin / this->npix]), 0.0f,
      this->npix * this->npix * this->nmaxhr, device);

  //same with textures
  /*
   rotate3d(this->d_ftlgskern->getData(),this->copyParams,this->d_spotarray,this->channelDesc,
   this->npix,this->npix,this->d_azimuth->getData(),0.0f,this->npix*this->npix*this->nvalid,device);
   */

  //prepare for wfs code
  carma_fft(this->d_ftlgskern->getData(), this->d_ftlgskern->getData(), 1,
      *this->ftlgskern_plan);

  return EXIT_SUCCESS;
}

int sutra_lgs::load_kernels(float *lgskern, carma_device *device) {
  this->d_lgskern->host2device(lgskern);
  cfillrealp(this->d_ftlgskern->getData(), this->d_lgskern->getData(),
      /*this->d_ftlgskern->getNbElem()*/this->npix*this->npix*this->nmaxhr, device);
  carma_fft(this->d_ftlgskern->getData(), this->d_ftlgskern->getData(), 1,
      *this->ftlgskern_plan);

  return EXIT_SUCCESS;
}

