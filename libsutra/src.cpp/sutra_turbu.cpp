#include <sutra_turbu.h>
#include <algorithm>

//#include <culapackdevice.h>

/*  ____                             ____        __ _       _ _   _
 * / ___|  ___ _ __ ___  ___ _ __   |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __
 * \___ \ / __| '__/ _ \/ _ \ '_ \  | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \
 *  ___) | (__| | |  __/  __/ | | | | |_| |  __/  _| | | | | | |_| | (_) | | | |
 * |____/ \___|_|  \___|\___|_| |_| |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
 */

sutra_tscreen::sutra_tscreen(carma_context *context, long size, long size2,
    float r0, float altitude, float windspeed, float winddir, float deltax,
    float deltay, int device) {
  this->current_context = context;
  this->screen_size = size;
  this->amplitude = powf(r0, -5.0f / 6.0f);
  // ajust amplitude so that phase screens are generated in microns
  // r0 has been given @0.5µm
  this->amplitude *= 0.5;
  this->amplitude /= (2 * 3.14159265);

  this->altitude = altitude;
  this->windspeed = windspeed;
  this->winddir = winddir;
  this->accumx = 0.0f;
  this->accumy = 0.0f;
  this->deltax = deltax;
  this->deltay = deltay;
  this->device = device;

  this->norm_vk = 0;
  this->d_tscreen_c = 0;

  cout << "r0^-5/6 :" << this->amplitude << endl;

  this->d_tscreen = new sutra_phase(current_context, this->screen_size);
  this->channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,
      cudaChannelFormatKindFloat);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;
  this->d_tscreen_o = new carma_obj<float>(current_context, dims_data2);
  this->d_B = new carma_obj<float>(current_context, dims_data2);

  dims_data2[2] = size2;
  this->d_A = new carma_obj<float>(current_context, dims_data2);

  long *dims_data = new long[2];
  dims_data[0] = 1;
  dims_data[1] = size2;
  this->d_istencilx = new carma_obj<unsigned int>(current_context, dims_data);
  this->d_istencily = new carma_obj<unsigned int>(current_context, dims_data);
  this->d_z = new carma_obj<float>(current_context, dims_data);

  dims_data[1] = this->screen_size;
  this->d_noise = new carma_obj<float>(current_context, dims_data);
  this->d_ytmp = new carma_obj<float>(current_context, dims_data);

  this->vk_on = false;
}

sutra_tscreen::~sutra_tscreen() {
  delete this->d_tscreen;
  delete this->d_tscreen_o;
  delete this->d_A;
  delete this->d_B;
  delete this->d_istencilx;
  delete this->d_istencily;
  delete this->d_z;
  delete this->d_noise;
  delete this->d_ytmp;
  if (vk_on)
    delete d_tscreen_c;
  //delete current_context;
}

int sutra_tscreen::init_screen(float *h_A, float *h_B,
    unsigned int *h_istencilx, unsigned int *h_istencily, int seed) {
  // initial memcopies
  this->d_A->host2device(h_A);
  this->d_B->host2device(h_B);
  this->d_istencilx->host2device(h_istencilx);
  this->d_istencily->host2device(h_istencily);

  // init random noise
  if (this->d_noise->is_rng_init() == false)
    this->d_noise->init_prng_host(seed);
  this->d_noise->prng_host('N');

  return EXIT_SUCCESS;
}

int sutra_tscreen::init_vk(int seed, int pupd) {
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;
  this->d_tscreen_c = new carma_obj<cuFloatComplex>(current_context,
      dims_data2);
  cufftHandle *plan = this->d_tscreen_c->getPlan();
  cufftSafeCall(
      cufftPlan2d(plan, this->screen_size, this->screen_size, CUFFT_C2C));

  if (this->d_tscreen_o->is_rng_init() == false)
    this->d_tscreen_o->init_prng_host(seed);

  this->norm_vk = pow(pupd * pow(this->amplitude, 6.0f / 5.0f), 5.0f / 6.0f)
      * 0.5f / (2.0f * 3.14159);

  this->vk_on = true;

  return EXIT_SUCCESS;
}

int sutra_tscreen::generate_vk(float l0, int nalias) {
  this->d_tscreen_o->prng_host('N');

  cuFloatComplex *data = this->d_tscreen_c->getData();
  cutilSafeCall(
      cudaMemset(data, 0,
          this->screen_size * this->screen_size * sizeof(cuFloatComplex)));

  float k0 = (l0 == 0. ? 0.0f : this->screen_size / l0);
  int block_size = 8;

  float *data_o = this->d_tscreen_o->getData();
  gene_vonkarman(data, data_o, k0, nalias, this->screen_size, this->screen_size,
      block_size);

  roll<cuFloatComplex>(data, this->screen_size, this->screen_size, current_context->get_device(device));

  carma_fft(data, data, 1, *this->d_tscreen_c->getPlan());

  cgetrealp(this->d_tscreen->d_screen->getData(), this->d_tscreen_c->getData(),
      this->d_tscreen->d_screen->getNbElem(), current_context->get_device(device));

  norm_pscreen(data_o, this->d_tscreen->d_screen->getData(), this->screen_size,
      this->screen_size, this->norm_vk, this->current_context->get_device(device));

  this->d_tscreen->d_screen->copyFrom(data_o,
      this->d_tscreen->d_screen->getNbElem());

  return EXIT_SUCCESS;
}

int sutra_tscreen::extrude(int dir)
// dir =1 moving in x
    {
  int x0, Ncol, NC, N;
  NC = screen_size;

  if (dir == 1) { // adding a column to the left
    fillindx(this->d_z->getData(), this->d_tscreen->d_screen->getData(),
        (int *) this->d_istencilx->getData(), this->d_z->getNbElem(),
        current_context->get_device(device));
    x0 = this->screen_size - 1; //not in stencil
  } else {
    fillindx(this->d_z->getData(), this->d_tscreen->d_screen->getData(),
        (int *) this->d_istencily->getData(), this->d_z->getNbElem(),
        current_context->get_device(device));
    x0 = this->screen_size * (this->screen_size - 1);
  }

  addai<float>(this->d_z->getData(), this->d_tscreen->d_screen->getData(), x0, -1.0f,
      this->d_z->getNbElem(), current_context->get_device(device));

  this->d_ytmp->gemv('n', 1.0f, this->d_A, this->d_A->getDims(1), this->d_z, 1,
      0.0f, 1);

  this->d_noise->prng_host('N');

  this->d_ytmp->gemv('n', this->amplitude, this->d_B, this->d_B->getDims(1),
      this->d_noise, 1, 1.0f, 1);

  addai<float>(this->d_ytmp->getData(), this->d_tscreen->d_screen->getData(), x0, 1.0f,
      this->d_ytmp->getNbElem(), current_context->get_device(device));

  if (dir == 1) {
    x0 = 1;
    Ncol = this->screen_size - 1;
    N = this->screen_size * (this->screen_size - 1);
  } else {
    x0 = this->screen_size;
    Ncol = this->screen_size;
    N = this->screen_size * (this->screen_size - 1);
  }

  getarr2d(this->d_tscreen_o->getData(), this->d_tscreen->d_screen->getData(),
      x0, Ncol, NC, N, current_context->get_device(device));

  if (dir == 1)
    x0 = 0;
  else
    x0 = 0;

  fillarr2d(this->d_tscreen->d_screen->getData(), this->d_tscreen_o->getData(),
      x0, Ncol, NC, N, current_context->get_device(device));

  if (dir == 1) {
    x0 = this->screen_size - 1;
    Ncol = 1;
    N = this->screen_size;
  } else {
    x0 = this->screen_size * (this->screen_size - 1);
    Ncol = this->screen_size;
    N = this->screen_size;
  }

  fillarr2d(this->d_tscreen->d_screen->getData(), this->d_ytmp->getData(), x0,
      Ncol, NC, N, current_context->get_device(device));

  return EXIT_SUCCESS;
}

/*     _   _                         ____        __ _       _ _   _
 *    / \ | |_ _ __ ___   ___  ___  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __
 *   / _ \| __| '_ ` _ \ / _ \/ __| | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \
 *  / ___ \ |_| | | | | | (_) \__ \ | |_| |  __/  _| | | | | | |_| | (_) | | | |
 * /_/   \_\__|_| |_| |_|\___/|___/ |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
 */

sutra_atmos::sutra_atmos(carma_context *context, int nscreens, float *r0,
    long *size, long *size2, float *altitude, float *windspeed, float *winddir,
    float *deltax, float *deltay, float *pupil, int device) {
  this->nscreens = nscreens;
  //this->r0       = r0;
  this->current_context = context;
  this->r0 = 0;

  for (int i = 0; i < nscreens; i++) {
    d_screens.insert(
        pair<float, sutra_tscreen *>(altitude[i],
            new sutra_tscreen(context, size[i], size2[i], r0[i], altitude[i],
                windspeed[i], winddir[i], deltax[i], deltay[i], device)));
  }

  if (d_screens.find(0.0f) != d_screens.end()) {
    long msize = d_screens[0.0f]->screen_size;
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = msize;
    dims_data2[2] = msize;
    this->d_pupil = new carma_obj<float>(context, dims_data2);
    this->d_pupil->host2device(pupil);
  } else
    cout << "Could not find 0km altitude screen" << endl;

}

sutra_atmos::~sutra_atmos() {
  for (map<float, sutra_tscreen *>::iterator it = d_screens.begin();
      it != d_screens.end(); ++it) {
    delete it->second;
    it->second = 0;
  }
  d_screens.clear();
  //d_screens.erase(d_screens.begin(),d_screens.end());
}

int sutra_atmos::init_screen(float altitude, float *h_A, float *h_B,
    unsigned int *h_istencilx, unsigned int *h_istencily, int seed) {
  d_screens[altitude]->init_screen(h_A, h_B, h_istencilx, h_istencily, seed);

  return EXIT_SUCCESS;
}

int sutra_atmos::move_atmos() {

  map<float, sutra_tscreen *>::iterator p;
  p = this->d_screens.begin();

  while (p != this->d_screens.end()) {
    p->second->accumx += p->second->deltax;
    p->second->accumy += p->second->deltay;

    int deltax = (int) p->second->accumx;
    int deltay = (int) p->second->accumy;
    for (int cc = 0; cc < deltax; cc++)
      p->second->extrude(1);
    p->second->accumx -= deltax;
    for (int cc = 0; cc < deltay; cc++)
      p->second->extrude(0);
    p->second->accumy -= deltay;
    p++;
  }
  return EXIT_SUCCESS;
}

/*
 void checkCulaStatus(culaStatus status)
 {
 if(!status)
 return;

 if(status == culaArgumentError)
 printf("Invalid value for parameter %d\n", culaGetErrorInfo());
 else if(status == culaDataError)
 printf("Data error (%d)\n", culaGetErrorInfo());
 else if(status == culaBlasError)
 printf("Blas error (%d)\n", culaGetErrorInfo());
 else if(status == culaRuntimeError)
 printf("Runtime error (%d)\n", culaGetErrorInfo());
 else
 printf("%s\n", culaGetStatusString(status));

 culaShutdown();
 exit(EXIT_FAILURE);
 }

 int _initAB(long size, void *h_zz, void *h_xz, void *h_txz, void *h_xx, float *d_A, float *d_B,float *tmpout)
 {

 culaStatus status;
 status = culaInitialize();
 checkCulaStatus(status);

 carma_matmult *mmLar1;
 carma_matmult *mmLar2;
 carma_matmult *mmMed1;
 carma_matmult *mmMed2;
 carma_matmult *mmSma;

 long size_data[6];
 size_data[0] = size_data[1] = size_data[2] = 2*size;
 size_data[3] = size_data[4] = size_data[5] = 2*size;
 mmLar1 = new carma_matmult(size_data,'F');
 mmLar2 = new carma_matmult(size_data,'F');

 size_data[0] = size; size_data[4] = size;
 mmMed1 = new carma_matmult(size_data,'F');

 size_data[3] = size_data[5] = size;
 mmMed2 = new carma_matmult(size_data,'F');

 size_data[1] = size_data[2] = size;
 mmSma = new carma_matmult(size_data,'F');

 char jobu = 'A';
 char jobvt = 'A';

 float* S = NULL;
 float* UU = NULL;

 cudaMalloc((void**)&S, 2*size*sizeof(float));
 cudaMalloc((void**)&UU, 2*size*2*size*sizeof(float));

 checkCulaStatus(status);

 //s = s1 = SVdec(zz,uuu,vt);
 mmLar1->host2device(0,h_zz);  

 status = culaDeviceSgesvd(jobu,jobvt,2*size,2*size,(float *)mmLar1->data[0],2*size, S, 
 UU, 2*size,(float *)mmLar1->data[0],2*size);

 checkCulaStatus(status);
 cudaFree(UU);

 //in this case h_zz is a symmetric matrix. so u = v. 
 //apparently in this case the lapack lib returns 2 different versions of u and v. 
 //indeed in the svd problem, the solution is not unique. here we can use both, 
 // we chose vt which seems to give the best results.
 // however : the single precision is not enough to get good estimation of bbt 
 // (see below). so this code is not working. need to buy the premium version of 
 // cula to get double precision on the svd


 cutilSafeCall(cudaMemcpy(mmLar2->data[1], mmLar1->data[0], 2*size*2*size*sizeof(float), 
 cudaMemcpyDeviceToDevice));
 // s1(0)=1;s1 = 1./s1;s1(0)=0; // the null eignevalue is not inverted
 cutilSafeCall(cudaMemset((float *)mmLar1->data[1],0,2*size*2*size * sizeof(float)));

 getinvCU((float *)mmLar1->data[1],S,2*size);

 cudaFree(S);

 //zz_1 =  (uuu*s1(-,))(,+) * vt(+,);   // inversion, with the null eigenvalue left to 0
 mmLar1->compute('t','n',1.0f,0.0f,mmLar2->data[0]);  

 mmLar2->compute('n','n',1.0f,0.0f,mmMed1->data[1]);  

 //A = xz(,+) * zz_1(+,);
 mmMed1->host2device(0,h_xz);  

 mmMed1->compute('n','n',1.0f,0.0f,mmMed2->data[0]);  


 cutilSafeCall(cudaMemcpy(d_A, mmMed2->data[0], size*2*size*sizeof(float), 
 cudaMemcpyDeviceToDevice));

 //bbt = xx - A(,+)*xz(,+);
 mmMed2->host2device(2,h_xx);  
 mmMed2->host2device(1,h_txz);  

 mmMed2->compute('n','n',-1.0f,1.0f,mmMed2->data[2]);  

 cutilSafeCall(cudaMemcpy(tmpout, mmMed2->data[2], size*size*sizeof(float), 
 cudaMemcpyDeviceToHost));
 float* VT = NULL;
 cudaMalloc((void**)&S, size*sizeof(float));
 cudaMalloc((void**)&VT, size*size*sizeof(float));

 //l = SVdec(bbt,uu);
 status = culaDeviceSgesvd(jobu, jobvt, size, size,(float *)mmMed2->data[2] , size, S, 
 (float *)mmSma->data[0], size, VT, size);

 checkCulaStatus(status);


 //fprintf(stderr,"%f %f %f \n",S[0],S[1],S[2]);

 //B = uu*sqrt(l)(-,);
 getsqrtCU((float *)mmSma->data[1],S,size);

 mmSma->compute('n','n',1.0f,0.0f,mmSma->data[2]);  

 cutilSafeCall(cudaMemcpy(d_B, mmSma->data[2], size*size*sizeof(float), 
 cudaMemcpyDeviceToDevice));

 cudaFree(VT);
 cudaFree(S);

 cublasStatus bstatus;
 int i;
 for( i=0; i<3; i++) {
 bstatus = cublasFree(mmLar1->data[i]);
 if (bstatus != CUBLAS_STATUS_SUCCESS) {
 fprintf (stderr, "!!!! Error mmLar1 matrix %d\n", i);
 }
 }
 for( i=0; i<3; i++) {
 bstatus = cublasFree(mmLar2->data[i]);
 if (bstatus != CUBLAS_STATUS_SUCCESS) {
 fprintf (stderr, "!!!! Error mmLar2 matrix %d\n", i);
 }
 }
 for( i=0; i<3; i++) {
 bstatus = cublasFree(mmMed1->data[i]);
 if (bstatus != CUBLAS_STATUS_SUCCESS) {
 fprintf (stderr, "!!!! Error mmMed1 matrix %d\n", i);
 }
 }
 for( i=0; i<3; i++) {
 bstatus = cublasFree(mmMed2->data[i]);
 if (bstatus != CUBLAS_STATUS_SUCCESS) {
 fprintf (stderr, "!!!! Error mmMed2 matrix %d\n", i);
 }
 }
 for( i=0; i<3; i++) {
 bstatus = cublasFree(mmSma->data[i]);
 if (bstatus != CUBLAS_STATUS_SUCCESS) {
 fprintf (stderr, "!!!! Error mmSma matrix %d\n", i);
 }
 }
 //delete mmLar1;
 //delete mmLar2;
 //delete mmMed1;
 //delete mmMed2;
 //delete mmSma;

 culaShutdown();
 
 return EXIT_SUCCESS;
 }

 */
