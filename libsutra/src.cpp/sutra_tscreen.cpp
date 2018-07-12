#include <sutra_tscreen.h>

//   ██████╗ ██████╗ ███╗   ██╗███████╗████████╗██████╗ ██╗   ██╗
//  ██╔════╝██╔═══██╗████╗  ██║██╔════╝╚══██╔══╝██╔══██╗██║   ██║
//  ██║     ██║   ██║██╔██╗ ██║███████╗   ██║   ██████╔╝██║   ██║
//  ██║     ██║   ██║██║╚██╗██║╚════██║   ██║   ██╔══██╗██║   ██║
//  ╚██████╗╚██████╔╝██║ ╚████║███████║   ██║   ██║  ██║╚██████╔╝██╗
//   ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═╝
//

sutra_tscreen::sutra_tscreen(carma_context *context, long size,
                             long stencilSize, float r0, float altitude,
                             float windspeed, float winddir, float deltax,
                             float deltay, int device) {
  this->current_context = context;
  this->screen_size = size;
  this->r0 = r0;
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
  this->current_context->set_activeDevice(device, 1);

  this->norm_vk = 0;
  this->d_tscreen_c = 0;

  std::cout << "r0^-5/6 :" << this->amplitude << std::endl;

  this->d_tscreen = new sutra_phase(current_context, this->screen_size);
  this->channelDesc =
      cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;
  this->d_tscreen_o = new carma_obj<float>(current_context, dims_data2);
  this->d_B = new carma_obj<float>(current_context, dims_data2);

  dims_data2[2] = stencilSize;
  this->d_A = new carma_obj<float>(current_context, dims_data2);

  long dims_data[2];
  dims_data[0] = 1;
  dims_data[1] = stencilSize;
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
  if (vk_on) delete d_tscreen_c;
  // delete current_context;
}

//  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
//  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
//  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
//  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
//  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
//  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
//

int sutra_tscreen::init_screen(float *h_A, float *h_B,
                               unsigned int *h_istencilx,
                               unsigned int *h_istencily, int seed) {
  this->current_context->set_activeDevice(device, 1);
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

int sutra_tscreen::refresh_screen() {
  this->current_context->set_activeDevice(device, 1);
  this->d_tscreen->d_screen->reset();
  this->accumx = 0.0f;
  this->accumy = 0.0f;

  for (int i = 0; i < 2 * this->screen_size; i++) {
    if (this->deltax > 0) {
      this->extrude(1);
    } else {
      this->extrude(-1);
    }
  }
  return EXIT_SUCCESS;
}

int sutra_tscreen::set_seed(int seed) {
  this->current_context->set_activeDevice(device, 1);
  this->d_noise->init_prng_host(seed);
  this->d_noise->prng_host('N');

  return EXIT_SUCCESS;
}

int sutra_tscreen::extrude(int dir) {
  // dir =1 moving in x

  this->current_context->set_activeDevice(device, 1);
  int x0, Ncol, NC, N;
  NC = screen_size;

  if (dir == 1 || dir == -1) {  // adding a column to the left
    fillindx(this->d_z->getData(), this->d_tscreen->d_screen->getData(),
             (int *)this->d_istencilx->getData(), this->d_z->getNbElem(),
             current_context->get_device(device));
    if (dir == 1)
      x0 = this->screen_size - 1;  // not in stencil
    else
      x0 = this->screen_size * (this->screen_size - 1);
  } else {
    fillindx(this->d_z->getData(), this->d_tscreen->d_screen->getData(),
             (int *)this->d_istencily->getData(), this->d_z->getNbElem(),
             current_context->get_device(device));
    if (dir == 2)
      x0 = this->screen_size * (this->screen_size - 1);
    else
      x0 = this->screen_size - 1;
  }

  addai<float>(this->d_z->getData(), this->d_tscreen->d_screen->getData(), x0,
               -1.0f, this->d_z->getNbElem(),
               current_context->get_device(device));

  this->d_ytmp->gemv('n', 1.0f, this->d_A, this->d_A->getDims(1), this->d_z, 1,
                     0.0f, 1);

  this->d_noise->prng_host('N');

  this->d_ytmp->gemv('n', this->amplitude, this->d_B, this->d_B->getDims(1),
                     this->d_noise, 1, 1.0f, 1);

  addai<float>(this->d_ytmp->getData(), this->d_tscreen->d_screen->getData(),
               x0, 1.0f, this->d_ytmp->getNbElem(),
               current_context->get_device(device));

  if (dir == 1 || dir == -1) {
    if (dir == 1)
      x0 = 1;
    else
      x0 = 0;
    Ncol = this->screen_size - 1;
    N = this->screen_size * (this->screen_size - 1);
  } else {
    if (dir == 2) x0 = this->screen_size;
    if (dir == -2) x0 = 0;
    Ncol = this->screen_size;
    N = this->screen_size * (this->screen_size - 1);
  }

  getarr2d(this->d_tscreen_o->getData(), this->d_tscreen->d_screen->getData(),
           x0, Ncol, NC, N, current_context->get_device(device));

  if (dir > 0) x0 = 0;
  if (dir == -1) x0 = 1;
  if (dir == -2) x0 = this->screen_size;

  fillarr2d(this->d_tscreen->d_screen->getData(), this->d_tscreen_o->getData(),
            x0, Ncol, NC, N, current_context->get_device(device));

  if (dir == 1 || dir == -1) {
    if (dir == 1)
      x0 = this->screen_size - 1;
    else
      x0 = 0;
    Ncol = 1;
    N = this->screen_size;
  } else {
    if (dir == 2) x0 = this->screen_size * (this->screen_size - 1);
    if (dir == -2) x0 = 0;
    Ncol = this->screen_size;
    N = this->screen_size;
  }

  fillarr2d(this->d_tscreen->d_screen->getData(), this->d_ytmp->getData(), x0,
            Ncol, NC, N, dir, current_context->get_device(device));

  return EXIT_SUCCESS;
}

//  ██╗   ██╗███╗   ██╗██╗   ██╗███████╗███████╗██████╗
//  ██║   ██║████╗  ██║██║   ██║██╔════╝██╔════╝██╔══██╗
//  ██║   ██║██╔██╗ ██║██║   ██║███████╗█████╗  ██║  ██║
//  ██║   ██║██║╚██╗██║██║   ██║╚════██║██╔══╝  ██║  ██║
//  ╚██████╔╝██║ ╚████║╚██████╔╝███████║███████╗██████╔╝
//   ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚══════╝╚═════╝
//

int sutra_tscreen::init_vk(int seed, int pupd) {
  this->current_context->set_activeDevice(device, 1);
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;
  this->d_tscreen_c =
      new carma_obj<cuFloatComplex>(current_context, dims_data2);
  cufftHandle *plan = this->d_tscreen_c->getPlan();
  carmafftSafeCall(
      cufftPlan2d(plan, this->screen_size, this->screen_size, CUFFT_C2C));

  if (this->d_tscreen_o->is_rng_init() == false)
    this->d_tscreen_o->init_prng_host(seed);

  this->norm_vk = pow(pupd * pow(this->amplitude, 6.0f / 5.0f), 5.0f / 6.0f) *
                  0.5f / (2.0f * 3.14159);

  this->vk_on = true;

  delete[] dims_data2;
  return EXIT_SUCCESS;
}

int sutra_tscreen::generate_vk(float l0, int nalias) {
  this->current_context->set_activeDevice(device, 1);
  this->d_tscreen_o->prng_host('N');

  cuFloatComplex *data = this->d_tscreen_c->getData();
  carmaSafeCall(cudaMemset(
      data, 0, this->screen_size * this->screen_size * sizeof(cuFloatComplex)));

  float k0 = (l0 == 0. ? 0.0f : this->screen_size / l0);
  int block_size = 8;

  float *data_o = this->d_tscreen_o->getData();
  gene_vonkarman(data, data_o, k0, nalias, this->screen_size, this->screen_size,
                 block_size);

  roll<cuFloatComplex>(data, this->screen_size, this->screen_size,
                       current_context->get_device(device));

  carma_fft(data, data, 1, *this->d_tscreen_c->getPlan());

  cgetrealp(this->d_tscreen->d_screen->getData(), this->d_tscreen_c->getData(),
            this->d_tscreen->d_screen->getNbElem(),
            current_context->get_device(device));

  norm_pscreen(data_o, this->d_tscreen->d_screen->getData(), this->screen_size,
               this->screen_size, this->norm_vk,
               this->current_context->get_device(device));

  this->d_tscreen->d_screen->copyFrom(data_o,
                                      this->d_tscreen->d_screen->getNbElem());

  return EXIT_SUCCESS;
}
