#include <sutra_atmos.h>
#include <algorithm>

//   ██████╗ ██████╗ ███╗   ██╗███████╗████████╗██████╗ ██╗   ██╗
//  ██╔════╝██╔═══██╗████╗  ██║██╔════╝╚══██╔══╝██╔══██╗██║   ██║
//  ██║     ██║   ██║██╔██╗ ██║███████╗   ██║   ██████╔╝██║   ██║
//  ██║     ██║   ██║██║╚██╗██║╚════██║   ██║   ██╔══██╗██║   ██║
//  ╚██████╗╚██████╔╝██║ ╚████║███████║   ██║   ██║  ██║╚██████╔╝██╗
//   ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═╝
//
sutra_atmos::sutra_atmos(carma_context *context, int nscreens, float global_r0,
                         float *r0, long *size, long *stencilSize,
                         float *altitude, float *windspeed, float *winddir,
                         float *deltax, float *deltay, int device) {
  this->nscreens = nscreens;
  // this->r0       = r0;
  this->current_context = context;
  this->r0 = global_r0;

  for (int i = 0; i < nscreens; i++) {
    d_screens.insert(pair<float, sutra_tscreen *>(
        altitude[i], new sutra_tscreen(context, size[i], stencilSize[i], r0[i],
                                       altitude[i], windspeed[i], winddir[i],
                                       deltax[i], deltay[i], device)));
  }
}

sutra_atmos::~sutra_atmos() {
  for (map<float, sutra_tscreen *>::iterator it = d_screens.begin();
       it != d_screens.end(); ++it) {
    delete it->second;
    it->second = 0;
  }
  d_screens.clear();
  // d_screens.erase(d_screens.begin(),d_screens.end());
}

//  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
//  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
//  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
//  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
//  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
//  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
//

int sutra_atmos::init_screen(float altitude, float *h_A, float *h_B,
                             unsigned int *h_istencilx,
                             unsigned int *h_istencily, int seed) {
  if (this->d_screens.count(altitude)) {
    d_screens[altitude]->init_screen(h_A, h_B, h_istencilx, h_istencily, seed);
    d_screens[altitude]->refresh_screen();
  } else
    DEBUG_TRACE("Screen not found at this altitude");

  return EXIT_SUCCESS;
}

int sutra_atmos::refresh_screen(float altitude) {
  if (this->d_screens.count(altitude))
    this->d_screens[altitude]->refresh_screen();
  else
    DEBUG_TRACE("Screen not found at this altitude");
  return EXIT_SUCCESS;
}

int sutra_atmos::add_screen(float alt, long size, long stencilSize,
                            float r0_thislayer, float windspeed, float winddir,
                            float deltax, float deltay, int device) {
  if (d_screens.find(alt) != d_screens.end()) {
    std::cout << "There is already a screen at this altitude" << std::endl;
    std::cout << "No screen created" << std::endl;
    return EXIT_FAILURE;
  }

  sutra_tscreen *screen =
      new sutra_tscreen(current_context, size, stencilSize, r0_thislayer, alt,
                        windspeed, winddir, deltax, deltay, device);
  this->r0 = powf(powf(r0, -5.0f / 3.0f) + powf(r0_thislayer, -5.0f / 3.0f),
                  -3.0f / 5.0f);
  d_screens.insert(pair<float, sutra_tscreen *>(alt, screen));
  nscreens++;
  return EXIT_SUCCESS;
}

int sutra_atmos::del_screen(const float alt) {
  if (d_screens.find(alt) == d_screens.end()) {
    std::cout << "There is no screen at this altitude" << std::endl;
    std::cout << "No screen erased" << std::endl;
    return EXIT_FAILURE;
  }
  nscreens--;
  this->r0 =
      powf(powf(r0, -5.0f / 3.0f) - powf(d_screens[alt]->r0, -5.0f / 3.0f),
           -3.0f / 5.0f);
  delete d_screens[alt];
  d_screens.erase(alt);
  return EXIT_SUCCESS;
}

int sutra_atmos::move_atmos() {
  map<float, sutra_tscreen *>::iterator p;
  p = this->d_screens.begin();

  while (p != this->d_screens.end()) {
    p->second->accumx += p->second->deltax;
    p->second->accumy += p->second->deltay;

    int deltax = (int)p->second->accumx;
    int deltay = (int)p->second->accumy;
    int cx = deltax > 0 ? 1 : -1;
    int cy = deltay > 0 ? 1 : -1;
    for (int cc = 0; cc < cx * deltax; cc++) p->second->extrude(1 * cx);
    p->second->accumx -= deltax;
    for (int cc = 0; cc < cy * deltay; cc++) p->second->extrude(2 * cy);
    p->second->accumy -= deltay;
    p++;
  }
  return EXIT_SUCCESS;
}

int sutra_atmos::set_global_r0(float r0) {
  // this->amplitude = powf(r0, -5.0f / 6.0f)
  float scaling = powf(r0 / this->r0, -5.0f / 6.0f);
  for (map<float, sutra_tscreen *>::iterator it = d_screens.begin();
       it != d_screens.end(); ++it) {
    it->second->r0 *= this->r0 / r0;
    it->second->amplitude *= scaling;
  }
  this->r0 = r0;
  return EXIT_SUCCESS;
}

int sutra_atmos::set_seed(float alt, float seed) {
  if (this->d_screens.count(alt))
    this->d_screens[alt]->set_seed(seed);
  else
    DEBUG_TRACE("Screen not found at this altitude");

  return EXIT_SUCCESS;
}
