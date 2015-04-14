/*! @file yoga_ao.cpp

 @brief Code for functions callable by the Yorick interpreter

 see yoga_ao_api.h for a complete description

 Copyright (c) 2011, Damien Gratadour & Arnaud Sevin

 This file is part of YoGA, the Yorick with GPU Acceleration plugin

 This program is free software; you can redistribute it and/or  modify it
 under the terms of the GNU General Public License  as  published  by the
 Free Software Foundation; either version 3 of the License,  or  (at your
 option) any later version.
 
 This program is distributed in the hope  that  it  will  be  useful, but
 WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sutra_turbu.h>
#include <sutra_target.h>
#include <sutra_phase.h>
#include <sutra_wfs_geom.h>
#include <sutra_wfs_pyr_pyr4.h>
#include <sutra_wfs_pyr_roof.h>
#include <sutra_wfs_sh.h>
#include <sutra_acquisim.h>
#include <sutra_rtc.h>
#include <sutra_dm.h>
#include <sutra_telemetry.h>
#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_cog.h>
#include <sutra_centroider_corr.h>
#include <sutra_centroider_pyr.h>
#include <sutra_centroider_roof.h>
#include <sutra_centroider_tcog.h>
#include <sutra_centroider_wcog.h>
#include <carma.h>
#include <sstream>
#include <iomanip>
#include <sutra_aotemplate.h>
#include <yoga_api.h>
#include <yoga_ao_api.h>

/**
 * @typedef Yorick API atmos userobj
 */
static y_userobj_t yAtmos = { const_cast<char*>("yAtmos Object"), &atmos_free,
                              &atmos_print, 0, 0, 0 };

/**
 * @typedef Yorick API tscreen userobj
 */
static y_userobj_t yTscreen = { const_cast<char*>("yTscreen Object"),
                                &tscreen_free, 0, 0, 0, 0 };

static y_userobj_t ySource = { const_cast<char*>("ySource Object"),
                               &source_free, 0, 0, 0, 0 };

static y_userobj_t yTarget = { const_cast<char*>("yTarget Object"),
                               &target_free, &target_print, 0, 0, 0 };

static y_userobj_t yPhase = { const_cast<char*>("yPhase Object"), &phase_free,
                              &phase_print, &phase_eval, 0, 0 };
/*
 static y_userobj_t yWfs = { const_cast<char*>("yWfs Object"), &wfs_free,
 &wfs_print, 0, 0, 0 };
 */
static y_userobj_t yTelemetry = { const_cast<char*>("yTelemetry Object"),
                                  &telemetry_free, &telemetry_print, 0, 0, 0 };

static y_userobj_t ySensors = { const_cast<char*>("ySensors Object"),
                                &sensors_free, &sensors_print, 0, 0, 0 };

static y_userobj_t yAcquisim = { const_cast<char*>("yAcquisim Object"),
                                 &acquisim_free, &acquisim_print, 0, 0, 0 };

static y_userobj_t yDMs = { const_cast<char*>("yDMs Object"), &dms_free,
                            &dms_print, 0, 0, 0 };

static y_userobj_t yRTC = { const_cast<char*>("yRTC Object"), &rtc_free,
                            &rtc_print, 0, 0, 0 };

static y_userobj_t yAotemplate =
    { const_cast<char*>("yAotemplate Object"), &aotemplate_free,
      &aotemplate_print, 0, 0, 0 };

/*
 *         _   _ _ _ _   _
 *   _   _| |_(_) (_) |_(_) ___  ___
 *  | | | | __| | | | __| |/ _ \/ __|
 *  | |_| | |_| | | | |_| |  __/\__ \
 *   \__,_|\__|_|_|_|\__|_|\___||___/
 *
 */

extern "C" {

  /*
   *      _   _
   *     / \ | |_ _ __ ___   ___  ___
   *    / _ \| __| '_ ` _ \ / _ \/ __|
   *   / ___ \ |_| | | | | | (_) \__ \
   *  /_/   \_\__|_| |_| |_|\___/|___/
   *
   */

  atmos_struct*
  yoga_ao_getyAtmos(int argc, int pos) {
    return (atmos_struct *) yget_obj(argc - pos, &yAtmos);
  }

  void atmos_free(void *obj) {
    atmos_struct *handler = (atmos_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_atmos *atmos_obj_handler = (sutra_atmos *) (handler->sutra_atmos);
      delete atmos_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void atmos_print(void *obj) {
    atmos_struct *handler = (atmos_struct *) obj;
    sutra_atmos *atmos_obj_handler = (sutra_atmos *) (handler->sutra_atmos);
    map<float, sutra_tscreen *>::iterator p;
    p = atmos_obj_handler->d_screens.begin();
    cout << "YoGA Atmos Object : " << endl;
    cout << "Contains " << atmos_obj_handler->nscreens
         << " turbulent screen(s) : " << endl;
    int i = 0;
    sutra_tscreen *tmp;
    cout << "Screen # " << " | " << "alt.(m)" << " | " << "speed (m/s)" << " | "
         << "dir.(deg)" << " | " << "r0 (pix)" << " | " << "deltax" << " | "
         << "deltay" << endl;
    while (p != atmos_obj_handler->d_screens.end()) {
      tmp = p->second;
      cout.precision(5);
      cout << setw(9) << i + 1 << " | " << setw(7) << tmp->altitude << " | "
           << setw(11) << tmp->windspeed << " | " << setw(9) << tmp->winddir
           << " | " << setw(8) << pow(tmp->amplitude, -6.0f / 5.0f) << " | "
           << setw(6) << tmp->deltax << " | " << setw(6) << tmp->deltay << endl;
      p++;
      i++;
    }
  }

  void Y_yoga_atmos(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();
      int odevice = activeDevice;
      int nscreens = ygets_i(argc - 1);
      float *r0 = ygeta_f(argc - 2, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens r0");

      long *size;
      size = ygeta_l(argc - 3, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens size");

      long *size2;
      size2 = ygeta_l(argc - 4, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens size");

      float *alt;
      alt = ygeta_f(argc - 5, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens alt");

      float *wspeed;
      wspeed = ygeta_f(argc - 6, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens speed");

      float *wdir;
      wdir = ygeta_f(argc - 7, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens dir");

      float *deltax;
      deltax = ygeta_f(argc - 8, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens deltax");

      float *deltay;
      deltay = ygeta_f(argc - 9, &ntot, dims);
      if (ntot != nscreens)
        y_error("wrong dimension for screens deltay");

      //if (ntot != nscreens) y_error("wrong dimension for screens speed");

      if (argc > 9)
        odevice = ygets_i(argc - 10);
      activeDevice = context_handle->set_activeDevice(odevice,1);

      atmos_struct *handle = (atmos_struct *) ypush_obj(&yAtmos,
                                                        sizeof(atmos_struct));
      handle->device = odevice;

      handle->sutra_atmos = new sutra_atmos(context_handle, nscreens,
                                            (float *) r0, (long *) size,
                                            (long *) size2, (float *) alt,
                                            (float *) wspeed, (float *) wdir,
                                            (float *) deltax, (float *) deltay,
                                            odevice);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_atmos construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  /*
   *     _   ____
   *    | |_/ ___|  ___ _ __ ___  ___ _ __
   *    | __\___ \ / __| '__/ _ \/ _ \ '_ \
   *    | |_ ___) | (__| | |  __/  __/ | | |
   *     \__|____/ \___|_|  \___|\___|_| |_|
   *
   */

  tscreen_struct*
  yoga_ao_getyTscreen(int argc, int pos) {
    return (tscreen_struct *) yget_obj(argc - pos, &yTscreen);
  }

  void tscreen_free(void *obj) {
    tscreen_struct *handler = (tscreen_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_tscreen *tscreen_obj_handler =
          (sutra_tscreen *) (handler->sutra_tscreen);
      delete tscreen_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void tscreen_print(void *obj) {
    tscreen_struct *handler = (tscreen_struct *) obj;
    sutra_tscreen *tscreen_obj_handler =
        (sutra_tscreen *) (handler->sutra_tscreen);
    cout << "Yoga tScreen Object : " << endl;
    cout << "Screen @ alt. " << tscreen_obj_handler->altitude << " m |"
         << " wind speed : " << tscreen_obj_handler->windspeed << " m/s |"
         << " wind dir. : " << tscreen_obj_handler->winddir << " deg |"
         << "r_0 : " << pow(tscreen_obj_handler->amplitude, -6. / 5.) << " m"
         << endl;
  }

  void Y_init_tscreen(int argc)
// here we init the atmos structure
// args are : the sutra_atmos, nscreen A, B, istencilx,istencily and seed
      {
    long ntot;
    long dims[Y_DIMSIZE];
    atmos_struct *handler = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
    sutra_atmos *atmos_obj_handler = (sutra_atmos *) (handler->sutra_atmos);

    carma_context *context_handle = _getCurrentContext();
    handler->device = context_handle->get_activeDevice();

    float altitude = ygets_f(argc - 2);

    float *a = ygeta_f(argc - 3, &ntot, dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_A->getNbElem())
      y_error("wrong size for A");
    float *b = ygeta_f(argc - 4, &ntot, dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_B->getNbElem())
      y_error("wrong size for B");
    unsigned int *istencilx = (unsigned int *) ygeta_i(argc - 5, &ntot, dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_istencilx->getNbElem())
      y_error("wrong size for istencilx");
    unsigned int *istencily = (unsigned int *) ygeta_i(argc - 6, &ntot, dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_istencily->getNbElem())
      y_error("wrong size for istencily");

    int seed = ygets_i(argc - 7);

    atmos_obj_handler->init_screen(altitude, a, b, istencilx, istencily, seed);
  }

  void Y_get_tscreen(int argc) {

    if (yarg_subroutine())
      y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDeviceForCpy(handle->device,1);
      float alt = ygets_f(argc - 2);
      sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
      caObjS *carma_obj_handler =
          (caObjS *) (atmos_handler->d_screens[alt]->d_tscreen->d_screen);
      float *data = ypush_f((long*) carma_obj_handler->getDims());
      carma_obj_handler->device2host(data);
    }
  }

  void Y_set_tscreen(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    if (!yarg_subroutine())
      y_error("can only be called as a subroutine");
    else {
      atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDeviceForCpy(handle->device,1);
      float alt = ygets_f(argc - 2);
      float *data = ygeta_f(argc - 3, &ntot, dims);
      sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
      caObjS *carma_obj_handler =
          (caObjS *) (atmos_handler->d_screens[alt]->d_tscreen->d_screen);
      carma_obj_handler->host2device(data);
    }
  }

  void Y_extrude_tscreen(int argc) {

    if (yarg_subroutine()) {
      atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDeviceForCpy(handle->device,1);
      float alt = ygets_f(argc - 2);
      sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
      sutra_tscreen *tscreen_handler =
          (sutra_tscreen *) (atmos_handler->d_screens.find(alt)->second);
      int dir;
      if (argc > 2)
        dir = ygets_i(argc - 3);
      else
        dir = 1;
      tscreen_handler->extrude(dir);
    } else
      y_error("can only be called as a subroutine");
  }

  void Y_get_tscreen_config(int argc) {

    if (yarg_subroutine())
      y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDeviceForCpy(handle->device,1);
      float alt = ygets_f(argc - 2);
      sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
      sutra_tscreen *tscreen_handler =
          (sutra_tscreen *) (atmos_handler->d_screens[alt]);
      char *type_data = ygets_q(argc - 3);
      if (strcmp(type_data, "A") == 0) {
        float *data = ypush_f((long*) tscreen_handler->d_A->getDims());
        tscreen_handler->d_A->device2host(data);
      } else if (strcmp(type_data, "B") == 0) {
        float *data = ypush_f((long*) tscreen_handler->d_B->getDims());
        tscreen_handler->d_B->device2host(data);
      } else if (strcmp(type_data, "istx") == 0) {
        int *data = ypush_i((long*) tscreen_handler->d_istencilx->getDims());
        tscreen_handler->d_istencilx->device2host((unsigned int *) data);
      } else if (strcmp(type_data, "isty") == 0) {
        int *data = ypush_i((long*) tscreen_handler->d_istencily->getDims());
        tscreen_handler->d_istencily->device2host((unsigned int *) data);
      } else if (strcmp(type_data, "values") == 0) {
        int *data = ypush_i(
            (long*) tscreen_handler->d_tscreen->d_screen->getDims());
        carmaSafeCall(
            cudaMemcpy(
                (void ** )&data,
                tscreen_handler->d_tscreen->d_screen->getValues(),
                sizeof(unsigned int) * tscreen_handler->d_tscreen->d_screen->getNbElem(),
                cudaMemcpyDeviceToHost));
      }
    }
  }

  void Y_get_tscreen_update(int argc) {

    if (yarg_subroutine())
      y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDeviceForCpy(handle->device,1);
      float alt = ygets_f(argc - 2);
      sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
      sutra_tscreen *tscreen_handler =
          (sutra_tscreen *) (atmos_handler->d_screens[alt]);
      float *data = ypush_f((long*) tscreen_handler->d_ytmp->getDims());
      tscreen_handler->d_ytmp->device2host(data);
    }
  }

  void Y_tscreen_initvk(int argc) {

    atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handle->device,1);
    float alt = ygets_f(argc - 2);
    int pupd = ygets_l(argc - 3);
    long seed = 1234;
    if (argc > 3)
      seed = ygets_l(argc - 4);
    sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
    sutra_tscreen *tscreen_handler =
        (sutra_tscreen *) (atmos_handler->d_screens[alt]);
    tscreen_handler->init_vk(seed, pupd);
  }

  void Y_tscreen_genevk(int argc) {

    atmos_struct *handle = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handle->device,1);
    float alt = ygets_f(argc - 2);
    float l0 = 0.0f;
    if (argc > 2)
      l0 = ygets_f(argc - 3);
    long nalias = 0;
    if (argc > 3)
      l0 = ygets_l(argc - 4);
    sutra_atmos *atmos_handler = (sutra_atmos *) handle->sutra_atmos;
    sutra_tscreen *tscreen_handler =
        (sutra_tscreen *) (atmos_handler->d_screens[alt]);
    tscreen_handler->generate_vk(l0, nalias);
  }

  /*
   *     ____
   *    / ___|  ___  _   _ _ __ ___ ___
   *    \___ \ / _ \| | | | '__/ __/ _ \
   *     ___) | (_) | |_| | | | (_|  __/
   *    |____/ \___/ \__,_|_|  \___\___|
   *
   */

  source_struct*
  yoga_ao_getySource(int argc, int pos) {
    return (source_struct *) yget_obj(argc - pos, &ySource);
  }

  void source_free(void *obj) {
    source_struct *handler = (source_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_source *source_obj_handler =
          (sutra_source *) (handler->sutra_source);
      delete source_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void source_print(void *obj) {
    source_struct *handler = (source_struct *) obj;
    sutra_source *source_obj_handler = (sutra_source *) (handler->sutra_source);
    cout << "Yoga Source Object : " << endl;
    cout << "Source @ pos. " << source_obj_handler->tposx << "\" : "
         << source_obj_handler->tposy << " \"  |" << " mag : "
         << source_obj_handler->mag << " |" << " lambda : "
         << source_obj_handler->lambda << " micron " << endl;
  }

  /*
   *     _____                    _
   *    |_   _|_ _ _ __ __ _  ___| |_
   *      | |/ _` | '__/ _` |/ _ \ __|
   *      | | (_| | | | (_| |  __/ |_
   *      |_|\__,_|_|  \__, |\___|\__|
   *                   |___/
   */

  target_struct*
  yoga_ao_getyTarget(int argc, int pos) {
    return (target_struct *) yget_obj(argc - pos, &yTarget);
  }

  void target_free(void *obj) {
    target_struct *handler = (target_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_target *target_obj_handler =
          (sutra_target *) (handler->sutra_target);
      delete target_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void target_print(void *obj) {
    target_struct *handler = (target_struct *) obj;
    sutra_target *target_obj_handler = (sutra_target *) (handler->sutra_target);
    cout << "Yoga Target Object : " << endl;
    cout << "Contains " << target_obj_handler->ntargets << " target(s) : "
         << endl;
    cout << "Source #" << " | " << "position(\")" << " | " << " Mag " << " | "
         << "Lambda (mic.)" << endl;
    for (int i = 0; i < target_obj_handler->ntargets; i++) {
      cout << setw(8) << i + 1 << " | " << setw(4)
           << target_obj_handler->d_targets.at(i)->tposx << " , " << setw(4)
           << left << target_obj_handler->d_targets.at(i)->tposy << " | "
           << right << setw(5) << target_obj_handler->d_targets.at(i)->mag
           << " | " << target_obj_handler->d_targets.at(i)->lambda << endl;
    }
  }

  void Y_yoga_target(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    try {
      carma_context *context_handle = _getCurrentContext();
      int odevice = context_handle->get_activeDevice();
      ;
      int ntargets = ygets_i(argc - 1);

      float *xpos;
      xpos = ygeta_f(argc - 2, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for xpos");

      float *ypos;
      ypos = ygeta_f(argc - 3, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens ypos");

      float *lambda;
      lambda = ygeta_f(argc - 4, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens lambda");

      float *mag;
      mag = ygeta_f(argc - 5, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens mag");

      float zerop;
      zerop = ygets_f(argc - 6);

      long *sizes;
      sizes = ygeta_l(argc - 7, &ntot, dims);

      float *pup;
      pup = ygeta_f(argc - 8, &ntot, dims);

      int Npts;
      Npts = ygets_i(argc - 9);

      if (argc > 9)
        odevice = ygets_i(argc - 10);

      target_struct *handle = (target_struct *) ypush_obj(
          &yTarget, sizeof(target_struct));
      handle->device = odevice;

      handle->sutra_target = new sutra_target(context_handle, ntargets, xpos,
                                              ypos, lambda, mag, zerop, sizes,
                                              pup, Npts, odevice);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_target construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_target_addlayer(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int ntarget = ygets_i(argc - 2);
    char *type = ygets_q(argc - 3);
    float alt = ygets_f(argc - 4);
    float xoff = ygets_f(argc - 5);
    float yoff = ygets_f(argc - 6);

    target_handler->d_targets.at(ntarget)->add_layer(type, alt, xoff, yoff);
  }

  void Y_target_init_strehlmeter(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    int ntarget = ygets_i(argc - 2);

    target_handler->d_targets.at(ntarget)->init_strehlmeter();
  }

  void Y_target_reset_strehlmeter(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    int ntarget = ygets_i(argc - 2);

    target_handler->d_targets.at(ntarget)->reset_strehlmeter();
  }

  void Y_target_atmostrace(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    atmos_struct *handler_a = (atmos_struct *) yget_obj(argc - 3, &yAtmos);
    sutra_atmos *atmos_handler = (sutra_atmos *) (handler_a->sutra_atmos);
    target_handler->d_targets.at(ntarget)->raytrace(atmos_handler);
  }

  void Y_target_getimage(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    char *type_im = ygets_q(argc - 3);

    int puponly = 0;
    if (argc > 3)
      puponly = (int) ygets_l(argc - 4);

    float *data = ypush_f(
        (long*) target_handler->d_targets[ntarget]->d_image->getDims());

    target_handler->d_targets[ntarget]->comp_image(puponly);
    carma_obj<float> *d_img;
    if (strcmp(type_im, "se") == 0) {
      d_img = target_handler->d_targets[ntarget]->d_image;
    } else if (strcmp(type_im, "le") == 0) {
      d_img = target_handler->d_targets[ntarget]->d_leimage;
    } else {
      y_error("type unknown");
      return;
    }
    carma_obj<float> tmp_img(context_handle, d_img->getDims());
    float flux = target_handler->d_targets[ntarget]->zp
        * powf(10, -0.4 * target_handler->d_targets[ntarget]->mag);
    roll_mult<float>(
        tmp_img.getData(),
        d_img->getData(),
        d_img->getDims(1),
        d_img->getDims(2),
        flux,
        target_handler->d_targets[ntarget]->current_context->get_device(
            target_handler->d_targets[ntarget]->device));
    tmp_img.device2host(data);
  }

  void Y_target_resetphase(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    carmaSafeCall(
        cudaMemset(
            target_handler->d_targets.at(ntarget)->d_phase->d_screen->getData(),
            0,
            sizeof(float) * target_handler->d_targets.at(ntarget)->d_phase->d_screen->getNbElem()));
  }

  void Y_target_getphase(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    float *data =
        ypush_f(
            (long*) target_handler->d_targets.at(ntarget)->d_phase->d_screen->getDims());
    target_handler->d_targets.at(ntarget)->d_phase->d_screen->device2host(data);
  }

  void Y_target_getphasetele(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    float *data =
        ypush_f(
            (long*) target_handler->d_targets.at(ntarget)->phase_telemetry->getDims());
    target_handler->d_targets.at(ntarget)->phase_telemetry->fill_into(data);
  }

  void Y_target_getamplipup(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    long * ndims_data = new long[4];
    ndims_data[0] = 3;
    ndims_data[1] = 2;
    const long * ndims_obj =
        target_handler->d_targets.at(ntarget)->d_amplipup->getDims();
    memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 2);
    float *data = ypush_f(ndims_data);
    target_handler->d_targets.at(ntarget)->d_amplipup->device2host(
        (cuFloatComplex *) data);
    delete ndims_data;
  }

  void Y_target_getstrehl(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    long * ndims_data = new long[2];
    ndims_data[0] = 1;
    ndims_data[1] = 4;
    float * strehl = ypush_f(ndims_data);
    target_handler->d_targets.at(ntarget)->comp_image(0);
    target_handler->d_targets.at(ntarget)->comp_strehl();
    strehl[0] = target_handler->d_targets.at(ntarget)->strehl_se;
    strehl[1] = target_handler->d_targets.at(ntarget)->strehl_le;
    strehl[2] = target_handler->d_targets.at(ntarget)->phase_var;
    strehl[3] = target_handler->d_targets.at(ntarget)->phase_var_avg
        / (target_handler->d_targets.at(ntarget)->phase_var_count + 0.00001);
    delete ndims_data;
  }

  /*
   *     ____  _
   *    |  _ \| |__   __ _ ___  ___
   *    | |_) | '_ \ / _` / __|/ _ \
   *    |  __/| | | | (_| \__ \  __/
   *    |_|   |_| |_|\__,_|___/\___|
   *
   */

  phase_struct*
  yoga_ao_getyPhase(int argc, int pos) {
    return (phase_struct *) yget_obj(argc - pos, &yPhase);
  }

  void phase_free(void *obj) {
    phase_struct *handler = (phase_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_phase *phase_obj_handler = (sutra_phase *) (handler->sutra_phase);
      delete phase_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void phase_print(void *obj) {
    phase_struct *handler = (phase_struct *) obj;
    sutra_phase *phase_obj_handler = (sutra_phase *) (handler->sutra_phase);
    cout << "Yoga phase Object : " << phase_obj_handler->screen_size << "x"
         << phase_obj_handler->screen_size << endl;
  }

  void phase_eval(void *obj, int n) {
    phase_struct *handler = (phase_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);
    sutra_phase *phase_handler = (sutra_phase *) (handler->sutra_phase);
    float *data = ypush_f((long*) phase_handler->d_screen->getDims());
    phase_handler->d_screen->device2host(data);
  }

  void Y_yoga_phase(int argc) {
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();
      int odevice = activeDevice;
      long size = ygets_l(argc - 1);

      if (argc > 1)
        odevice = ygets_i(argc - 2);
      activeDevice = context_handle->set_activeDevice(odevice,1);

      phase_struct *handle = (phase_struct *) ypush_obj(&yPhase,
                                                        sizeof(phase_struct));
      handle->device = odevice;

      handle->sutra_phase = new sutra_phase(context_handle, size);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_phase construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_phase_set(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    phase_struct *handle = (phase_struct *) yget_obj(argc - 1, &yPhase);
    sutra_phase *phase_handler = (sutra_phase *) handle->sutra_phase;
    float *a = ygeta_f(argc - 2, &ntot, dims);
    if (ntot != phase_handler->d_screen->getNbElem())
      y_error("wrong size for array");

    phase_handler->d_screen->host2device(a);

  }

  /*
   *    __        _______ ____
   *    \ \      / /  ___/ ___|
   *     \ \ /\ / /| |_  \___ \
   *      \ V  V / |  _|  ___) |
   *       \_/\_/  |_|   |____/
   *
   */
  /*
   wfs_struct*
   yoga_ao_getyWfs(int argc, int pos) {
   return (wfs_struct *) yget_obj(argc - pos, &yWfs);
   }

   void wfs_free(void *obj) {
   wfs_struct *handler = (wfs_struct *) obj;
   carma_context *context_handle = _getCurrentContext();
   context_handle->set_activeDevice(handler->device,1);
   try {
   sutra_wfs *wfs_obj_handler = (sutra_wfs *) (handler->sutra_wfs);
   delete wfs_obj_handler;
   } catch (string &msg) {
   y_error(msg.c_str());
   } catch (char const * msg) {
   y_error(msg);
   }
   }

   void wfs_print(void *obj) {
   //wfs_struct *handler = (wfs_struct *)obj;
   //sutra_wfs *wfs_obj_handler = (sutra_wfs *)(handler->sutra_wfs);
   cout << "Yoga wfs Object : " << endl;
   }

   void Y_yoga_wfs(int argc)
   //long nxsub, long nvalid, long npix, long nrebin, long nfft
   {
   try {
   carma_context *context_handle = _getCurrentContext();
   int activeDevice = context_handle->get_activeDevice();
   int odevice = activeDevice;
   long nxsub = ygets_l(argc - 1);
   long nvalid = ygets_l(argc - 2);
   long npix = ygets_l(argc - 3);
   long nphase = ygets_l(argc - 4);
   long nrebin = ygets_l(argc - 5);
   long nfft = ygets_l(argc - 6);
   long ntot = ygets_l(argc - 7);
   long npup = ygets_l(argc - 8);
   float pdiam = ygets_f(argc - 9);
   float nphot = ygets_f(argc - 10);
   int lgs = ygets_l(argc - 11);

   if (argc > 11)
   odevice = ygets_i(argc - 12);
   activeDevice = context_handle->set_activeDevice(odevice,1);

   wfs_struct *handle = (wfs_struct *) ypush_obj(&yWfs, sizeof(wfs_struct));
   handle->device = odevice;

   handle->sutra_wfs = new sutra_wfs(context_handle, "sh", nxsub, nvalid, npix,
   nphase, nrebin, nfft, ntot, npup, pdiam, nphot, lgs, odevice);

   } catch (string &msg) {
   y_error(msg.c_str());
   } catch (char const * msg) {
   y_error(msg);
   } catch (...) {
   stringstream buf;
   buf << "unknown error with sutra_wfs construction in " << __FILE__ << "@"
   << __LINE__ << endl;
   y_error(buf.str().c_str());
   }
   }

   void Y_wfs_initgs(int argc) {
   wfs_struct *handle = (wfs_struct *) yget_obj(argc - 1, &yWfs);
   float xpos = ygets_f(argc - 2);
   float ypos = ygets_f(argc - 3);
   float lambda = ygets_f(argc - 4);
   float mag = ygets_f(argc - 5);
   long size = ygets_l(argc - 6);
   float noise = -1;
   long seed = 1234;
   if (argc > 6)
   noise = ygets_f(argc - 7);
   if (argc > 7)
   seed = ygets_l(argc - 8);

   carma_context *context_handle = _getCurrentContext();
   context_handle->set_activeDevice(handle->device,1);

   sutra_wfs *wfs_handler = (sutra_wfs *) handle->sutra_wfs;
   wfs_handler->wfs_initgs(xpos, ypos, lambda, mag, size, noise, seed);
   }
   */
  /*
   *     _       _                     _
   *    | |_ ___| | ___ _ __ ___   ___| |_ _ __ _   _
   *    | __/ _ \ |/ _ \ '_ ` _ \ / _ \ __| '__| | | |
   *    | ||  __/ |  __/ | | | | |  __/ |_| |  | |_| |
   *     \__\___|_|\___|_| |_| |_|\___|\__|_|   \__, |
   *                                            |___/
   */

  telemetry_struct*
  yoga_ao_getyTelemetry(int argc, int pos) {
    return (telemetry_struct *) yget_obj(argc - pos, &yTelemetry);
  }

  void telemetry_free(void *obj) {
    telemetry_struct *handler = (telemetry_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_telemetry *telemetry_obj_handler =
          (sutra_telemetry *) (handler->sutra_telemetry);
      delete telemetry_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void telemetry_print(void *obj) {
    telemetry_struct *handler = (telemetry_struct *) obj;
    sutra_telemetry *telemetry_obj_handler =
        (sutra_telemetry *) (handler->sutra_telemetry);
    cout << "Yoga Telemetry Object : " << endl;
    cout << "nb objs : " << telemetry_obj_handler->get_nbObjs() << endl;
    cout << "nb streams : " << telemetry_obj_handler->get_nbStreams() << endl;
  }

  void Y_yoga_telemetry(int argc) {
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();
      int odevice = activeDevice;

      //      long *size_obj   = ygeta_l(argc-1, &ntot,dims);
      //      carma_host_obj<float> *yho_tmp = new carma_host_obj<float>(size_obj);
      //      int nbStreams = ygets_i(argc-2);
      //
      telemetry_struct *handle = (telemetry_struct *) ypush_obj(
          &yTelemetry, sizeof(telemetry_struct));
      handle->device = odevice;
      handle->sutra_telemetry = new sutra_telemetry();

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_telemetry construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  /*
   *     ____
   *    / ___|  ___ _ __  ___  ___  _ __ ___
   *    \___ \ / _ \ '_ \/ __|/ _ \| '__/ __|
   *     ___) |  __/ | | \__ \ (_) | |  \__ \
   *    |____/ \___|_| |_|___/\___/|_|  |___/
   *
   */

  sensors_struct*
  yoga_ao_getySensors(int argc, int pos) {
    return (sensors_struct *) yget_obj(argc - pos, &ySensors);
  }

  void sensors_free(void *obj) {
    sensors_struct *handler = (sensors_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_sensors *sensors_obj_handler =
          (sutra_sensors *) (handler->sutra_sensors);
      delete sensors_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void sensors_print(void *obj) {
    sensors_struct *handler = (sensors_struct *) obj;
    sutra_sensors *sensors_obj_handler =
        (sutra_sensors *) (handler->sutra_sensors);
    cout << "Yoga sensors Object" << endl;

    cout << "Contains " << sensors_obj_handler->nsensors() << " WFS(s) : "
         << endl;
    cout << "WFS #" << " | " << "Nsubaps" << " | " << "Nvalid" << " | "
         << "Npix" << " | " << "Nphase" << " | " << "Nfft" << " | " << "Nrebin"
         << " | " << "Ntot" << " | " << "Npup" << endl;

    for (size_t i = 0; i < sensors_obj_handler->nsensors(); i++) {
      cout << setw(5) << i + 1 << " | " << setw(3)
           << sensors_obj_handler->d_wfs.at(i)->nxsub << "x" << left << setw(3)
           << sensors_obj_handler->d_wfs.at(i)->nxsub << " | " << right
           << setw(6) << sensors_obj_handler->d_wfs.at(i)->nvalid << " | "
           << setw(4) << sensors_obj_handler->d_wfs.at(i)->npix << " | "
           << setw(6) << sensors_obj_handler->d_wfs.at(i)->nphase << " | "
           << setw(4) << sensors_obj_handler->d_wfs.at(i)->nfft << " | "
           << setw(6) << sensors_obj_handler->d_wfs.at(i)->nrebin << " | "
           << setw(4) << sensors_obj_handler->d_wfs.at(i)->ntot << " | "
           << setw(4) << sensors_obj_handler->d_wfs.at(i)->npup << endl;
    }
  }

  void Y_yoga_sensors(int argc)
//long *nxsub,long *nvalid,long *npix,long *nrebin,long *nfft
      {
    long ntot;
    long dims[Y_DIMSIZE];

    try {
      int nsensors = ygets_i(argc - 1);
      ystring_t *type_data = ygeta_q(argc - 2, &ntot, dims);
      if (strcmp(type_data[0], "geo") == 0) {
        long *nxsub = ygeta_l(argc - 3, &ntot, dims);
        long *nvalid = ygeta_l(argc - 4, &ntot, dims);
        long *nphase = ygeta_l(argc - 5, &ntot, dims);
        long npup = ygets_l(argc - 6);
        float *pdiam = ygeta_f(argc - 7, &ntot, dims);

        carma_context *context_handle = _getCurrentContext();
        int odevice = context_handle->get_activeDevice();
        if (argc > 7)
          odevice = ygets_i(argc - 8);

        sensors_struct *handle = (sensors_struct *) ypush_obj(
            &ySensors, sizeof(sensors_struct));
        handle->device = odevice;

        //done inside sutra_wfs constructor
        //odevice = context_handle->set_activeDevice(odevice);

        handle->sutra_sensors = new sutra_sensors(context_handle, nsensors,
                                                  nxsub, nvalid, nphase, npup,
                                                  pdiam, odevice);

      } else {
        long *nxsub = ygeta_l(argc - 3, &ntot, dims);
        long *nvalid = ygeta_l(argc - 4, &ntot, dims);
        long *npix = ygeta_l(argc - 5, &ntot, dims);
        long *nphase = ygeta_l(argc - 6, &ntot, dims);
        long *nrebin = ygeta_l(argc - 7, &ntot, dims);
        long *nfft = ygeta_l(argc - 8, &ntot, dims);
        long *ntota = ygeta_l(argc - 9, &ntot, dims);
        long *npup = ygeta_l(argc - 10, &ntot, dims);
        float *pdiam = ygeta_f(argc - 11, &ntot, dims);
        float *nphot = ygeta_f(argc - 12, &ntot, dims);
        int *lgs = ygeta_i(argc - 13, &ntot, dims);

        carma_context *context_handle = _getCurrentContext();
        int odevice = context_handle->get_activeDevice();
        if (argc > 13)
          odevice = ygets_i(argc - 14);

        sensors_struct *handle = (sensors_struct *) ypush_obj(
            &ySensors, sizeof(sensors_struct));
        handle->device = odevice;

        //done inside sutra_wfs constructor
        //odevice = context_handle->set_activeDevice(odevice);
        sutra_sensors *sensors = new sutra_sensors(context_handle, type_data,
                                                   nsensors, nxsub, nvalid,
                                                   npix, nphase, nrebin, nfft,
                                                   ntota, npup, pdiam, nphot,
                                                   lgs, odevice);
        sensors->allocate_buffers();
        handle->sutra_sensors = (void*) sensors;
      }
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_sensors construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_sensors_initgs(int argc) {
    long int ntot;
    long dims[Y_DIMSIZE];

    sensors_struct *handle = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) handle->sutra_sensors;
    float *xpos = ygeta_f(argc - 2, &ntot, dims);
    if (ntot != abs(sensors_handler->nsensors()))
      y_error("wrong dimension for xpos");
    float *ypos = ygeta_f(argc - 3, &ntot, dims);
    if (ntot != abs(sensors_handler->nsensors()))
      y_error("wrong dimension for ypos");
    float *lambda = ygeta_f(argc - 4, &ntot, dims);
    if (ntot != abs(sensors_handler->nsensors()))
      y_error("wrong dimension for lambda");
    float *mag = ygeta_f(argc - 5, &ntot, dims);
    if (ntot != abs(sensors_handler->nsensors()))
      y_error("wrong dimension for mag");
    float zerop = ygets_f(argc - 6);
    long *size = ygeta_l(argc - 7, &ntot, dims);
    if (ntot != abs(sensors_handler->nsensors()))
      y_error("wrong dimension for size");
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handle->device,1);

    if (argc > 8) {
      float *noise = ygeta_f(argc - 8, &ntot, dims);
      long *seed = ygeta_l(argc - 9, &ntot, dims);
      sensors_handler->sensors_initgs(xpos, ypos, lambda, mag, zerop, size,
                                      noise, seed);
    } else if (argc > 7) {
      float *noise = ygeta_f(argc - 8, &ntot, dims);
      sensors_handler->sensors_initgs(xpos, ypos, lambda, mag, zerop, size,
                                      noise);
    } else {
      sensors_handler->sensors_initgs(xpos, ypos, lambda, mag, zerop, size);

    }
  }

  void Y_sensors_addlayer(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);
    char *type = ygets_q(argc - 3);
    float alt = ygets_f(argc - 4);
    float xoff = ygets_f(argc - 5);
    float yoff = ygets_f(argc - 6);

    sensors_handler->d_wfs.at(nsensor)->d_gs->add_layer(type, alt, xoff, yoff);
  }

  void Y_sensors_initarr(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);
    string type = sensors_handler->d_wfs.at(nsensor)->type;
    if (type == "sh") {
      int *phasemap = ygeta_i(argc - 3, &ntot, dims);
      int *hrmap = ygeta_i(argc - 4, &ntot, dims);
      int *binmap = ygeta_i(argc - 5, &ntot, dims);
      float *offsets = ygeta_f(argc - 6, &ntot, dims);
      float *pupil = ygeta_f(argc - 7, &ntot, dims);
      float *fluxPerSub = ygeta_f(argc - 8, &ntot, dims);
      int *validsubsx = ygeta_i(argc - 9, &ntot, dims);
      int *validsubsy = ygeta_i(argc - 10, &ntot, dims);
      int *istart = ygeta_i(argc - 11, &ntot, dims);
      int *jstart = ygeta_i(argc - 12, &ntot, dims);
      float *kernel = ygeta_f(argc - 13, &ntot, dims);

      sutra_wfs_sh *wfs =
          dynamic_cast<sutra_wfs_sh *>(sensors_handler->d_wfs.at(nsensor));
      wfs->wfs_initarrays(phasemap, hrmap, binmap, offsets, pupil, fluxPerSub,
                          validsubsx, validsubsy, istart, jstart,
                          (cuFloatComplex*) kernel);
    }
    if (type == "roof") {
      float *halfxy = ygeta_f(argc - 3, &ntot, dims);
      float *offsets = ygeta_f(argc - 4, &ntot, dims);
      float *focmask = ygeta_f(argc - 5, &ntot, dims);
      float *pupil = ygeta_f(argc - 6, &ntot, dims);
      int *cx = ygeta_i(argc - 7, &ntot, dims);
      int *cy = ygeta_i(argc - 8, &ntot, dims);
      float *sincar = ygeta_f(argc - 9, &ntot, dims);
      int *phasemap = ygeta_i(argc - 10, &ntot, dims);
      int *validx = ygeta_i(argc - 11, &ntot, dims);
      int *validy = ygeta_i(argc - 12, &ntot, dims);

      sutra_wfs_pyr_roof *wfs =
          dynamic_cast<sutra_wfs_pyr_roof *>(sensors_handler->d_wfs.at(nsensor));
      wfs->wfs_initarrays((cuFloatComplex*) halfxy, (cuFloatComplex*) offsets,
                          focmask, pupil, cx, cy, sincar, phasemap, validx,
                          validy);
    }
    if (type == "pyr") {
      float *halfxy = ygeta_f(argc - 3, &ntot, dims);
      float *offsets = ygeta_f(argc - 4, &ntot, dims);
      float *focmask = ygeta_f(argc - 5, &ntot, dims);
      float *pupil = ygeta_f(argc - 6, &ntot, dims);
      int *cx = ygeta_i(argc - 7, &ntot, dims);
      int *cy = ygeta_i(argc - 8, &ntot, dims);
      float *sincar = ygeta_f(argc - 9, &ntot, dims);
      int *phasemap = ygeta_i(argc - 10, &ntot, dims);
      int *validx = ygeta_i(argc - 11, &ntot, dims);
      int *validy = ygeta_i(argc - 12, &ntot, dims);

      sutra_wfs_pyr_pyr4 *wfs =
          dynamic_cast<sutra_wfs_pyr_pyr4 *>(sensors_handler->d_wfs.at(nsensor));
      wfs->wfs_initarrays((cuFloatComplex*) halfxy, (cuFloatComplex*) offsets,
                          focmask, pupil, cx, cy, sincar, phasemap, validx,
                          validy);
    }
    if (sensors_handler->d_wfs.at(nsensor)->type == "geo") {
      int *phasemap = ygeta_i(argc - 3, &ntot, dims);
      float *offsets = ygeta_f(argc - 4, &ntot, dims);
      float *pupil = ygeta_f(argc - 5, &ntot, dims);
      float *fluxPerSub = ygeta_f(argc - 6, &ntot, dims);
      int *validsubsx = ygeta_i(argc - 7, &ntot, dims);
      int *validsubsy = ygeta_i(argc - 8, &ntot, dims);

      sutra_wfs_geom *wfs =
          dynamic_cast<sutra_wfs_geom *>(sensors_handler->d_wfs.at(nsensor));
      wfs->wfs_initarrays(phasemap, offsets, pupil, fluxPerSub, validsubsx,
                          validsubsy);
    }
  }

  void Y_sensors_loadkernels(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    float *kernels = ygeta_f(argc - 3, &ntot, dims);

    sensors_handler->d_wfs.at(nsensor)->load_kernels(kernels);
  }

  void Y_sensors_initlgs(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    int nprof = ygets_i(argc - 3);
    float hg = ygets_f(argc - 4);
    float h0 = ygets_f(argc - 5);
    float deltah = ygets_f(argc - 6);
    float pixsize = ygets_f(argc - 7);
    float *doffaxis = ygeta_f(argc - 8, &ntot, dims);
    float *prof1d = ygeta_f(argc - 9, &ntot, dims);
    float *profcum = ygeta_f(argc - 10, &ntot, dims);
    float *beam = ygeta_f(argc - 11, &ntot, dims);
    float *ftbeam = ygeta_f(argc - 12, &ntot, dims);
    float *azimuth = ygeta_f(argc - 13, &ntot, dims);

    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_init(
        nprof, hg, h0, deltah, pixsize, doffaxis, prof1d, profcum, beam,
        (cuFloatComplex *) ftbeam, azimuth);
    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_update(
        context_handle->get_device(handler->device));
    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_makespot(
        context_handle->get_device(handler->device), 0);
  }

  void Y_sensors_updatelgsprof(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    float *prof1d = ygeta_f(argc - 3, &ntot, dims);
    float *profcum = ygeta_f(argc - 4, &ntot, dims);
    float hg = ygets_f(argc - 5);
    float h0 = ygets_f(argc - 6);
    float deltah = ygets_f(argc - 7);

    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->load_prof(prof1d, profcum,
                                                               hg, h0, deltah);
    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_update(
        context_handle->get_device(handler->device));
    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_makespot(
        context_handle->get_device(handler->device), 0);
  }

  void Y_sensors_updatelgs(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_update(
        context_handle->get_device(handler->device));
    sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->lgs_makespot(
        context_handle->get_device(handler->device), 0);
  }

  void Y_sensors_fillbinimage(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    int nsensor = ygets_i(argc - 2);

    /*
     if(argc>2) { // With a telemetry object
     carma_host_struct *handle_obj = (carma_host_struct *)yget_obj(piargs[0],&carma_host_caObj);telemetry_struct *handler2 = (telemetry_struct *)yget_obj(argc-3,&yTelemetry);
     sutra_telemetry *telemetry_handler = (sutra_telemetry *)(handler2->sutra_telemetry);

     string type_obj = string(ygets_q(argc-4));

     sensors_handler->d_wfs.at(nsensor)->comp_image(telemetry_handler, type_obj, nsensor);
     } else {
     sensors_handler->d_wfs.at(nsensor)->comp_image();
     }
     */
    if (sensors_handler->d_wfs.at(nsensor)->type != "sh") {
      y_error("wfs should be a SH");
    }
    sutra_wfs_sh *wfs = dynamic_cast<sutra_wfs_sh *>(sensors_handler->d_wfs.at(
        nsensor));
    wfs->fill_binimage(0);
    //sensors_handler->d_wfs.at(nsensor)->streams->wait_all_streams();
  }

  void Y_sensors_compimg(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    int nsensor = ygets_i(argc - 2);

    /*
     if(argc>2) { // With a telemetry object
     carma_host_struct *handle_obj = (carma_host_struct *)yget_obj(piargs[0],&carma_host_caObj);telemetry_struct *handler2 = (telemetry_struct *)yget_obj(argc-3,&yTelemetry);
     sutra_telemetry *telemetry_handler = (sutra_telemetry *)(handler2->sutra_telemetry);

     string type_obj = string(ygets_q(argc-4));

     sensors_handler->d_wfs.at(nsensor)->comp_image(telemetry_handler, type_obj, nsensor);
     } else {
     sensors_handler->d_wfs.at(nsensor)->comp_image();
     }
     */
    sensors_handler->d_wfs.at(nsensor)->comp_image();
    //sensors_handler->d_wfs.at(nsensor)->streams->wait_all_streams();
  }

  void Y_sensors_getimg(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    float *data = ypush_f(
        (long*) sensors_handler->d_wfs.at(nsensor)->d_binimg->getDims());
    /*
     if(argc>2) { // With a telemetry object
     telemetry_struct *handler2 = (telemetry_struct *)yget_obj(argc-3,&yTelemetry);
     sutra_telemetry *telemetry_handler = (sutra_telemetry *)(handler2->sutra_telemetry);

     string type_obj = string(ygets_q(argc-4));

     telemetry_handler->fill_into(type_obj, nsensor, data);
     telemetry_handler->wait_obj(type_obj, nsensor);
     } else
     */
    if (sensors_handler->d_wfs.at(nsensor)->type == "sh")
      dynamic_cast<sutra_wfs_sh*>(sensors_handler->d_wfs.at(nsensor))->fill_binimage(
          0);
    sensors_handler->d_wfs.at(nsensor)->d_binimg->device2host(data);
  }

  void Y_sensors_getdata(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    char *type_data = ygets_q(argc - 3);
    if (strcmp(type_data, "amplipup") == 0) {
      if (sensors_handler->d_wfs.at(nsensor)->type == "sh") {
        long *ndims_data = new long[5];
        ndims_data[0] = 4;
        ndims_data[1] = 2;
        const long *ndims_obj =
            sensors_handler->d_wfs.at(nsensor)->d_camplipup->getDims();
        memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
        float *data = ypush_f(ndims_data);
        sensors_handler->d_wfs.at(nsensor)->d_camplipup->device2host(
            (cuFloatComplex *) data);
      }
      if (sensors_handler->d_wfs.at(nsensor)->type == "pyr") {
        long *ndims_data = new long[4];
        ndims_data[0] = 3;
        ndims_data[1] = 2;
        const long *ndims_obj =
            sensors_handler->d_wfs.at(nsensor)->d_camplipup->getDims();
        memcpy(&(ndims_data[2]), &(ndims_obj[1]), sizeof(long) * 2);
        float *data = ypush_f(ndims_data);
        sensors_handler->d_wfs.at(nsensor)->d_camplipup->device2host(
            (cuFloatComplex *) data);
      }
    }
    if (strcmp(type_data, "amplifoc") == 0) {
      if (sensors_handler->d_wfs.at(nsensor)->type == "sh") {
        long *ndims_data = new long[5];
        ndims_data[0] = 4;
        ndims_data[1] = 2;
        const long *ndims_obj =
            sensors_handler->d_wfs.at(nsensor)->d_camplifoc->getDims();
        memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
        float *data = ypush_f(ndims_data);
        sensors_handler->d_wfs.at(nsensor)->d_camplifoc->device2host(
            (cuFloatComplex *) data);
      }
      if ((sensors_handler->d_wfs.at(nsensor)->type == "pyr") || (sensors_handler->d_wfs.at(
          nsensor)->type
                                                                  == "roof")) {
        long *ndims_data = new long[4];
        ndims_data[0] = 3;
        ndims_data[1] = 2;
        const long *ndims_obj =
            sensors_handler->d_wfs.at(nsensor)->d_camplifoc->getDims();
        memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 2);
        float *data = ypush_f(ndims_data);
        sensors_handler->d_wfs.at(nsensor)->d_camplifoc->device2host(
            (cuFloatComplex *) data);
      }
    }
    if (strcmp(type_data, "hrimg") == 0) {
      if ((sensors_handler->d_wfs.at(nsensor)->type == "pyr") || (sensors_handler->d_wfs.at(
          nsensor)->type
                                                                  == "roof")) {
        sutra_wfs_pyr *wfs =
            dynamic_cast<sutra_wfs_pyr *>(sensors_handler->d_wfs.at(nsensor));

        float *data = ypush_f((long*) wfs->d_hrimg->getDims());
        wfs->d_hrimg->device2host(data);
      }
    }
    if (strcmp(type_data, "fttotim") == 0) {
      long *ndims_data = new long[5];
      ndims_data[0] = 4;
      ndims_data[1] = 2;
      const long *ndims_obj =
          sensors_handler->d_wfs.at(nsensor)->d_fttotim->getDims();
      memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_fttotim->device2host(
          (cuFloatComplex *) data);
    }
    if (strcmp(type_data, "bincube") == 0) {
      float *data = ypush_f(
          (long*) sensors_handler->d_wfs.at(nsensor)->d_bincube->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_bincube->device2host(data);
    }
    if (strcmp(type_data, "phase") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->device2host(
          data);
    }
    if (strcmp(type_data, "pupil") == 0) {
      float *data = ypush_f(
          (long*) sensors_handler->d_wfs.at(nsensor)->d_pupil->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_pupil->device2host(data);
    }
    if (strcmp(type_data, "imgtele") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->image_telemetry->getDims());
      if (sensors_handler->d_wfs.at(nsensor)->type == "sh")
        dynamic_cast<sutra_wfs_sh*>(sensors_handler->d_wfs.at(nsensor))->fill_binimage(
            1);

      sensors_handler->d_wfs.at(nsensor)->image_telemetry->fill_into(data);
    }
    if (strcmp(type_data, "phasetele") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->d_gs->phase_telemetry->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_gs->phase_telemetry->fill_into(
          data);
    }
    if (strcmp(type_data, "lgskern") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_lgskern->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_lgskern->device2host(
          data);
    }
    if (strcmp(type_data, "ftlgskern") == 0) {
      long *ndims_data = new long[5];
      ndims_data[0] = 4;
      ndims_data[1] = 2;
      const long *ndims_obj =
          sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_ftlgskern->getDims();
      memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_ftlgskern->device2host(
          (cuFloatComplex *) data);
    }
    if (strcmp(type_data, "subsum") == 0) {
      float *data = ypush_f(
          (long*) sensors_handler->d_wfs.at(nsensor)->d_subsum->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_subsum->device2host(data);
    }
    if (strcmp(type_data, "psum") == 0) {
      if ((sensors_handler->d_wfs.at(nsensor)->type == "pyr") || (sensors_handler->d_wfs.at(
          nsensor)->type
                                                                  == "roof")) {
        sutra_wfs_pyr *wfs =
            dynamic_cast<sutra_wfs_pyr *>(sensors_handler->d_wfs.at(nsensor));

        float *data = ypush_f((long*) wfs->d_psum->getDims());
        wfs->d_psum->device2host(data);
      }
    }
    if (strcmp(type_data, "slopes") == 0) {
      float *data = ypush_f(
          (long*) sensors_handler->d_wfs.at(nsensor)->d_slopes->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_slopes->device2host(data);
    }
    if (strcmp(type_data, "offsets") == 0) {
      float *data = ypush_f(
          (long*) sensors_handler->d_wfs.at(nsensor)->d_offsets->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_offsets->device2host(data);
    }
    if (strcmp(type_data, "prof1d") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_prof1d->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_prof1d->device2host(
          data);
    }
    if (strcmp(type_data, "profcum") == 0) {
      float *data =
          ypush_f(
              (long*) sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_profcum->getDims());
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_profcum->device2host(
          data);
    }
    if (strcmp(type_data, "prof2d") == 0) {
      long *ndims_data = new long[4];
      ndims_data[0] = 3;
      ndims_data[1] = 2;
      const long *ndims_obj =
          sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_prof2d->getDims();
      memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 2);
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_prof2d->device2host(
          (cuFloatComplex *) data);
    }
    if (strcmp(type_data, "ftbeam") == 0) {
      long *ndims_data = new long[3];
      ndims_data[0] = 2;
      ndims_data[1] = 2;
      const long *ndims_obj =
          sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_ftbeam->getDims();
      ndims_data[2] = ndims_obj[1];
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_lgs->d_ftbeam->device2host(
          (cuFloatComplex *) data);
    }
    if (strcmp(type_data, "poffsets") == 0) {
      if ((sensors_handler->d_wfs.at(nsensor)->type == "pyr") || (sensors_handler->d_wfs.at(
          nsensor)->type
                                                                  == "roof")) {
        sutra_wfs_pyr *wfs =
            dynamic_cast<sutra_wfs_pyr *>(sensors_handler->d_wfs.at(nsensor));
        long *ndims_data = new long[4];
        ndims_data[0] = 3;
        ndims_data[1] = 2;
        const long *ndims_obj = wfs->d_poffsets->getDims();
        memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 2);
        float *data = ypush_f(ndims_data);
        wfs->d_poffsets->device2host((cuFloatComplex *) data);
      }
    }
    if (strcmp(type_data, "phalfxy") == 0) {
      if ((sensors_handler->d_wfs.at(nsensor)->type == "pyr") || (sensors_handler->d_wfs.at(
          nsensor)->type
                                                                  == "roof")) {
        sutra_wfs_pyr *wfs =
            dynamic_cast<sutra_wfs_pyr *>(sensors_handler->d_wfs.at(nsensor));

        long *ndims_data = new long[4];
        ndims_data[0] = 3;
        ndims_data[1] = 2;
        const long *ndims_obj = wfs->d_phalfxy->getDims();
        memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 2);
        float *data = ypush_f(ndims_data);
        wfs->d_phalfxy->device2host((cuFloatComplex *) data);
      }
    }

  }

  void Y_sensors_setphase(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (ntot != sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->getNbElem()) {
      stringstream buf;
      buf << "wrong size for phase array" << endl;
      y_error(buf.str().c_str());
    }

    sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->host2device(
        data);
  }

  void Y_sensors_resetphase(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    carmaSafeCall(
        cudaMemset(
            sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->getData(),
            0,
            sizeof(float) * sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->getNbElem()));
  }

  /*
   *                        _     _
   *   __ _  ___ __ _ _   _(_)___(_)_ __ ___
   *  / _` |/ __/ _` | | | | / __| | '_ ` _ \
   * | (_| | (_| (_| | |_| | \__ \ | | | | | |
   *  \__,_|\___\__, |\__,_|_|___/_|_| |_| |_|
   *               |_|
   */

  acquisim_struct*
  yoga_ao_getyAcquisim(int argc, int pos) {
    return (acquisim_struct *) yget_obj(argc - pos, &yAcquisim);
  }

  void acquisim_free(void *obj) {
    acquisim_struct *handler = (acquisim_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_acquisim *acquisim_obj_handler =
          (sutra_acquisim *) (handler->sutra_acquisim);
      delete acquisim_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void acquisim_print(void *obj) {
    //wfs_struct *handler = (wfs_struct *)obj;
    //sutra_wfs *wfs_obj_handler = (sutra_wfs *)(handler->sutra_wfs);
    cout << "Yoga acquisim Object : " << endl;
  }

  void Y_yoga_acquisim(int argc)
//long nxsub, long nvalid, long npix, long nrebin, long nfft
      {
    try {
      sensors_struct *sensors_handle = (sensors_struct *) yget_obj(argc - 1,
                                                                   &ySensors);
      sutra_sensors *sensors_handler =
          (sutra_sensors *) sensors_handle->sutra_sensors;
      int wfs_num = ygets_i(argc - 2);

      acquisim_struct *handle = (acquisim_struct *) ypush_obj(
          &yAcquisim, sizeof(acquisim_struct));
      handle->device = sensors_handle->device;

      handle->sutra_acquisim = new sutra_acquisim(sensors_handler, wfs_num);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_acquisim construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_acquisim_fillbcube(int argc)
//long nxsub, long nvalid, long npix, long nrebin, long nfft
      {
    try {
      acquisim_struct *handle = (acquisim_struct *) yget_obj(argc - 1,
                                                             &yAcquisim);
      sutra_acquisim *acquisim_handler =
          (sutra_acquisim *) handle->sutra_acquisim;
      long ntot;
      long dims[Y_DIMSIZE];
      float *a = ygeta_f(argc - 2, &ntot, dims);
      acquisim_handler->comp_image(dims, a);
      //acquisim_handler->comp_image_tele(dims,a);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_wfs construction in " << __FILE__ << "@"
          << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_acquisim_fillbcube_2D(int argc)
//long nxsub, long nvalid, long npix, long nrebin, long nfft
      {
    try {
      acquisim_struct *handle = (acquisim_struct *) yget_obj(argc - 1,
                                                             &yAcquisim);
      sutra_acquisim *acquisim_handler =
          (sutra_acquisim *) handle->sutra_acquisim;
      long ntot;
      long dims[Y_DIMSIZE];
      float *a = ygeta_f(argc - 2, &ntot, dims);
      long dims2[Y_DIMSIZE];
      int *num_ssp = ygeta_i(argc - 3, &ntot, dims2);
      acquisim_handler->comp_image_2D(dims, a, num_ssp);
      //acquisim_handler->comp_image_tele(dims,a);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_wfs construction in " << __FILE__ << "@"
          << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  /*
   *     ____  __  __
   *    |  _ \|  \/  |___
   *    | | | | |\/| / __|
   *    | |_| | |  | \__ \
 *    |____/|_|  |_|___/
   *
   */

  dms_struct*
  yoga_ao_getyDMs(int argc, int pos) {
    return (dms_struct *) yget_obj(argc - pos, &yDMs);
  }

  void dms_free(void *obj) {
    dms_struct *handler = (dms_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_dms *dms_obj_handler = (sutra_dms *) (handler->sutra_dms);
      delete dms_obj_handler;
    } catch (string& msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void dms_print(void *obj) {
    dms_struct *handler = (dms_struct *) obj;
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    cout << "Yoga DMs Object" << endl;

    cout << "Contains " << dms_handler->d_dms.size() << " DMs : " << endl;
    cout << "DM #" << " | " << "Type " << " | " << "  Alt  " << " | " << "Nact"
         << " | " << "Dim" << endl;
    vector<sutra_dm *>::iterator p;
    vector<type_screen>::iterator q;
    p = dms_handler->d_dms.begin();
    q = dms_handler->d_type.begin();
    int cpt = 0;
    while (p != dms_handler->d_dms.end()) {
      sutra_dm *dm = *p;
      cpt++;
      string types = (*q).first;
      float alt = (*q).second;
      long dim = dm->dim;
      long ninflu = dm->ninflu;
      cout << setw(4) << cpt << " | " << setw(5) << types << " | " << setw(7)
           << alt << " | " << setw(4) << ninflu << " | " << setw(4) << dim
           << endl;
      p++;
      q++;
    }
  }

  void Y_yoga_dms(int argc) {
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();

      long ndm = ygets_l(argc - 1);
      dms_struct *handle = (dms_struct *) ypush_obj(&yDMs, sizeof(dms_struct));
      handle->device = activeDevice;
      handle->sutra_dms = new sutra_dms(ndm);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_dms construction in " << __FILE__ << "@"
          << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_yoga_addpzt(int argc) {
    carma_context *context_handle = _getCurrentContext();
    int activeDevice = context_handle->get_activeDevice();
    int odevice = activeDevice;
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    long dim = ygets_l(argc - 3);
    long ninflu = ygets_l(argc - 4);
    long influsize = ygets_l(argc - 5);
    long ninflupos = ygets_l(argc - 6);
    long n_npts = ygets_l(argc - 7);
    float push4imat = ygets_f(argc - 8);

    if (argc > 8)
      odevice = ygets_i(argc - 9);

    activeDevice = context_handle->set_activeDevice(odevice,1);

    if (dms_handler->add_dm(context_handle, "pzt", alt, dim, ninflu, influsize,
                            ninflupos, n_npts, push4imat, odevice)
        == EXIT_FAILURE) {
      stringstream buf;
      buf << "unknown error with yoga_addpzt function in " << __FILE__ << "@"
          << __LINE__ << " DM(ptz," << alt << ") already exist" << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_yoga_rmdm(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    dms_handler->remove_dm(type, alt);
  }

  void Y_yoga_loadpzt(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    float *influ = ygeta_f(argc - 3, &ntot, dims);
    int *influpos = ygeta_i(argc - 4, &ntot, dims);
    int *npoints = ygeta_i(argc - 5, &ntot, dims);
    int *istart = ygeta_i(argc - 6, &ntot, dims);
    int *xoff = ygeta_i(argc - 7, &ntot, dims);
    int *yoff = ygeta_i(argc - 8, &ntot, dims);
    float *kern = ygeta_f(argc - 9, &ntot, dims);

    int inddm = dms_handler->get_inddm("pzt", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_loadpzt function in " << __FILE__ << "@"
          << __LINE__ << " DM(pzt," << alt << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->pzt_loadarrays(influ, influpos, npoints, istart,
                                              xoff, yoff, kern);
  }

  void Y_yoga_addkl(int argc) {
    carma_context *context_handle = _getCurrentContext();
    int activeDevice = context_handle->get_activeDevice();
    int odevice = activeDevice;
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    long dim = ygets_l(argc - 3);
    long ninflu = ygets_l(argc - 4);
    long influsize = ygets_l(argc - 5);
    long nr = ygets_l(argc - 6);
    long np = ygets_l(argc - 7);
    float push4imat = ygets_f(argc - 8);

    if (argc > 8)
      odevice = ygets_i(argc - 9);

    activeDevice = context_handle->set_activeDevice(odevice,1);

    if (dms_handler->add_dm(context_handle, "kl", alt, dim, ninflu, influsize,
                            nr, np, push4imat, odevice)
        == EXIT_FAILURE) {
      stringstream buf;
      buf << "unknown error with yoga_addkl function in " << __FILE__ << "@"
          << __LINE__ << " DM(kl," << alt << ") already exist" << endl;
      y_error(buf.str().c_str());
    }
  }
  void Y_yoga_loadkl(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    float *rabas = ygeta_f(argc - 3, &ntot, dims);
    float *azbas = ygeta_f(argc - 4, &ntot, dims);
    int *ord = ygeta_i(argc - 5, &ntot, dims);
    float *cr = ygeta_f(argc - 6, &ntot, dims);
    float *cp = ygeta_f(argc - 7, &ntot, dims);

    int inddm = dms_handler->get_inddm("kl", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_loadkl function in " << __FILE__ << "@"
          << __LINE__ << " DM(kl," << alt << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->kl_loadarrays(rabas, azbas, ord, cr, cp);
  }

  void Y_yoga_rmkl(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);

    int inddm = dms_handler->get_inddm("kl", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_loadkl function in " << __FILE__ << "@"
          << __LINE__ << " DM(kl," << alt << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->remove_dm("kl", alt);
  }

  void Y_yoga_getkl(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    long nkl = ygets_l(argc - 3);

    int inddm = dms_handler->get_inddm("kl", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_loadkl function in " << __FILE__ << "@"
          << __LINE__ << " DM(kl," << alt << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    int xoff =
        (int) ((dms_handler->d_dms[inddm]->d_shape->d_screen->getDims()[1]
            - dms_handler->d_dms[inddm]->d_kl->dim)
               / 2.0f);
    int yoff = xoff;

    dms_handler->d_dms[inddm]->d_kl->do_compute(
        dms_handler->d_dms[inddm]->d_shape->d_screen->getData(), nkl,
        dms_handler->d_dms[inddm]->d_shape->d_screen->getDims()[1], xoff, yoff);

    float *data = ypush_f(
        (long*) dms_handler->d_dms[inddm]->d_shape->d_screen->getDims());
    dms_handler->d_dms[inddm]->d_shape->d_screen->device2host(data);
  }

  void Y_yoga_addtt(int argc) {
    carma_context *context_handle = _getCurrentContext();
    int activeDevice = context_handle->get_activeDevice();

    int odevice = activeDevice;
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    float alt = ygets_f(argc - 2);
    long dim = ygets_l(argc - 3);
    float push4imat = ygets_f(argc - 4);
    if (argc > 4)
      odevice = ygets_i(argc - 5);

    activeDevice = context_handle->set_activeDevice(odevice,1);

    if (dms_handler->add_dm(context_handle, "tt", alt, dim, 2, dim, 1, 1,
                            push4imat, odevice)
        == EXIT_FAILURE) {
      stringstream buf;
      buf << "unknown error with yoga_addtt function in " << __FILE__ << "@"
          << __LINE__ << " DM(tt," << alt << ") already exist" << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_yoga_loadtt(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
    float alt = ygets_f(argc - 2);
    float *influ = ygeta_f(argc - 3, &ntot, dims);

    int inddm = dms_handler->get_inddm("tt", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_loadtt function in " << __FILE__ << "@"
          << __LINE__ << " DM(tt," << alt << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->d_influ->host2device(influ);
  }

  void Y_yoga_setcomm(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);
    float *comm = ygeta_f(argc - 4, &ntot, dims);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_setcomm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->d_comm->host2device(comm);
  }

  void Y_yoga_getcomm(int argc) {

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getcomm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    float *comm = ypush_f(
        const_cast<long*>(dms_handler->d_dms[inddm]->d_comm->getDims()));

    dms_handler->d_dms[inddm]->d_comm->device2host(comm);
  }

  void Y_yoga_getmapactu(int argc) {

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getmapactu function in " << __FILE__
          << "@" << __LINE__ << " DM(" << type << "," << alt
          << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    float *mapactu = ypush_f(
        const_cast<long*>(dms_handler->d_dms[inddm]->d_mapactu->getDims()));

    dms_handler->d_dms[inddm]->d_mapactu->device2host(mapactu);
  }

  void Y_yoga_getpos(int argc) {

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getpos function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    int *pos = ypush_i(
        const_cast<long*>(dms_handler->d_dms[inddm]->d_pos->getDims()));

    dms_handler->d_dms[inddm]->d_pos->device2host(pos);
  }

  void Y_yoga_shapedm(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_shapedm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->comp_shape();
  }

  void Y_yoga_resetdm(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_resetdm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->reset_shape();
  }

  void Y_yoga_oneactu(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);
    long nactu = ygets_l(argc - 4);
    float ampli = ygets_f(argc - 5);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_oneactu function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->comp_oneactu(nactu, ampli);
  }

  void Y_yoga_getdm(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getdm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    float *data = ypush_f(
        (long*) dms_handler->d_dms[inddm]->d_shape->d_screen->getDims());
    dms_handler->d_dms[inddm]->d_shape->d_screen->device2host(data);
  }

  void Y_yoga_setdm(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);
    float *data = ygeta_f(argc - 4, &ntot, dims);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_setdm function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->d_shape->d_screen->host2device(data);
  }

  void Y_yoga_computeKLbasis(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);
    float *xpos = ygeta_f(argc - 4, &ntot, dims);
    float *ypos = ygeta_f(argc - 5, &ntot, dims);
    int *indx_pup = ygeta_i(argc - 6, &ntot, dims);
    long dim = ygets_l(argc - 7);
    float norm = ygets_f(argc - 8);
    float ampli = ygets_f(argc - 9);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_computeKLbasis function in " << __FILE__
          << "@" << __LINE__ << " DM(" << type << "," << alt
          << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->compute_KLbasis(xpos, ypos, indx_pup, dim, norm,
                                               ampli);

  }

  void Y_yoga_setcomkl(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);
    float *comvec = ygeta_f(argc - 4, &ntot, dims);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_setcomkl function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    dms_handler->d_dms[inddm]->set_comkl(comvec);

  }

  void Y_yoga_getKLbasis(int argc) {
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getKLbasis function in " << __FILE__
          << "@" << __LINE__ << " DM(" << type << "," << alt
          << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }
    float *data = ypush_f(
        (long*) dms_handler->d_dms[inddm]->d_KLbasis->getDims());
    dms_handler->d_dms[inddm]->d_KLbasis->device2host(data);
  }

  void Y_target_dmtrace(int argc) {
    target_struct *handler = (target_struct *) yget_obj(argc - 1, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int ntarget = ygets_i(argc - 2);

    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);

    int rst = 0;
    if (argc > 3)
      rst = ygets_i(argc - 4);

    target_handler->d_targets.at(ntarget)->raytrace(dms_handler, rst);
  }

  void Y_dms_getdata(int argc) {
    dms_struct *dhandler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (dhandler->sutra_dms);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(dhandler->device,1);

    char *type = ygets_q(argc - 2);
    float alt = ygets_f(argc - 3);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getdata function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    char *type_data = ygets_q(argc - 4);
    if (strcmp(type_data, "xoff") == 0) {
      int *data = ypush_i((long*) dms_handler->d_dms[inddm]->d_xoff->getDims());
      dms_handler->d_dms[inddm]->d_xoff->device2host(data);
    }
    if (strcmp(type_data, "yoff") == 0) {
      int *data = ypush_i((long*) dms_handler->d_dms[inddm]->d_yoff->getDims());
      dms_handler->d_dms[inddm]->d_yoff->host2device(data);
    }
    if (strcmp(type_data, "rabas") == 0) {
      float *data = ypush_f(
          (long*) dms_handler->d_dms[inddm]->d_kl->d_rabas->getDims());
      dms_handler->d_dms[inddm]->d_kl->d_rabas->host2device(data);
    }
    if (strcmp(type_data, "azbas") == 0) {
      float *data = ypush_f(
          (long*) dms_handler->d_dms[inddm]->d_kl->d_azbas->getDims());
      dms_handler->d_dms[inddm]->d_kl->d_azbas->host2device(data);
    }
  }

  void Y_dms_comp_shape(int argc) {
    dms_struct *dhandler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (dhandler->sutra_dms);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(dhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    if (yarg_typeid(argc - 2) == Y_FLOAT) {
      float *com = ygeta_f(argc - 2, &ntot, dims);

      carma_obj<float> d_com(context_handle, const_cast<const long*>(dims),
                             com);

      vector<sutra_dm *>::iterator p;
      p = dms_handler->d_dms.begin();
      int idx = 0;
      while (p != dms_handler->d_dms.end()) {
        sutra_dm *dm = *p;
        dm->comp_shape(d_com[idx]);
        idx += dm->ninflu;
        p++;
      }
    }

  }

  /*
   *     ____ _____ ____
   *    |  _ \_   _/ ___|
   *    | |_) || || |
   *    |  _ < | || |___
   *    |_| \_\|_| \____|
   *
   */

  rtc_struct*
  yoga_ao_getyRTC(int argc, int pos) {
    return (rtc_struct *) yget_obj(argc - pos, &yRTC);
  }

  void rtc_free(void *obj) {
    rtc_struct *handler = (rtc_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_rtc *rtc_obj_handler = (sutra_rtc *) (handler->sutra_rtc);
      delete rtc_obj_handler;
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void rtc_print(void *obj) {
    rtc_struct *handler = (rtc_struct *) obj;
    sutra_rtc *rtc_handler = (sutra_rtc *) (handler->sutra_rtc);
    cout << "Yoga RTC Object" << endl;

    cout << "Contains " << rtc_handler->d_centro.size() << " Centroider(s) : "
         << endl;
    cout << "Centro #" << " | " << "Type " << " | " << "nwfs" << " | "
         << "Nvalid" << endl;

    for (size_t idx = 0; idx < rtc_handler->d_centro.size(); idx++) {
      cout << setw(8) << idx + 1 << " | " << setw(5)
           << rtc_handler->d_centro.at(idx)->get_type() << " | " << setw(4)
           << rtc_handler->d_centro.at(idx)->nwfs + 1 << " | "
           << rtc_handler->d_centro.at(idx)->nvalid << endl;
    }

    cout << "Contains " << rtc_handler->d_control.size() << " Controller(s) : "
         << endl;
    cout << "Control #" << " | " << "Type " << " | " << "Nslope" << " | "
         << "Nactu" << endl;

    for (size_t idx = 0; idx < rtc_handler->d_control.size(); idx++) {
      cout << setw(9) << idx + 1 << " | " << setw(5)
           << rtc_handler->d_control.at(idx)->get_type() << " | " << setw(6)
           << rtc_handler->d_control.at(idx)->nslope() << " | " << setw(5)
           << rtc_handler->d_control.at(idx)->nactu() << endl;
    }
  }

  void Y_yoga_rtc(int argc) {
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();
      if (argc == 1) {
        activeDevice = ygets_i(argc - 1);
        context_handle->set_activeDevice(activeDevice,1);
      }
      rtc_struct *handle = (rtc_struct *) ypush_obj(&yRTC, sizeof(rtc_struct));
      handle->use_brama = 0;
      handle->device = activeDevice;
      handle->sutra_rtc = new sutra_rtc(context_handle);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_rtc construction in " << __FILE__ << "@"
          << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_rtc_addcentro(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    sensors_struct *shandler = (sensors_struct *) yget_obj(argc - 2, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (shandler->sutra_sensors);
    ;
    long nwfs = ygets_l(argc - 3);
    long nvalid = ygets_l(argc - 4);
    char *type_centro = ygets_q(argc - 5);
    float offset = ygets_f(argc - 6);
    float scale = ygets_f(argc - 7);

    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    carma_context *context_handle = _getCurrentContext();
    int activeDevice =
        context_handle->set_activeDeviceForCpy(rhandler->device,1);

    rtc_handler->add_centroider(sensors_handler, nwfs, nvalid, offset, scale,
                                activeDevice, type_centro);
  }

  void Y_rtc_addcontrol(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long nactu = ygets_l(argc - 2);
    float delay = ygets_f(argc - 3);
    char *type_control = ygets_q(argc - 4);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 5, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);
    ystring_t *type_dmseen = ygeta_q(argc - 6, &ntot, dims);
    float *alt = ygeta_f(argc - 7, &ntot, dims);
    int ndm = ygets_i(argc - 8);

    carma_context *context_handle = _getCurrentContext();
    int activeDevice =
        context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (argc > 8) {
      long Nphi = ygets_l(argc - 9);
      rtc_handler->add_controller_geo(nactu, Nphi, delay, activeDevice,
                                      dms_handler, type_dmseen, alt, ndm);
    }
    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    else
      rtc_handler->add_controller(nactu, delay, activeDevice, type_control,
                                  dms_handler, type_dmseen, alt, ndm);
  }

  void Y_rtc_rmcontrol(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    rtc_handler->rm_controller();
  }

  void Y_rtc_setthresh(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);
    float thresh = ygets_f(argc - 3);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    if (rtc_handler->d_centro.at(ncentro)->is_type("tcog")) {
      CAST(sutra_centroider_tcog *, centroider_tcog,
           rtc_handler->d_centro.at(ncentro));
      centroider_tcog->set_threshold(thresh);
    }
  }

  void Y_rtc_setnmax(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);
    int nmax = ygets_f(argc - 3);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);
    if (rtc_handler->d_centro.at(ncentro)->is_type("bpcog")) {
      CAST(sutra_centroider_bpcog *, centroider_bpcog,
           rtc_handler->d_centro.at(ncentro));
      centroider_bpcog->set_nmax(nmax);
    }
  }

  void Y_rtc_docentroids(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    if (argc > 1) {
      long ncontrol = ygets_l(argc - 2);
      rtc_handler->do_centroids(ncontrol);
    } else
      rtc_handler->do_centroids();

    /* now done inside sutra_rtc::do_centroids
     sutra_context *context_handle = _getCurrentContext();
     int activeDevice = context_handle->set_activeDeviceForCpy(rhandler->device,1);
     */
  }

  void Y_rtc_docentroids_geom(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    if (argc > 1) {
      long ncontrol = ygets_l(argc - 2);
      rtc_handler->do_centroids_geom(ncontrol);
    } else
      rtc_handler->do_centroids(0);

    /* now done inside sutra_rtc::do_centroids
     sutra_context *context_handle = _getCurrentContext();
     int activeDevice = context_handle->set_activeDeviceForCpy(rhandler->device,1);
     */
  }

  void Y_rtc_doimat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);

    int geom = -1;
    if (argc > 3)
      geom = ygets_i(argc - 4);
    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (geom > -1)
      rtc_handler->do_imat_geom(ncontrol, dms_handler, geom);
    else
      rtc_handler->do_imat(ncontrol, dms_handler);
  }

  void Y_rtc_setgain(int argc) {
    // retreive pointer to persistent object pointing on c++ class
    // inside the yorick interpreter
    // this is the first expected element on the stack => argc-1
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    // retreive the pointer to sutra_rtc from the previous
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    // retreive the index of the desired controller in the vector of controllers
    long ncontrol = ygets_l(argc - 2);

    // get context and check if proper device is selected
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    // retreive gain and set it in the corresponding sutra_controller
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      float gain = ygets_f(argc - 3);
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->set_gain(gain);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      float gain = ygets_f(argc - 3);
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->set_gain(gain);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("cured") == 0) {
      float gain = ygets_f(argc - 3);
      CAST(sutra_controller_cured *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_gain(gain);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("geo") == 0) {
      float gain = ygets_f(argc - 3);
      CAST(sutra_controller_geo *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_gain(gain);
    }
    if ((rtc_handler->d_control.at(ncontrol)->get_type().compare("kalman_CPU")
        == 0)
        || (rtc_handler->d_control.at(ncontrol)->get_type().compare(
            "kalman_GPU")
            == 0)
        || (rtc_handler->d_control.at(ncontrol)->get_type().compare(
            "kalman_uninitialized")
            == 0)) {
      //Pour kalman gain correspond a k_W
      float gain = ygets_f(argc - 3);
      CAST(sutra_controller_kalman *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_gain(gain);
    }
  }

  void Y_rtc_getmgain(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_gain->getDims());
      control->d_gain->device2host(data);
    } else {
      y_error("Controller needs to be ls or mv\n");
    }
  }

  void Y_rtc_loadmgain(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      float *mgain = ygeta_f(argc - 3, &ntot, dims);
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->set_mgain(mgain);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv")
        == 0) {
      float *mgain = ygeta_f(argc - 3, &ntot, dims);
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->set_mgain(mgain);
    } else {
      y_error("Controller needs to be ls or mv\n");
    }
  }

  void Y_rtc_loadOpenLoopSlp(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      float *ol_slopes = ygeta_f(argc - 3, &ntot, dims);
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->loadOpenLoopSlp(ol_slopes);
    } else {
      y_error("Controller needs to be ls \n");
    }
  }

  void Y_rtc_loadnoisemat(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      float *noise = ygeta_f(argc - 3, &ntot, dims);
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->load_noisemat(noise);
    } else {
      y_error("Controller needs to be mv\n");
    }
  }

  void Y_rtc_getnoisemat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_noisemat->getDims());
      control->d_noisemat->device2host(data);
    } else {
      y_error("Controller needs to be mv\n");
    }
  }

  void Y_rtc_getslpol(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (control->is_modopti == 1) {
        float *data = ypush_f((long*) control->d_slpol->getDims());
        control->d_slpol->device2host(data);
      } else {
        y_error("Modal Optimization not initialized");
      }
    } else {
      y_error("Controller needs to be ls\n");
    }
  }

  void Y_rtc_getHcor(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (control->is_modopti == 1) {
        float *data = ypush_f((long*) control->d_Hcor->getDims());
        control->d_Hcor->device2host(data);
      } else {
        y_error("Modal Optimization not initialized");
      }
    } else {
      y_error("Controller needs to be ls\n");
    }
  }

  void Y_rtc_getS2M(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (control->is_modopti == 1) {
        float *data = ypush_f((long*) control->d_S2M->getDims());
        control->d_S2M->device2host(data);
      } else {
        y_error("Modal Optimization not initialized");
      }
    } else {
      y_error("Controller needs to be ls\n");
    }
  }

  void Y_rtc_getM2V(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (control->is_modopti == 1) {
        float *data = ypush_f((long*) control->d_M2V->getDims());
        control->d_M2V->device2host(data);
      } else {
        y_error("Modal Optimization not initialized");
      }
    } else {
      y_error("Controller needs to be ls\n");
    }
  }

  void Y_rtc_setdelay(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      float delay = ygets_f(argc - 3);
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->set_delay(delay);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      float delay = ygets_f(argc - 3);
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->set_delay(delay);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("cured") == 0) {
      float delay = ygets_f(argc - 3);
      CAST(sutra_controller_cured *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_delay(delay);
    }
  }

  void Y_rtc_getimat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_imat->getDims());
      control->d_imat->device2host(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv")
        == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_imat->getDims());
      control->d_imat->device2host(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("cured")
        == 0) {
      CAST(sutra_controller_cured *, control,
           rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_imat->getDims());
      control->d_imat->device2host(data);
    }
  }

  void Y_rtc_setimat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->d_imat->host2device(data);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_imat->host2device(data);
    }
  }

  void Y_rtc_setcentroids(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->d_centroids->host2device(data);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_centroids->host2device(data);
    }
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("cured") == 0) {
      CAST(sutra_controller_cured *, control,
           rtc_handler->d_control.at(ncontrol));
      control->d_centroids->host2device(data);
    }
  }

  void Y_rtc_getcentroids(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    sutra_wfs *wfs_handler = 0L;

    if (ncontrol >= abs(rtc_handler->d_control.size())) {
      if (argc < 4) {
        fprintf(
            stderr,
            "Controller not initialized on the GPU, you have to specify the WFS\n");
        return;
      }
      sensors_struct *shandler = (sensors_struct *) yget_obj(argc - 3,
                                                             &ySensors);
      sutra_sensors *sensor_handler = (sutra_sensors *) shandler->sutra_sensors;
      long nwfs = ygets_l(argc - 4);
      wfs_handler = sensor_handler->d_wfs.at(nwfs);
    }
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (ncontrol >= abs(rtc_handler->d_control.size())) {

      float *data = ypush_f((long*) wfs_handler->d_slopes->getDims());
      carma_obj<float> d_tmp(context_handle, wfs_handler->d_subsum->getDims());
      carma_obj<float> d_data(context_handle, wfs_handler->d_slopes->getDims());
      rtc_handler->d_centro.at(ncontrol)->get_cog(d_tmp, d_data);
      d_data.device2host(data);
    } else {
      float *data = ypush_f(
          (long*) rtc_handler->d_control.at(ncontrol)->d_centroids->getDims());
      rtc_handler->d_control.at(ncontrol)->d_centroids->device2host(data);
    }
  }

  void Y_rtc_getcom(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    float *data = ypush_f(
        (long*) rtc_handler->d_control.at(ncontrol)->d_com->getDims());
    rtc_handler->d_control.at(ncontrol)->d_com->device2host(data);
  }

  void Y_rtc_getcmat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_cmat->getDims());
      control->d_cmat->device2host(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv")
        == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_cmat->getDims());
      control->d_cmat->device2host(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare(
        "generic")
               == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      float *data = ypush_f((long*) control->d_cmat->getDims());
      control->d_cmat->device2host(data);
    } else {
      y_error("Controller needs to be ls or mv\n");
    }
  }
  void Y_rtc_setcmat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->d_cmat->host2device(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv")
        == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_cmat->host2device(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare(
        "generic")
               == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      control->d_cmat->host2device(data);
    } else {
      y_error("Controller needs to be ls or mv\n");
    }
  }
  void Y_rtc_setdecayFactor(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("generic") == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_decayFactor(data);
    } else
      y_error("Controller needs to be generic\n");
  }

  void Y_rtc_setmatE(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("generic") == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_matE(data);
    } else
      y_error("Controller needs to be generic\n");
  }

  void Y_rtc_setmgain(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("generic") == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_mgain(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls")
        == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->d_gain->host2device(data);
    } else if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls")
        == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_gain->host2device(data);
    } else
      y_error("Controller needs to be generic, ls or mv \n");
  }

  void Y_rtc_setcommandlaw(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    ystring_t law = ygets_q(argc - 3);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("generic") == 0) {
      CAST(sutra_controller_generic *, control,
           rtc_handler->d_control.at(ncontrol));
      control->set_commandlaw(law);
    } else
      y_error("Controller needs to be generic\n");
  }

  void Y_rtc_setcovmat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_covmat->host2device(data);
    } else {
      y_error("Controller needs to be mv\n");
    }
  }
  void Y_rtc_setCphim(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_Cphim->host2device(data);
    } else {
      y_error("Controller needs to be mv\n");
    }
  }

  void Y_rtc_setCmm(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    long ntot;
    long dims[Y_DIMSIZE];
    float *data = ygeta_f(argc - 3, &ntot, dims);
    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("mv") == 0) {
      CAST(sutra_controller_mv *, control, rtc_handler->d_control.at(ncontrol));
      control->d_Cmm->host2device(data);
    } else {
      y_error("Controller needs to be mv\n");
    }
  }

  void Y_rtc_setperturbcom(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float *perturb = ygeta_f(argc - 3, &ntot, dims);
    int dim = ygets_i(argc - 4);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    rtc_handler->d_control.at(ncontrol)->set_perturbcom(perturb, dim);

  }

  void Y_rtc_setopenloop(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    int openloop = ygets_i(argc - 3);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    rtc_handler->d_control.at(ncontrol)->set_openloop(openloop);

  }

  void Y_rtc_buildcmat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    long nfilt = ygets_l(argc - 3);

    long filt_tt = 0;
    if (argc > 3)
      filt_tt = ygets_l(argc - 4);
    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (filt_tt > 0)
        control->build_cmat(nfilt, true);
      else
        control->build_cmat(nfilt);
    }
  }

  void Y_rtc_init_proj(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);
    int *indx_dm = ygeta_i(argc - 4, &ntot, dims);
    float *unitpervolt = ygeta_f(argc - 5, &ntot, dims);
    int *indx_pup = ygeta_i(argc - 6, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("geo") == 0) {
      CAST(sutra_controller_geo *, control,
           rtc_handler->d_control.at(ncontrol));
      control->init_proj_sparse(dms_handler, indx_dm, unitpervolt, indx_pup);
    } else
      y_error("**** ERROR : init_proj only for controller type geo **** \n");
  }

  void Y_rtc_initModalOpti(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float nmodes = ygets_f(argc - 3);
    float nrec = ygets_f(argc - 4);
    float *M2V = ygeta_f(argc - 5, &ntot, dims);
    float gmin = ygets_f(argc - 6);
    float gmax = ygets_f(argc - 7);
    int ngain = ygets_i(argc - 8);
    float Fs = ygets_f(argc - 9);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->init_modalOpti(nmodes, nrec, M2V, gmin, gmax, ngain, Fs);
    } else
      y_error(
          "**** ERROR : Modal Optimization only for controller type ls **** \n");

  }

  void Y_rtc_modalControlOptimization(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      if (control->is_modopti)
        control->modalControlOptimization();
      else
        y_error("**** ERROR : Modal Optimization not initialized **** \n");
    } else
      y_error(
          "**** ERROR : Modal Optimization only for controller type ls **** \n");
  }

  void Y_rtc_imatsvd(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));

      if (control->svdec_imat() == EXIT_FAILURE)
        y_error(
            "***** ERROR : sutra controller has no SVD implementation *****\n");
    }
  }

  void Y_rtc_framedelay(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    int ncontrol = ygets_i(argc - 2);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
      CAST(sutra_controller_ls *, control, rtc_handler->d_control.at(ncontrol));
      control->frame_delay();
    }
  }

  void Y_rtc_applycontrol(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    rtc_handler->apply_control(ncontrol, dms_handler);

  }

  void Y_rtc_docontrol(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    rtc_handler->do_control(ncontrol);

  }

  void Y_rtc_docontrol_geo(int argc) {
    //long ntot;
    //long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);
    target_struct *handler = (target_struct *) yget_obj(argc - 4, &yTarget);
    sutra_target *target_handler = (sutra_target *) (handler->sutra_target);
    long ntarget = ygets_l(argc - 5);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    CAST(sutra_controller_geo *, control, rtc_handler->d_control.at(ncontrol));
    control->comp_dphi(target_handler->d_targets.at(ntarget));
    rtc_handler->do_control(ncontrol);
  }

  void Y_controller_setdata(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    int ncontrol = ygets_i(argc - 2);

    char *type_data = ygets_q(argc - 3);
    long ntot;
    long dims[Y_DIMSIZE];

    if (strcmp(type_data, "U") == 0) {
      float *data = ygeta_f(argc - 4, &ntot, dims);
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        control->d_U->host2device(data);
      }
    } else if (strcmp(type_data, "eigenvals") == 0) {
      float *data = ygeta_f(argc - 4, &ntot, dims);
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        control->h_eigenvals->fill_from(data);
        control->d_eigenvals->host2device(data);
      }
    }
  }

  void Y_controller_getdata(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    int ncontrol = ygets_i(argc - 2);

    char *type_data = ygets_q(argc - 3);
    if (strcmp(type_data, "U") == 0) {
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        float *data = ypush_f((long*) control->d_U->getDims());
        control->d_U->device2host(data);
      }
    } else if (strcmp(type_data, "eigenvals") == 0) {
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        float *data = ypush_f((long*) control->h_eigenvals->getDims());
        control->h_eigenvals->fill_into(data);
      }
    } else if (strcmp(type_data, "cenbuff") == 0) {
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        float *data = ypush_f((long*) control->d_cenbuff->getDims());
        control->d_cenbuff->device2host(data);
      }
    } else if (strcmp(type_data, "err") == 0) {
      if (rtc_handler->d_control.at(ncontrol)->get_type().compare("ls") == 0) {
        CAST(sutra_controller_ls *, control,
             rtc_handler->d_control.at(ncontrol));
        float *data = ypush_f((long*) control->d_err->getDims());
        control->d_err->device2host(data);
      }
    } else if (strcmp(type_data, "com") == 0) {
      float *data = ypush_f(
          (long*) rtc_handler->d_control.at(ncontrol)->d_com->getDims());
      rtc_handler->d_control.at(ncontrol)->d_com->device2host(data);
    } else if (strcmp(type_data, "centroids") == 0) {
      float *data = ypush_f(
          (long*) rtc_handler->d_control.at(ncontrol)->d_centroids->getDims());
      rtc_handler->d_control.at(ncontrol)->d_centroids->device2host(data);
    }
  }

  void Y_controller_initcured(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    int nxsubs = ygets_l(argc - 3);
    int *isvalid = ygeta_i(argc - 4, &ntot, dims);
    int ndivs = ygets_l(argc - 5);
    int tt = ygets_l(argc - 6);

    //rtc_struct *handle    =(rtc_struct *)ypush_obj(&yRTC, sizeof(rtc_struct));

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    if (rtc_handler->d_control.at(ncontrol)->get_type().compare("cured") == 0) {
      CAST(sutra_controller_cured *, control,
           rtc_handler->d_control.at(ncontrol));
      control->init_cured(nxsubs, isvalid, ndivs, tt);
    }
  }

  /*
   *     ____  _
   *    / ___|| | ___  _ __   ___  ___
   *    \___ \| |/ _ \| '_ \ / _ \/ __|
   *     ___) | | (_) | |_) |  __/\__ \
   *    |____/|_|\___/| .__/ \___||___/
   *                  |_|
   *
   */

  void Y_slopes_geom(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    int nsensor = ygets_i(argc - 2);

    int type = ygets_i(argc - 3);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    if (sensors_handler->d_wfs.at(nsensor)->type != "sh") {
      y_error("wfs should be a SH");
    }
    sutra_wfs_sh *wfs = dynamic_cast<sutra_wfs_sh *>(sensors_handler->d_wfs.at(
        nsensor));
    wfs->slopes_geom(type);
  }

  void Y_sensors_getslopes(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    float *data = ypush_f(
        (long*) sensors_handler->d_wfs.at(nsensor)->d_slopes->getDims());
    sensors_handler->d_wfs.at(nsensor)->d_slopes->device2host(data);
  }

  void Y_sensors_initnmax(int argc) {
    /*
     sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
     sutra_sensors *sensors_handler = (sutra_sensors *)(handler->sutra_sensors);

     carma_context *context_handle = _getCurrentContext();
     int activeDevice = context_handle->set_activeDevice(handler->device,1);

     int nsensor = ygets_i(argc-2);

     int nmax = ygets_i(argc-3);

     //sensors_handler->d_wfs.at(nsensor)->init_nmax(nmax);
     */
  }

  void Y_sensors_initweights(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);

    float *weights = ygeta_f(argc - 3, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    CAST(sutra_centroider_wcog *, centroider_wcog,
         rtc_handler->d_centro.at(ncentro));
    centroider_wcog->init_weights();
    centroider_wcog->load_weights(weights, dims[0]);
  }

  void Y_sensors_initbcube(int argc) {

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    int ncentro = ygets_i(argc - 2);

    if (rtc_handler->d_centro.at(ncentro)->get_type() == "corr") {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      centroider_corr->init_bincube();
    }
  }

  void Y_sensors_initcorr(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);

    float *weights = ygeta_f(argc - 3, &ntot, dims);
    int mydim = dims[0];
    float *corr_norm = ygeta_f(argc - 4, &ntot, dims);
    int sizex = ygets_i(argc - 5);
    int sizey = ygets_i(argc - 6);
    float *interpmat = ygeta_f(argc - 7, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    CAST(sutra_centroider_corr *, centroider_corr,
         rtc_handler->d_centro.at(ncentro));
    centroider_corr->init_corr(sizex, sizey, interpmat);
    centroider_corr->load_corr(weights, corr_norm, mydim);
  }

  void Y_sensors_loadcorrfnct(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);

    float *weights = ygeta_f(argc - 3, &ntot, dims);
    int mydim = dims[0];
    float *corr_norm = ygeta_f(argc - 4, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    CAST(sutra_centroider_corr *, centroider_corr,
         rtc_handler->d_centro.at(ncentro));
    centroider_corr->load_corr(weights, corr_norm, mydim);
  }

  void Y_sensors_loadweights(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 2, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 3);

    float *weights = ygeta_f(argc - 4, &ntot, dims);

    CAST(sutra_centroider_wcog *, centroider_wcog,
         rtc_handler->d_centro.at(ncentro));
    centroider_wcog->load_weights(weights, dims[0]);
  }

  void Y_sensors_compslopes(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    int ncentro = ygets_i(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(rhandler->device,1);

    //cout << ncentro << " " << (rtc_handler->d_centro.at(ncentro)->typec) << endl;
    if (argc > 4) {
      if (rtc_handler->d_centro.at(ncentro)->is_type("bpcog")) {
        CAST(sutra_centroider_bpcog *, centroider_bpcog,
             rtc_handler->d_centro.at(ncentro));
        centroider_bpcog->set_nmax(ygets_i(argc - 5) );
      } else if (rtc_handler->d_centro.at(ncentro)->is_type("tcog")) {
        CAST(sutra_centroider_tcog *, centroider_tcog,
             rtc_handler->d_centro.at(ncentro));
        centroider_tcog->set_threshold(ygets_f(argc - 5) );
      }

    }

    rtc_handler->d_centro.at(ncentro)->get_cog();

  }

  void Y_sensors_getnmax(int argc) {
    /*
     sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
     sutra_sensors *sensors_handler = (sutra_sensors *)(handler->sutra_sensors);

     carma_context *context_handle = _getCurrentContext();
     int activeDevice = context_handle->set_activeDeviceForCpy(handler->device,1);

     int nsensor = ygets_i(argc-2);

     //sensors_handler->d_wfs.at(nsensor)->get_nmax();
     */
  }

  void Y_centroider_getdata(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    int ncentro = ygets_i(argc - 2);

    char *type_data = ygets_q(argc - 3);
    if ((strcmp(type_data, "weights") == 0) && rtc_handler->d_centro.at(ncentro)->is_type(
        "wcog")) {
      CAST(sutra_centroider_wcog *, centroider_wcog,
           rtc_handler->d_centro.at(ncentro));
      float *data = ypush_f((long*) centroider_wcog->d_weights->getDims());
      centroider_wcog->d_weights->device2host(data);
    }
    if ((strcmp(type_data, "corrnorm") == 0) && rtc_handler->d_centro.at(
        ncentro)->is_type("corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      float *data = ypush_f((long*) centroider_corr->d_corrnorm->getDims());

      centroider_corr->d_corrnorm->device2host(data);
    }
    if ((strcmp(type_data, "corrfnct") == 0) && rtc_handler->d_centro.at(
        ncentro)->is_type("corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      long *ndims_data = new long[5];
      ndims_data[0] = 4;
      ndims_data[1] = 2;
      const long *ndims_obj = centroider_corr->d_corrfnct->getDims();
      memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
      float *data = ypush_f(ndims_data);
      centroider_corr->d_corrfnct->device2host((cuFloatComplex *) data);
    }
    if ((strcmp(type_data, "corrspot") == 0) && rtc_handler->d_centro.at(
        ncentro)->is_type("corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      long *ndims_data = new long[5];
      ndims_data[0] = 4;
      ndims_data[1] = 2;
      const long *ndims_obj = centroider_corr->d_corrspot->getDims();
      memcpy(&ndims_data[2], &(ndims_obj[1]), sizeof(long) * 3);
      float *data = ypush_f(ndims_data);
      centroider_corr->d_corrspot->device2host((cuFloatComplex *) data);
    }
    if ((strcmp(type_data, "corr") == 0) && rtc_handler->d_centro.at(ncentro)->is_type(
        "corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      float *data = ypush_f((long*) centroider_corr->d_corr->getDims());
      centroider_corr->d_corr->device2host(data);
    }
    if ((strcmp(type_data, "corrmax") == 0) && rtc_handler->d_centro.at(ncentro)->is_type(
        "corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      int *data = ypush_i((long*) centroider_corr->d_corrmax->getDims());
      centroider_corr->d_corrmax->device2host(data);
    }
    if ((strcmp(type_data, "matinterp") == 0) && rtc_handler->d_centro.at(
        ncentro)->is_type("corr")) {
      CAST(sutra_centroider_corr *, centroider_corr,
           rtc_handler->d_centro.at(ncentro));
      float *data = ypush_f((long*) centroider_corr->d_interpmat->getDims());
      centroider_corr->d_interpmat->device2host(data);
    }
  }

  /*
   *      ____ _       _           _
   *     / ___| | ___ | |__   __ _| |
   *    | |  _| |/ _ \| '_ \ / _` | |
   *    | |_| | | (_) | |_) | (_| | |
   *     \____|_|\___/|_.__/ \__,_|_|
   *
   */

  void Y_move_atmos(int argc) {
    atmos_struct *handler_a = (atmos_struct *) yget_obj(argc - 1, &yAtmos);
    sutra_atmos *atmos_handler = (sutra_atmos *) (handler_a->sutra_atmos);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(handler_a->device,1);

    atmos_handler->move_atmos();
  }

  void Y_sensors_trace(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForce(handler->device,1);

    int nsensor = ygets_i(argc - 2);

    char *type_trace = ygets_q(argc - 3);
    if (strcmp(type_trace, "all") == 0) {
      atmos_struct *handlera = (atmos_struct *) yget_obj(argc - 4, &yAtmos);
      sutra_atmos *atmos_handler = (sutra_atmos *) (handlera->sutra_atmos);
      dms_struct *handlerd = (dms_struct *) yget_obj(argc - 5, &yDMs);
      sutra_dms *dms_handler = (sutra_dms *) (handlerd->sutra_dms);
      sensors_handler->d_wfs.at(nsensor)->sensor_trace(atmos_handler,
                                                       dms_handler);
    } else if (strcmp(type_trace, "atmos") == 0) {
      atmos_struct *handlera = (atmos_struct *) yget_obj(argc - 4, &yAtmos);
      sutra_atmos *atmos_handler = (sutra_atmos *) (handlera->sutra_atmos);
      sensors_handler->d_wfs.at(nsensor)->sensor_trace(atmos_handler);
    } else if (strcmp(type_trace, "dm") == 0) {
      dms_struct *handlera = (dms_struct *) yget_obj(argc - 4, &yDMs);
      sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);
      int rst;
      if (argc > 4)
        rst = ygets_i(argc - 5);
      else
        rst = 0;
      sensors_handler->d_wfs.at(nsensor)->sensor_trace(dms_handler, rst);
    }

  }

  void Y_yoga_add_telemetry_obj(int argc) {
    telemetry_struct *handler = (telemetry_struct *) yget_obj(argc - 1,
                                                              &yTelemetry);
    sutra_telemetry *telemetry_handler =
        (sutra_telemetry *) (handler->sutra_telemetry);

    string type_obj = string(ygets_q(argc - 2));

    if (type_obj.compare("binimg") == 0) {
      sensors_struct *handler2 = (sensors_struct *) yget_obj(argc - 3,
                                                             &ySensors);
      sutra_sensors *sensors_handler =
          (sutra_sensors *) (handler2->sutra_sensors);

      for (size_t i = 0; i < sensors_handler->nsensors(); i++) {
        carma_host_obj<float> *yho = new carma_host_obj<float>(
            sensors_handler->d_wfs[i]->d_binimg->getDims(), MA_PAGELOCK, 1);
        telemetry_handler->add_obj(type_obj, i, yho);
      }
    }
  }

  void Y_yoga_add_telemetry_stream(int argc) {

    telemetry_struct *handler = (telemetry_struct *) yget_obj(argc - 1,
                                                              &yTelemetry);
    sutra_telemetry *telemetry_handler =
        (sutra_telemetry *) (handler->sutra_telemetry);

    int nb_stream = 1;
    if (argc > 1)
      nb_stream = ygets_i(argc - 2);

    telemetry_handler->add_stream(nb_stream);
  }

  /*
   *               _                       _       _
   *  __ _  ___   | |_ ___ _ __ ___  _ __ | | __ _| |_ ___
   * / _` |/ _ \  | __/ _ \ '_ ` _ \| '_ \| |/ _` | __/ _ \
   *| (_| | (_) | | ||  __/ | | | | | |_) | | (_| | ||  __/
   * \__,_|\___/   \__\___|_| |_| |_| .__/|_|\__,_|\__\___|
   *                                |_|
   *
   */

  aotemplate_struct*
  yoga_ao_getyAotemplate(int argc, int pos) {
    return (aotemplate_struct *) yget_obj(argc - pos, &yAotemplate);
  }

  void aotemplate_free(void *obj) {
    aotemplate_struct *handler = (aotemplate_struct *) obj;
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);
    try {
      sutra_aotemplate *aotemplate_obj_handler =
          (sutra_aotemplate *) (handler->sutra_aotemplate);
      delete aotemplate_obj_handler;
    } catch (string& msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    }
  }

  void aotemplate_print(void *obj) {
    aotemplate_struct *handler = (aotemplate_struct *) obj;
    sutra_aotemplate *aotemplate_handler =
        (sutra_aotemplate *) (handler->sutra_aotemplate);
    cout << "Yoga Aotemplate Object" << endl;

    cout << "Contains " << aotemplate_handler->d_data->getNbElem()
         << " elements " << endl;
  }

  void Y_yoga_aotemplate(int argc) {
    try {
      carma_context *context_handle = _getCurrentContext();
      int activeDevice = context_handle->get_activeDevice();
      int odevice = activeDevice;
      //char *type            = ygets_q(argc-2);
      long dim = ygets_l(argc - 1);
      if (argc > 1)
        odevice = ygets_i(argc - 2);

      activeDevice = context_handle->set_activeDevice(odevice,1);

      aotemplate_struct *handle = (aotemplate_struct *) ypush_obj(
          &yAotemplate, sizeof(aotemplate_struct));
      handle->device = activeDevice;
      handle->sutra_aotemplate = new sutra_aotemplate(context_handle, "test",
                                                      dim, odevice);

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_aotemplate construction in " << __FILE__
          << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_yoga_templatefill(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    aotemplate_struct *handler = (aotemplate_struct *) yget_obj(argc - 1,
                                                                &yAotemplate);
    sutra_aotemplate *aotemplate_handler =
        (sutra_aotemplate *) (handler->sutra_aotemplate);
    float *data;
    if (argc > 1) {
      data = ygeta_f(argc - 2, &ntot, dims);
      aotemplate_handler->fill_data(data);
    } else
      aotemplate_handler->fill_data();
  }

  void Y_yoga_templatecomp(int argc) {
    aotemplate_struct *handler = (aotemplate_struct *) yget_obj(argc - 1,
                                                                &yAotemplate);
    sutra_aotemplate *aotemplate_handler =
        (sutra_aotemplate *) (handler->sutra_aotemplate);
    aotemplate_handler->do_compute();
  }

  void Y_yoga_gettemplate(int argc) {
    aotemplate_struct *handler = (aotemplate_struct *) yget_obj(argc - 1,
                                                                &yAotemplate);
    sutra_aotemplate *aotemplate_handler =
        (sutra_aotemplate *) (handler->sutra_aotemplate);

    if (argc > 1) {
      char *type_data = ygets_q(argc - 2);
      if (strcmp(type_data, "data") == 0) {
        float *data = ypush_f((long*) aotemplate_handler->d_data->getDims());
        aotemplate_handler->d_data->device2host(data);
      }
      if (strcmp(type_data, "res") == 0) {
        float *data = ypush_f((long*) aotemplate_handler->d_res->getDims());
        aotemplate_handler->d_res->device2host(data);
      }
    } else {
      float *data = ypush_f((long*) aotemplate_handler->d_data->getDims());
      aotemplate_handler->d_data->device2host(data);
    }
  }
// Florian features
  void Y_rtc_buildcmatmv(int argc) {
//  long ntot;
//  long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float cond = ygets_f(argc - 3);
    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    controller->build_cmat(cond);
  }
  void Y_rtc_doimatkl(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    rtc_handler->do_imatkl(ncontrol, dms_handler);
  }
  void Y_rtc_doimatkl4pzt(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handlera = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handlera->sutra_dms);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    rtc_handler->do_imatkl4pzt(ncontrol, dms_handler);
  }

  // Florian features
  void Y_rtc_doCmm(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    sensors_struct *shandler = (sensors_struct *) yget_obj(argc - 3, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (shandler->sutra_sensors);
    atmos_struct *ahandler = (atmos_struct *) yget_obj(argc - 4, &yAtmos);
    sutra_atmos *atmos_handler = (sutra_atmos *) (ahandler->sutra_atmos);
    float diamTel = ygets_f(argc - 5);
    float cobs = ygets_f(argc - 6);
    float *L0 = ygeta_f(argc - 7, &ntot, dims);
    float *cn2 = ygeta_f(argc - 8, &ntot, dims);
    float *alphaX = ygeta_f(argc - 9, &ntot, dims);
    float *alphaY = ygeta_f(argc - 10, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    controller->compute_Cmm(atmos_handler, sensors_handler, L0, cn2, alphaX,
                            alphaY, diamTel, cobs);

  }

  void Y_rtc_doCphim(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    sensors_struct *shandler = (sensors_struct *) yget_obj(argc - 3, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (shandler->sutra_sensors);
    atmos_struct *ahandler = (atmos_struct *) yget_obj(argc - 4, &yAtmos);
    sutra_atmos *atmos_handler = (sutra_atmos *) (ahandler->sutra_atmos);
    dms_struct *handler = (dms_struct *) yget_obj(argc - 5, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    float *L0 = ygeta_f(argc - 6, &ntot, dims);
    float *cn2 = ygeta_f(argc - 7, &ntot, dims);
    float *alphaX = ygeta_f(argc - 8, &ntot, dims);
    float *alphaY = ygeta_f(argc - 9, &ntot, dims);
    float *X = ygeta_f(argc - 10, &ntot, dims);
    float *Y = ygeta_f(argc - 11, &ntot, dims);
    float *xactu = ygeta_f(argc - 12, &ntot, dims);
    float *yactu = ygeta_f(argc - 13, &ntot, dims);
    float diamTel = ygets_f(argc - 14);
    float k2 = ygets_f(argc - 15);
    float *Nact = ygeta_f(argc - 16, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    controller->compute_Cphim(atmos_handler, sensors_handler, dms_handler, L0,
                              cn2, alphaX, alphaY, X, Y, xactu, yactu, diamTel,
                              k2, Nact);

  }

  void Y_rtc_docovmat(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    dms_struct *handler = (dms_struct *) yget_obj(argc - 3, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    char *type = ygets_q(argc - 4);
    float alt = ygets_f(argc - 5);
    int *indx_pup = ygeta_i(argc - 6, &ntot, dims);
    long dim = ygets_l(argc - 7);
    float *xpos = ygeta_f(argc - 8, &ntot, dims);
    float *ypos = ygeta_f(argc - 9, &ntot, dims);
    long Nkl = ygets_l(argc - 10);
    float norm = ygets_f(argc - 11);
    float ampli = ygets_f(argc - 12);
    char *method = ygets_q(argc - 13);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);
    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);

    int inddm = dms_handler->get_inddm(type, alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_docovmat function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << type << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }
    controller->do_covmat(dms_handler->d_dms[inddm], method, indx_pup, dim,
                          xpos, ypos, Nkl, norm, ampli);

  }
  void Y_rtc_loadcovmat(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float *covmat = ygeta_f(argc - 3, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    controller->load_covmat(covmat);

  }
  void Y_rtc_loadklbasis(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float *klbasis = ygeta_f(argc - 3, &ntot, dims);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    controller->load_klbasis(klbasis);

  }
  void Y_rtc_getvoltage(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    float *data = ypush_f(
        (long*) rtc_handler->d_control[ncontrol]->d_voltage->getDims());
    rtc_handler->d_control[ncontrol]->d_voltage->device2host(data);
  }

  void Y_rtc_getproj(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_geo *, controller, rtc_handler->d_control[ncontrol]);
    float *data = ypush_f((long*) controller->d_proj->getDims());
    controller->d_proj->device2host(data);
  }

  void Y_rtc_getphi(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_geo *, controller, rtc_handler->d_control[ncontrol]);
    double *data = (double*) ypush_f((long*) controller->d_phi->getDims());
    controller->d_phi->device2host(data);
  }

  void Y_rtc_getCmm(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    float *data = ypush_f((long*) controller->d_Cmm->getDims());
    controller->d_Cmm->device2host(data);
  }
  void Y_rtc_getCphim(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    float *data = ypush_f((long*) controller->d_Cphim->getDims());
    controller->d_Cphim->device2host(data);
  }
  void Y_rtc_getcovmat(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    float *data = ypush_f((long*) controller->d_covmat->getDims());
    controller->d_covmat->device2host(data);
  }
  void Y_rtc_getklbasis(int argc) {
    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDeviceForCpy(rhandler->device,1);

    CAST(sutra_controller_mv *, controller, rtc_handler->d_control[ncontrol]);
    float *data = ypush_f((long*) controller->d_KLbasis->getDims());
    controller->d_KLbasis->device2host(data);
  }

  void Y_sensors_rmlayer(int argc) {
    sensors_struct *handler = (sensors_struct *) yget_obj(argc - 1, &ySensors);
    sutra_sensors *sensors_handler = (sutra_sensors *) (handler->sutra_sensors);

    carma_context *context_handle = _getCurrentContext();
    context_handle->set_activeDevice(handler->device,1);

    int nsensor = ygets_i(argc-2);
    char *type = ygets_q(argc - 3);
    float alt = ygets_f(argc-4);

    sensors_handler->d_wfs.at(nsensor)->d_gs->remove_layer(type, alt);
  }

// Florian features
  void Y_yoga_getflokl(int argc) {
    //long ntot;
    //long dims[Y_DIMSIZE];
    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    //char *type            = ygets_q(argc-2);
//  long nkl = ygets_l(argc - 2);
//  long dim = ygets_l(argc - 3);
//  float *covmat = ygeta_f(argc - 4, &ntot, dims);
//  float *filter = ygeta_f(argc - 5, &ntot, dims);
//  float *evals = ygeta_f(argc - 6, &ntot, dims);
//  float *bas = ygeta_f(argc - 7, &ntot, dims);
    float alt = ygets_f(argc-8);
    int inddm = dms_handler->get_inddm("kl", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getflokl function in " << __FILE__ << "@"
          << __LINE__ << " DM(" << "kl" << "," << alt << ") doesn't exist"
          << endl;
      y_error(buf.str().c_str());
    }

    cout << "flag 1" << endl;
    dms_handler->d_dms[inddm]->d_kl->get_flokl();
    cout << "flag 2" << endl;
    float *data = ypush_f(
        (long*) dms_handler->d_dms[inddm]->d_kl->d_bas->getDims());
    cout << "flag 3" << endl;
    dms_handler->d_dms[inddm]->d_kl->d_bas->device2host(data);
    cout << "flag 4" << endl;
  }

  void Y_yoga_floloadkl(int argc) {
    cout << "flag0" << endl;
    long ntot;
    long dims[Y_DIMSIZE];
    cout << "flag1" << endl;

    dms_struct *handler = (dms_struct *) yget_obj(argc - 1, &yDMs);
    cout << "flag2" << endl;
    sutra_dms *dms_handler = (sutra_dms *) (handler->sutra_dms);
    cout << "flag3" << endl;
    //char *type            = ygets_q(argc-2);
    //long  nkl             = ygets_l(argc-2);
    //cout << "flag4"<< endl;
    //long dim              = ygets_l(argc-3);
    //cout << "flag5"<< endl;
    float *covmat = ygeta_f(argc - 4, &ntot, dims);
    cout << "flag6" << endl;
    float *filter = ygeta_f(argc - 5, &ntot, dims);
    cout << "flag7" << endl;
    float *evals = ygeta_f(argc - 6, &ntot, dims);
    cout << "flag8" << endl;
    float *bas = ygeta_f(argc - 7, &ntot, dims);
    cout << "flag9" << endl;
    float alt = ygets_f(argc-8);
    cout << "flag10" << endl;
    int inddm = dms_handler->get_inddm("kl", alt);
    if (inddm < 0) {
      stringstream buf;
      buf << "unknown error with yoga_getflolaodkl function in " << __FILE__
          << "@" << __LINE__ << " DM(" << "kl" << "," << alt
          << ") doesn't exist" << endl;
      y_error(buf.str().c_str());
    }

    dms_handler->d_dms[inddm]->kl_floloadarrays(covmat, filter, evals, bas);
    cout << "flag11" << endl;
  }

//-----------------------------------------------------------------------------------------

  void Y_rtc_initkalman(int argc) {

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);
    float bruit = ygets_f(argc-3);
    long ntot, dims[Y_DIMSIZE];

    float* hD_Mo = ygeta_f(argc - 4, &ntot, dims);
    carma_host_obj<float> D_Mo(dims, hD_Mo);

    float* hN_Act = ygeta_f(argc - 5, &ntot, dims);
    carma_host_obj<float> N_Act(dims, hN_Act);

    float* hPROJ = ygeta_f(argc - 6, &ntot, dims);
    carma_host_obj<float> PROJ(dims, hPROJ);

    float* hSigmaV = ygeta_f(argc - 7, &ntot, dims);
    carma_host_obj<float> SigmaV(dims, hSigmaV);

    float* hatur = ygeta_f(argc - 8, &ntot, dims);
    carma_host_obj<float> atur(dims, hatur);

    float* hbtur = ygeta_f(argc - 9, &ntot, dims);
    carma_host_obj<float> btur(dims, hbtur);

    bool is_zonal = ygets_s(argc - 10);
    bool is_sparse = ygets_s(argc - 11);
    bool is_GPU = ygets_s(argc - 12);

    CAST(sutra_controller_kalman *, control,
         rtc_handler->d_control.at(ncontrol));
    if ((rtc_handler->d_control.at(ncontrol)->get_type().compare(
        "kalman_uninitialized")
         == 0)) {

      control->init_kalman(D_Mo, N_Act, PROJ, is_zonal, is_sparse, is_GPU);
    } else
      y_error("Controller needs to be either kalman_CPU or kalman_CPU\n");

    if ((rtc_handler->d_control.at(ncontrol)->get_type().compare("kalman_CPU")
        == 0)
        || (rtc_handler->d_control.at(ncontrol)->get_type().compare(
            "kalman_GPU")
            == 0)) {

      //bruit en arcsec^2
      control->calculate_gain(bruit, SigmaV, atur, btur);
    } else
      y_error(
          "Controller needs to be either kalman_CPU or Kalman_GPU.\nEnsure that COMPILATION_LAM was defined as an environment variable before compiling.\n");
  }

  void Y_rtc_kalmancalculategain(int argc) {

    long ntot, dims[Y_DIMSIZE];

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    float bruit = ygets_f(argc - 3);

    float* hSigmaV = ygeta_f(argc - 4, &ntot, dims);
    carma_host_obj<float> SigmaV(dims, hSigmaV);

    float* hatur = ygeta_f(argc - 5, &ntot, dims);
    carma_host_obj<float> atur(dims, hatur);

    float* hbtur = ygeta_f(argc - 6, &ntot, dims);
    carma_host_obj<float> btur(dims, hbtur);

    if ((rtc_handler->d_control.at(ncontrol)->get_type().compare("kalman_CPU")
        == 0)
        || (rtc_handler->d_control.at(ncontrol)->get_type().compare(
            "kalman_GPU")
            == 0)) {

      CAST(sutra_controller_kalman *, control,
           rtc_handler->d_control.at(ncontrol));
      control->calculate_gain(bruit, SigmaV, atur, btur);
    } else
      y_error("Controller needs to be either kalman_CPU or kalman_GPU\n");
  }

  void Y_rtc_kalmangettime(int argc) {

    rtc_struct *rhandler = (rtc_struct *) yget_obj(argc - 1, &yRTC);
    sutra_rtc *rtc_handler = (sutra_rtc *) (rhandler->sutra_rtc);
    long ncontrol = ygets_l(argc - 2);

    if ((rtc_handler->d_control.at(ncontrol)->get_type().compare("kalman_CPU")
        == 0)
        || (rtc_handler->d_control.at(ncontrol)->get_type().compare(
            "kalman_GPU")
            == 0)) {

      CAST(sutra_controller_kalman *, control,
           rtc_handler->d_control.at(ncontrol));
      double tps = control->gettime();
      double tps1 = control->gettime_op1();
      double tps2 = control->gettime_op2();
      double tps3 = control->gettime_op3();

      if (tps == 0)
        cout
            << "WARNING ! Elapsed time in Kalman is 0.00 s. Make sure that COMPILATION_LAM (the environment variable) contained the word \"temps_kalman\" when compiling. This is necessary to use rtc_kalmangettime."
            << endl;
      else
        cout << "temps boucle OA kalman = " << tps << " s" << endl;
      cout << "      - temps boucle OA kalman op1 = " << tps1 << " s" << endl;
      cout << "      - temps boucle OA kalman op2 = " << tps2 << " s" << endl;
      cout << "      - temps boucle OA kalman op3 = " << tps3 << " s" << endl;

    }
  }

}

