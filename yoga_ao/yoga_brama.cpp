#include <sutra_rtc_brama.h>
#include <sutra_target_brama.h>
#include <yoga_api.h>
#include <yoga_ao_api.h>

extern "C" {

  void Y_yoga_rtc_brama(int argc) {
    try {
      if (argc == 1) {
        rtc_struct* handle = yoga_ao_getyRTC(argc, 1);
        SCAST(sutra_rtc*, handle_rtc, handle->sutra_rtc);
        delete handle_rtc;

        carma_context *context_handle = _getCurrentContext();
        char brama_name[] = "yoga_rtc_brama";
        handle->sutra_rtc = new sutra_rtc_brama(context_handle, brama_name);
        handle->use_brama = 1;
      }
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

  void Y_rtc_publish(int argc) {
    try {
      rtc_struct* handle = yoga_ao_getyRTC(argc, 1);
      if (handle->use_brama == 0)
        y_error("DDS is not initialized");

      sutra_rtc_brama *rtc_handler = (sutra_rtc_brama *) (handle->sutra_rtc);

      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDevice(handle->device,1);

      rtc_handler->publish();
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with rtc_initDDS in " << __FILE__ << "@" << __LINE__
          << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_yoga_target_brama(int argc) {
    long ntot;
    long dims[Y_DIMSIZE];

    try {
      target_struct* handle = yoga_ao_getyTarget(argc, 1);
      SCAST(sutra_target*, handle_target, handle->sutra_target);
      delete handle_target;

      char brama_name[] = "yoga_target_brama";

      carma_context *context_handle = _getCurrentContext();
      int odevice = context_handle->get_activeDevice();

      int ntargets = ygets_i(argc - 2);

      float *xpos;
      xpos = ygeta_f(argc - 3, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for xpos");

      float *ypos;
      ypos = ygeta_f(argc - 4, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens ypos");

      float *lambda;
      lambda = ygeta_f(argc - 5, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens lambda");

      float *mag;
      mag = ygeta_f(argc - 6, &ntot, dims);
      if (ntot != ntargets)
        y_error("wrong dimension for screens mag");

      float zerop;
      zerop = ygets_f(argc - 7);

      long *sizes;
      sizes = ygeta_l(argc - 8, &ntot, dims);

      float *pup;
      pup = ygeta_f(argc - 9, &ntot, dims);

      int Npts;
      Npts = ygets_i(argc - 10);

      if (argc > 10) {
        odevice = ygets_i(argc - 11);
      }

      handle->device = odevice;

      handle->sutra_target = new sutra_target_brama(context_handle, brama_name,
                                                    -1, ntargets, xpos, ypos,
                                                    lambda, mag, zerop, sizes, pup,
                                                    Npts, odevice);
      handle->use_brama = 1;

    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with sutra_target_brama construction in "
          << __FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

  void Y_target_publish(int argc) {
    try {
      target_struct* handle = yoga_ao_getyTarget(argc, 1);
      if (handle->use_brama == 0)
        y_error("DDS is not initialized");

      sutra_target_brama *target_handler =
          (sutra_target_brama *) (handle->sutra_target);

      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDevice(handle->device,1);

      target_handler->publish();
    } catch (string &msg) {
      y_error(msg.c_str());
    } catch (char const * msg) {
      y_error(msg);
    } catch (...) {
      stringstream buf;
      buf << "unknown error with target_initDDS in " << __FILE__ << "@"
          << __LINE__ << endl;
      y_error(buf.str().c_str());
    }
  }

}
