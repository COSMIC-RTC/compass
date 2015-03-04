#include <sutra_rtc_brama.h>
#include <yoga_api.h>
#include <yoga_ao_api.h>

extern "C" {

void Y_yoga_rtc_brama(int argc) {
  try {
    if(argc==1){
      rtc_struct* handle = yoga_ao_getyRTC(argc, 1);
      SCAST(sutra_rtc*, handle_rtc,handle->sutra_rtc);
      delete handle_rtc;

      carma_context *context_handle = _getCurrentContext();
      handle->sutra_rtc = new sutra_rtc_brama(context_handle, "yoga_rtc_brama");
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
  rtc_struct* handle = yoga_ao_getyRTC(argc, 1);
  if(handle->use_brama==0)
    y_error("DDS is not initialized");

  sutra_rtc_brama *rtc_handler = (sutra_rtc_brama *) (handle->sutra_rtc);

  carma_context *context_handle = _getCurrentContext();
  context_handle->set_activeDevice(handle->device,1);

  rtc_handler->publish();
}

}
