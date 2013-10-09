#include "yoga_template.h"
#include <yoga_api.h>

extern "C" {

  void
  Y_yoga_immult(int argc)
  {
    if (yarg_subroutine()) {
      yoga_struct *handle_out = yoga_getobj(argc,1);
      yoga_struct *handle_in  = yoga_getobj(argc,2);
      
      yoga_context *context_handle = _getCurrentContext();
      context_handle->set_activeDevice(handle_in->device);
      if (handle_in->type == Y_FLOAT) {
	yObjS *yoga_out_handler = (yObjS *)(handle_out->yoga_object);
	yObjS *yoga_in_handler = (yObjS *)(handle_in->yoga_object);
	multim((float *)yoga_out_handler->getData(),(float *)yoga_in_handler->getData(),yoga_in_handler->getNbElem());
      } else y_error("double not implemented yet");
    } else y_error("can only be called as a subroutine");
  }

}
