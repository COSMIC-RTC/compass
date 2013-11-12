#include "yoga_template.h"
#include <yoga_api.h>

extern "C" {

  void
  Y_yoga_immult(int argc)
  {
    if (yarg_subroutine()) {
      yObj_struct *handle_out = (yObj_struct *) yget_obj(argc - 1, &yObj);
      yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 2, &yObj);
      
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDevice(handle_in->device);
      if (handle_in->type == Y_FLOAT) {
	caObjS *carma_out_handler = (caObjS *)(handle_out->carma_object);
	caObjS *carma_in_handler = (caObjS *)(handle_in->carma_object);
	multim((float *)carma_out_handler->getData(),(float *)carma_in_handler->getData(),carma_in_handler->getNbElem());
      } else y_error("double not implemented yet");
    } else y_error("can only be called as a subroutine");
  }

}
