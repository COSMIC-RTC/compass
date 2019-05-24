#include <wyrm>

#include "declare_name.hpp"

#include <carma.h>
#include "obj.hpp"
#include "obj_complex.hpp"
#include "obj_half.hpp"
#include "type_list.hpp"

void declare_carmaWrap_obj(py::module &mod) {
  apply<CarmaObjInterfacer, GenericTypeList<int, uint16_t, float, double>>(mod);
  apply<CarmaObjComplexInterfacer,
        GenericTypeList<cuFloatComplex, cuDoubleComplex>>(mod);
#ifdef CAN_DO_HALF
  apply<CarmaObjHalfInterfacer, GenericTypeList<half>>(mod);
#endif
}
