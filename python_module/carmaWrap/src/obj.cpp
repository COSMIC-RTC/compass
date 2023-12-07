// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.h>

#include "obj.hpp"
#include "obj_complex.hpp"
#include "obj_half.hpp"
#include "type_list.hpp"

void declare_carmaWrap_obj(py::module &mod) {
  apply<CarmaObjInterfacer, GenericTypeList<int32_t, uint16_t, uint32_t, float, double>>(mod);
  apply<CarmaObjComplexInterfacer, GenericTypeList<float, double>>(mod);
#ifdef CAN_DO_HALF
  apply<CarmaObjHalfInterfacer, GenericTypeList<half>>(mod);
#endif
}
