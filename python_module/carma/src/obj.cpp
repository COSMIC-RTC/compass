// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.hpp>

#include "obj.hpp"
#include "obj_complex.hpp"
#include "type_list.hpp"

void declare_carma_obj(py::module &mod) {
  apply<CarmaObjInterfacer, GenericTypeList<int32_t, uint16_t, uint32_t, float, double>>(mod);
  apply<CarmaObjComplexInterfacer, GenericTypeList<float, double>>(mod);
}
