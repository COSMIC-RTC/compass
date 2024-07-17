// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sparse_obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaSparseObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <carma.hpp>
#include <carma_sparse_obj.hpp>

#include "sparse_obj.hpp"

using TypeListSparseObj = GenericTypeList<float, double>;

void declare_carma_sparse_obj(py::module &mod) {
  apply<CarmaSparseObjInterfacer, TypeListSparseObj>(mod);
}
