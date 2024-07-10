// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sparse_obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaSparseObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma.hpp>
#include <carma_sparse_obj.hpp>

#include "sparse_obj.hpp"

using TypeListSparseObj = GenericTypeList<float, double>;

void declare_carmaWrap_sparse_obj(py::module &mod) {
  apply<CarmaSparseObjInterfacer, TypeListSparseObj>(mod);
}
