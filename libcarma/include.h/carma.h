// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma.h
//! \defgroup  libcarma Carma
//! \brief     Carma is a library that provides GPU acceleration
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _CARMA_H_
#define _CARMA_H_

#include "carma_context.h"
#include "carma_cublas.h"
#include "carma_cusparse.h"
#include "carma_exception.h"
// #include "carma_fft.h"
#include "carma_cusolver.h"
#include "carma_host_obj.h"
#include "carma_ipcs.h"
#include "carma_magma.h"
#include "carma_multithread.h"
#include "carma_obj.h"
// #include "carma_sparse_host_obj.h"
// #include "carma_sparse_obj.h"
#include "carma_streams.h"
#include "carma_timer.h"
#include "carma_utils.h"

#endif  // _CARMA_H_
