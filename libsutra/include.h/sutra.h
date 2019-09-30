// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra.h
//! \defgroup  libsutra Sutra
//! \brief     Sutra is a library that provides OA tools with GPU acceleration
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_H_
#define _SUTRA_H_

#include "carma_context.h"
#include "sutra_acquisim.h"
#include "sutra_aotemplate.h"
#include "sutra_atmos.h"
#include "sutra_centroider.h"
#include "sutra_centroider_bpcog.h"
#include "sutra_centroider_cog.h"
#include "sutra_centroider_corr.h"
#include "sutra_centroider_maskedPix.h"
#include "sutra_centroider_pyr.h"
#include "sutra_centroider_tcog.h"
#include "sutra_centroider_wcog.h"
#include "sutra_controller.h"
#include "sutra_controller_cured.h"
#include "sutra_controller_generic.h"
#include "sutra_controller_geo.h"
#include "sutra_controller_ls.h"
#include "sutra_controller_mv.h"
#include "sutra_controller_utils.h"
#include "sutra_dm.h"
#include "sutra_gamora.h"
#include "sutra_groot.h"
#include "sutra_kl.h"
#include "sutra_lgs.h"
#include "sutra_phase.h"
#include "sutra_roket.h"
#include "sutra_rtc.h"
#include "sutra_rtc_brahma.h"
#include "sutra_rtc_brahmaListenerImpl.h"
#include "sutra_rtc_cacao.h"
#include "sutra_sensors.h"
#include "sutra_source.h"
#include "sutra_target.h"
#include "sutra_target_brahma.h"
#include "sutra_target_brahmaListenerImpl.h"
#include "sutra_telemetry.h"
#include "sutra_telescope.h"
#include "sutra_tscreen.h"
#include "sutra_utils.h"
#include "sutra_wfs.h"
#include "sutra_wfs_geom.h"
#include "sutra_wfs_pyr_pyrhr.h"

#endif  // _SUTRA_H_
