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

//! \file      sutra_aotemplate.h
//! \ingroup   libsutra
//! \class     sutra_aotemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_AOTEMPLATE_H_
#define _SUTRA_AOTEMPLATE_H_

#include <sutra_wfs.h>

class sutra_aotemplate {
 public:
  int device;        // # device
  std::string type;  // a name for your data
  long dim;          // # of elements

  carma_obj<float> *d_data;  // the data
  carma_obj<float> *d_res;   // the result

  carma_context *current_context;  // the context in which it has been created

 public:
  sutra_aotemplate(carma_context *context, const char *type, long dim,
                   int device);
  sutra_aotemplate(const sutra_aotemplate &aotemplate);
  ~sutra_aotemplate();

  int fill_data(float *idata);
  int fill_data();
  int do_compute();
};
template <class T>
void comp_aotemplate(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif  // _SUTRA_AOTEMPLATE_H_
