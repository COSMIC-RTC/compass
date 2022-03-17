// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_exception.h
//! \ingroup   libcarma
//! \class     carma_exception
//! \brief     this class provides the exception to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License


#ifndef CARMA_EXCEPTION_H_
#define CARMA_EXCEPTION_H_

#include <iostream>
#include <string>
#include <sstream>

#define Carma_Error(s) CarmaException(s, __FILE__, __LINE__)

//! Carma exception throw by libcarma
/*!
 * \class CarmaException
 *
 * \brief Carma exception throw by libcarma
 *
 *  This class is used to throw a Carma exception (generally with libcarma)
 */

class CarmaException {
 private:
  std::string a_reason;  //!< a detailed description of the error
  std::string a_file;    //!< in which file this exception has been created
  unsigned int a_line;   //!< on which line this exception has been created
 public:
  /*!
   *  \brief Constructor
   *
   *  CarmaException Constructor
   *
   *  \param reason : a detailed description of the error
   *  \param file : which file this exception has been created
   *  \param line : which line this exception has been created
   */

  CarmaException(std::string reason, std::string file, unsigned int line)
      : a_reason(reason), a_file(file), a_line(line) {}

  /*!
   *  \brief Destructor
   *
   *  CarmaException Destructor
   */

  ~CarmaException() {}

  /*!
   *  \brief Format into a const char *
   *
   *  Format the Carma exception into a const char *
   *
   *  \return Formated exception
   */

  const char* show_reason() const {
    std::stringstream buf;
    buf << a_reason << " in " << a_file << "@" << a_line << std::endl;
    return buf.str().c_str();
  }

  /*!
   *  \brief Format into a string
   *
   *  Format the Carma exception into a string
   *
   *  \return Formated exception
   */

  std::string show_reason_str() const {
    std::stringstream buf;
    buf << a_reason << " in " << a_file << "@" << a_line << std::endl;
    return buf.str();
  }
};

#endif /* CARMA_EXCEPTION_H_ */
