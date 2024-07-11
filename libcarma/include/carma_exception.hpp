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

//! \file      carma_exception.hpp
//! \ingroup   libcarma
//! \class     carma_exception
//! \brief     this class provides the exception to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24


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
  uint32_t a_line;   //!< on which line this exception has been created
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

  CarmaException(std::string reason, std::string file, uint32_t line)
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

  std::string show_reason() const {
    std::stringstream buf;
    buf << a_reason << " in " << a_file << "@" << a_line << std::endl;
    return buf.str();
  }
};

#endif /* CARMA_EXCEPTION_H_ */
