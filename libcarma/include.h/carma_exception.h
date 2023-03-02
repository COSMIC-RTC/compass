// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_exception.h
//! \ingroup   libcarma
//! \class     carma_exception
//! \brief     this class provides the exception to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
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

  std::string show_reason() const {
    std::stringstream buf;
    buf << a_reason << " in " << a_file << "@" << a_line << std::endl;
    return buf.str();
  }
};

#endif /* CARMA_EXCEPTION_H_ */
