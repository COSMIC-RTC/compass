/*
 * carma_exception.h
 *
 *  Created on: Mar 18, 2013
 *      Author: sevin
 */

#ifndef CARMA_EXCEPTION_H_
#define CARMA_EXCEPTION_H_

#include <string>

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
  std::string aReason;  //!< a detailed description of the error
  std::string aFile;    //!< in which file this exception has been created
  unsigned int aLine;   //!< on which line this exception has been created
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
      : aReason(reason), aFile(file), aLine(line) {}

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

  const char* showReason() const {
    std::stringstream buf;
    buf << aReason << " in " << aFile << "@" << aLine << std::endl;
    return buf.str().c_str();
  }

  /*!
   *  \brief Format into a string
   *
   *  Format the Carma exception into a string
   *
   *  \return Formated exception
   */

  std::string showReasonStr() const {
    std::stringstream buf;
    buf << aReason << " in " << aFile << "@" << aLine << std::endl;
    return buf.str();
  }
};

#endif /* CARMA_EXCEPTION_H_ */
