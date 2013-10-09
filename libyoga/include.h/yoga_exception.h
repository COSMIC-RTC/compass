/*
 * yoga_exception.h
 *
 *  Created on: Mar 18, 2013
 *      Author: sevin
 */

#ifndef YOGA_EXCEPTION_H_
#define YOGA_EXCEPTION_H_


#include <string>

using namespace std;

  //! Yoga exception throw by FT_Check
  /*!
   * \class YogaException
   *
   * \brief Yoga exception throw by FT_Check
   *
   *  This class is used to throw a Yoga exception (generally with FT_Check)
   */

  class YogaException {
    private:
      string aReason;      //!< a detailed description of the error
      string aFile;        //!< in which file this exception has been created
      unsigned int aLine;  //!< on which line this exception has been created
    public:

      /*!
       *  \brief Constructor
       *
       *  YogaException Constructor
       *
       *  \param reason : a detailed description of the error
       *  \param file : which file this exception has been created
       *  \param line : which line this exception has been created
       */

      YogaException(string reason, string file, unsigned int line) :
          aReason(reason), aFile(file), aLine(line) {
      }

      /*!
       *  \brief Destructor
       *
       *  YogaException Destructor
       */

      ~YogaException() {
      }

      /*!
       *  \brief Format into a const char *
       *
       *  Format the Yoga exception into a const char *
       *
       *  \return Formated exception
       */

      const char *showReason() const {
        stringstream buf;
        buf << aReason << " in " << aFile << "@" << aLine << endl;
        return buf.str().c_str();
      }

      /*!
       *  \brief Format into a string
       *
       *  Format the Yoga exception into a string
       *
       *  \return Formated exception
       */

      string showReasonStr() const {
        stringstream buf;
        buf << aReason << " in " << aFile << "@" << aLine << endl;
        return buf.str();
      }
  };



#endif /* YOGA_EXCEPTION_H_ */
