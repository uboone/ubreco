/**
 * \file CMTException.h
 *
 * \ingroup Clusmtool
 * 
 * \brief Class def header for exception classes in CMTException
 *
 * @author kazuhiro
 */

/** \addtogroup CMTException

    @{*/
#ifndef RECOTOOL_CMTEXCEPTION_H
#define RECOTOOL_CMTEXCEPTION_H

#include <iostream>
#include <exception>

namespace clusmtool {
  /**
     \class CMTException
     Generic (base) exception class
  */
  class CMTException : public std::exception{

  public:

    CMTException(std::string msg="") : std::exception(), _msg(msg)
    {}

    virtual ~CMTException() throw(){};
    virtual const char* what() const throw() 
    {return _msg.c_str(); }

  private:

    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group 

