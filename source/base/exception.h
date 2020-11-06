/***********************************************/
/**
* @file exception.h
*
* @brief Exception handling.
*
* @author Torsten Mayer-Guerr
* @date 2001-06-16
*
*/
/***********************************************/

#ifndef __GROOPS_EXCEPTION__
#define __GROOPS_EXCEPTION__

#include <exception>
#include <string>

/***** CLASS ***********************************/

/** @brief Exception handling.
* @ingroup base */
class Exception : public std::exception
{
  std::string message;

public:
  /** @brief Exception handling.
  * @param msg Error message */
  Exception(const std::string &msg) noexcept : message(msg) {}

  /** @brief Error message.
  * Returns a C-style character string
  * describing the general cause of the current error. */
  const char *what() const noexcept {return message.c_str();}
};

/***********************************************/

#define _GROOPS_QUOTEME_(x) #x
#define _GROOPS_QUOTEME(x) _GROOPS_QUOTEME_(x)
#define _GROOPS_ERRORLINE (std::string("in ")+__FILE__+":"+_GROOPS_QUOTEME(__LINE__)+" ("+__func__+")")

/** @brief Rethrow an error message.
* Prepend the calling function and source file.
* Is to be used in the form:
* @code
* void func()
* {
*   try
*   {
*   }
*   catch(std::exception &e)
*   {
*     GROOPS_RETHROW(e)
*   }
* }
* @endcode
* @ingroup base */
#define GROOPS_RETHROW(e) throw(Exception(_GROOPS_ERRORLINE+"\n"+e.what()));

/** @brief Rethrow an error message with extra function text.
* @see GROOPS_RETHROW
* @ingroup base */
#define GROOPS_RETHROW_EXTRA(text, e) throw(Exception(_GROOPS_ERRORLINE+": "+text+"\n"+e.what()));

/***********************************************/

#endif /* __GROOPS_Exception__ */
