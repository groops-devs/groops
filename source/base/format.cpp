/***********************************************/
/**
* @file format.cpp
*
* @brief Format of numbers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2016-07-15
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/time.h"
#include "base/format.h"
#include <regex>

/***** FUNCTIONS *******************************/

std::string operator%(LongDouble value, const std::string &format)
{
  try
  {
    UInt   year=0, month=0, day=0, hour=0, min=0;
    Double sec=0;
    Int    secondPrecision = -1;

    // parse format string -> calculate new result string
    std::string result;
    std::string::size_type posFormat = 0;
    for(;;)
    {
      if(posFormat >= format.size())
        break;
      std::string::size_type posFormatOld = posFormat;
      posFormat = format.find('%', posFormatOld);
      result += format.substr(posFormatOld, posFormat-posFormatOld);
      if(posFormat == std::string::npos)
        break;
      posFormat++;
      if(posFormat>=format.size())
        throw(Exception("expecting qualifier after '%'"));

      // parse %[flags][width][.precision]specifier format
      std::stringstream ss;
      auto c = format.at(posFormat++);

      Bool hasSubSpecifiers = FALSE;
      // [flags]
      for(;;)
      {
        if(c=='0')      ss.fill('0');
        else if(c==' ') ss.fill(' ');
        else if(c=='-') ss.setf(std::ios_base::left, std::ios_base::adjustfield);
        else if(c=='#') throw(Exception("'#' not supported"));
        else if(c=='+') throw(Exception("'+' not supported"));
        else break;
        c = format.at(posFormat++);
        hasSubSpecifiers = TRUE;
      }

      // [width]
      if(std::isdigit(c))
      {
        const char *ptr1 = format.c_str()+posFormat-1;
        char *ptr2;
        ss.width(std::strtol(ptr1, &ptr2, 10));
        posFormat += ptr2-ptr1-1;
        c = format.at(posFormat++);
        hasSubSpecifiers = TRUE;
      }

      // [.precision]
      if(c=='.')
      {
        const char *ptr1 = format.c_str()+posFormat;
        char *ptr2;
        ss.precision(std::strtol(ptr1, &ptr2, 10));
        posFormat += ptr2-ptr1;
        c = format.at(posFormat++);
        hasSubSpecifiers = TRUE;
      }

      // specifier
      switch(c)
      {
        case '%':
        {
          ss<<'%';
          break;
        }
        case 'c': // character
        {
          ss<<static_cast<char>(std::round(value));
          break;
        }
        case 'i': // integer
        {
          ss<<static_cast<int>(std::round(value));
          break;
        }
        case 'f':  // float
        case 'e':
        case 'g':
        {
          if(c=='e') ss.setf(std::ios::scientific, std::ios::floatfield);
          if(c=='f') ss.setf(std::ios::fixed,      std::ios::floatfield);
          ss<<value;
          break;
        }
        case 'O': // day of year
        {
          ss<<mjd2time(value).dayOfYear();
          break;
        }
        case 'W':
        {
          ss<<static_cast<int>(std::floor(value)-44244)/7; //GPS week
          break;
        }
        case 'w':
        {
          ss<<static_cast<int>(std::floor(value)-44244)%7; // GPS day of week
          break;
        }

        // date/time preparation
        case 'y': case 'm': case 'd':
        case 'H': case 'M': case 'S':
        case 'D': case 'T': case 'Y':
        {
          // check (once) if complete format string contains %S with precision (e.g. %.2S), find maximum precision of any %S, and round time accordingly
          if(secondPrecision == -1)
          {
            secondPrecision = 0;
            static const std::regex regexPattern("%[[\\d #+\\-]*\\.(\\d+)[S]");
            std::sregex_iterator next(format.begin(), format.end(), regexPattern);
            std::sregex_iterator end;
            while(next != end)
            {
              std::smatch match = *next;
              secondPrecision = std::max(std::atoi(match[1].str().c_str()), secondPrecision);
              next++;
            }

            // round seconds to precision considering rollovers
            mjd2time(value).date(year, month, day, hour, min, sec);
            sec = std::round(sec*std::pow(10., secondPrecision))/std::pow(10., secondPrecision) + std::pow(10., -(secondPrecision+1));
            date2time(year, month, day, hour, min, sec).date(year, month, day, hour, min, sec);
          }

          if(c=='y')      {ss<<std::setw(4)<<std::setfill('0')<<year;}
          else if(c=='m') {ss<<std::setw(2)<<std::setfill('0')<<month;}
          else if(c=='d') {ss<<std::setw(2)<<std::setfill('0')<<day;}
          else if(c=='H') {ss<<std::setw(2)<<std::setfill('0')<<hour;}
          else if(c=='M') {ss<<std::setw(2)<<std::setfill('0')<<min;}
          else if(c=='S')
          {
            if(hasSubSpecifiers)
            {
              ss.setf(std::ios::fixed, std::ios::floatfield);
              ss<<sec;
            }
            else
              ss<<std::setw(2)<<std::setfill('0')<<std::round(sec);
          }
          else if(c=='D') {ss<<std::setw(4)<<std::setfill('0')<<year<<"-"<<std::setw(2)<<std::setfill('0')<<month<<"-"<<std::setw(2)<<std::setfill('0')<<day;}
          else if(c=='T') {ss<<std::setw(2)<<std::setfill('0')<<hour<<"-"<<std::setw(2)<<std::setfill('0')<<min<<"-"<<std::setw(2)<<std::setfill('0')<<std::round(sec);}
          else if(c=='Y') {ss<<std::setw(2)<<std::setfill('0')<<year%100;}; // two digit year
          break;
        }

        default:
          throw(Exception(std::string("unknown formater '%")+c+"'"));
      } // switch
      result += ss.str();
    } // for(;;)

    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/

std::string operator%(const Time &x, const std::string &format)
{
  LongDouble value = x.mjdMod();
  value += x.mjdInt();
  return value%format;
}

/*************************************************/
