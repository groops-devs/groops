/***********************************************/
/**
* @file fileSinex.h
*
* @brief SINEX file representation.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-05-15
*
*/
/***********************************************/

#ifndef __GROOPS_FILESINEX__
#define __GROOPS_FILESINEX__

#include "base/import.h"
#include "base/parameterName.h"
#include "inputOutput/fileName.h"
#include "inputOutput/logging.h"
#include "config/config.h"

/***** TYPES ***********************************/

class SinexBlock;
typedef std::shared_ptr<SinexBlock> SinexBlockPtr;

/***** CLASS ***********************************/

/** @brief SINEX file representation. */
class Sinex
{
public:
  std::string              header; /// SINEX file header line
  std::list<SinexBlockPtr> blocks; /// SINEX blocks

  SinexBlockPtr findBlock(const std::string &label);
  SinexBlockPtr addBlock(const std::string &label);

  // conversion
  static std::string resize(std::string str, UInt length) {str.resize(length, ' '); return str;}
  static std::string format(Double value, UInt length=6, UInt precision=4);
  static std::string time2str(Time time, Bool fourDigitYear=FALSE);
  static Time        str2time(const std::string &line, UInt pos, Bool fourDigitYear=FALSE);
};

/***********************************************/

/** @brief SINEX base block representation. */
class SinexBlock
{
public:
  std::string              label;
  std::vector<std::string> lines;
  std::stringstream        ss;

  std::ostream &operator<<(std::ostream  &(*pf)(std::ostream  &)) {pf(ss); return ss;} // sinex<<std::endl;
  template<typename T> std::ostream &operator<<(const T &t) {return ss<<t;}
};

/***** FUNCTIONS *******************************/

void writeFileSinex(const FileName &fileName, const Sinex &sinex);
void readFileSinex(const FileName &fileName, Sinex &sinex);

template<> Bool readConfig(Config &config, const std::string &name, Sinex &sinex, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

#endif /* __GROOPS__ */
