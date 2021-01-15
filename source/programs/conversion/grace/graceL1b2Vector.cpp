/***********************************************/
/**
* @file graceL1b2Vector.cpp
*
* @brief Read GRACE L1B data.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads vector orientation data (positions of intruments in the satellite frame) from the GRACE SDS format.
The \configFile{outputfileVector}{matrix} is a $(3n\times1)$ matrix containing $(x,y,z)$ for each record.
The GRACE SDS format is described in "GRACE Level 1B Data Product User Handbook JPL D-22027"
given at \url{http://podaac.jpl.nasa.gov/grace/documentation.html}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Vector
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Vector, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Matrix)

/***********************************************/

void GraceL1b2Vector::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameIn;

    readConfig(config, "outputfileVector", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",        fileNameIn,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    UInt numberOfRecords;
    FileInGrace file(fileNameIn, numberOfRecords);

    Vector vec(3*numberOfRecords);

    for(UInt i=0; i<numberOfRecords; i++)
    {
      Int32    seconds;
      Byte     GRACE_id, qualflg;
      Double   mag, x, y, z;

      file>>seconds>>GRACE_id>>mag>>x>>y>>z>>FileInGrace::flag(qualflg);

      const Time time = mjd2time(51544.5) + seconds2time(seconds);
      logInfo<<time.dateTimeStr()<<", flag = "<<Int(qualflg)<<Log::endl;
      vec(3*i+0) = mag*x;
      vec(3*i+1) = mag*y;
      vec(3*i+2) = mag*z;
    }

    // =============================================

    if(!fileNameOut.empty())
    {
      logInfo<<"write data to <"<<fileNameOut<<">"<<Log::endl;
      writeFileMatrix(fileNameOut, vec);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
