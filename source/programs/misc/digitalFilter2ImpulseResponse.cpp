/***********************************************/
/**
* @file digitalFilter2ImpulseResponse.cpp
*
* @brief Impulse response of a filter cascade.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-02
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Impulse response of a \configClass{digitalFilter}{digitalFilterType} cascade.
The impulse response is computed by filtering a sequence with \config{length} samples and a unit impulse at index \config{pulseLag}.

The \configFile{outputfileResponse}{matrix} is a matrix with the time stamp (zero at \config{pulseLag})
in the first column and the impulse response $h_k$ in the second column.

See also \program{DigitalFilter2FrequencyResponse}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Impulse response of a filter cascade.
* @ingroup programsGroup */
class DigitalFilter2ImpulseResponse
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(DigitalFilter2ImpulseResponse, SINGLEPROCESS, "impulse response of a filter cascade", Misc)

/***********************************************/

void DigitalFilter2ImpulseResponse::run(Config& config)
{
  try
  {
    FileName         fileNameOut;
    DigitalFilterPtr filter;
    UInt             length, lag;
    Double           sampling;

    readConfig(config, "outputfileResponse", fileNameOut, Config::MUSTSET,  "",    "columns: time [seconds], response");
    readConfig(config, "digitalFilter",      filter,      Config::MUSTSET,  "", "");
    readConfig(config, "length",             length,      Config::DEFAULT,  "201", "length of the impulse response");
    readConfig(config, "pulseLag",           lag,         Config::DEFAULT,  "100", "start of the pulse in the data series");
    readConfig(config, "sampling",           sampling,    Config::DEFAULT,  "1.0", "[seconds]");
    if(isCreateSchema(config)) return;

    logStatus<<"compute impulse response"<<Log::endl;
    Matrix A(length, 2);
    for(UInt i=0; i<length; i++)
      A(i,0) = sampling*i-sampling*lag;
    A(lag,1) = 1.;
    copy(filter->filter(A.column(1)), A.column(1));

    logStatus<<"write response to <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
