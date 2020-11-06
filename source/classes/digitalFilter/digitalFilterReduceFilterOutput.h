/***********************************************/
/**
* @file DigitalFilterReduceFilterOutput.h
*
* @brief Remove filtered result from signal.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-09
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERONEMINUS__
#define __GROOPS_DIGITALFILTERONEMINUS__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterReduceFilterOutput = R"(
\subsection{ReduceFilterOutput}
Removes the filtered signal from the input, i.e. the input is passed
through a \configClass{digitalFilter}{digitalFilterType} with a frequency response of $1-H(f)$.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Remove filtered result from signal.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterReduceFilterOutput : public DigitalFilterBase
{
  DigitalFilterPtr filter2;

public:
  DigitalFilterReduceFilterOutput(Config &config);

  Matrix filter(const_MatrixSliceRef input) const;
  std::vector<std::complex<Double>> frequencyResponse(UInt length) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterReduceFilterOutput::DigitalFilterReduceFilterOutput(Config &config)
{
  try
  {
    readConfig(config, "filter", filter2, Config::MUSTSET, "", "remove filter output from input signal");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix DigitalFilterReduceFilterOutput::filter(const_MatrixSliceRef input) const
{
  try
  {
    return input - filter2->filter(input);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<std::complex<Double>> DigitalFilterReduceFilterOutput::frequencyResponse(UInt length) const
{
  auto H = filter2->frequencyResponse(length);
  for(UInt i=0; i<H.size(); i++)
    H.at(i) = 1.-H.at(i);
  return H;
}

/***********************************************/

#endif
