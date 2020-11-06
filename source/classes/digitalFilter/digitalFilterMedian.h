/***********************************************/
/**
* @file digitalFilterMedian.h
*
* @brief Moving median filter implementation.
*
* @author Andreas Kvas
* @date 2017-02-05
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERMEDIAN__
#define __GROOPS_DIGITALFILTERMEDIAN__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterMedian = R"(
\subsection{Median}
Moving median filter of length $n$. The filter output at epoch $k$ is the median of the set start at $k-n/2$ to $k+n/2$.
The filter length $n$ should be uneven to avoid a phase shift.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Moving median filter implementation.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterMedian : public DigitalFilterBase
{
  UInt    length;
  PadType padType;

public:
  DigitalFilterMedian(Config &config);

  Matrix filter(const_MatrixSliceRef input) const;
  std::vector< std::complex<Double> > frequencyResponse(UInt length) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterMedian::DigitalFilterMedian(Config &config)
{
  try
  {
    readConfig(config, "length",  length,  Config::MUSTSET, "", "length of the moving window [epochs]");
    readConfig(config, "padType", padType, Config::MUSTSET, "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix DigitalFilterMedian::filter(const_MatrixSliceRef input) const
{
  try
  {
    Matrix padded = pad(input, length/2, padType);
    Matrix output(padded.rows(), padded.columns());

    for(UInt k=0; k<output.columns(); k++)
      for(UInt i=0; i<output.rows(); i++)
      {
        const UInt start = std::max(i,length/2)-length/2;
        const UInt end   = std::min(padded.rows(), i+(length+1)/2);
        output(i,k) = median(padded.slice(start, k, end-start, 1));
      }

    return trim(output, length/2, padType);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector< std::complex<Double> > DigitalFilterMedian::frequencyResponse(UInt /*length*/) const
{
  throw(Exception("Median filter cannot be represented by frequency response."));
}

/***********************************************/

#endif
