/***********************************************/
/**
* @file fileMeanPolarMotion.h
*
* @brief Read/write MeanPolarMotion.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#ifndef __GROOPS_FILEMEANPOLARMOTION__
#define __GROOPS_FILEMEANPOLARMOTION__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_MeanPolarMotion
static const char *docstringMeanPolarMotion = R"(
The mean pole of the Earth rotation is represented by a polynomial in a time interval.

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="meanPolarMotion" version="20200123">
  <meanPolarMotion>
    <intervalCount>2</intervalCount>
    <interval>
      <timeStart>42778.0000000000000000</timeStart>
      <degree>3</degree>
      <xp>5.59741000000000e-02</xp>
      <xp>1.82430000000000e-03</xp>
      <xp>1.84130000000000e-04</xp>
      <xp>7.02400000000000e-06</xp>
      <yp>3.46346000000000e-01</yp>
      <yp>1.78960000000000e-03</yp>
      <yp>-1.07290000000000e-04</yp>
      <yp>-9.08000000000000e-07</yp>
    </interval>
    <interval>
      <timeStart>55197.0000000000000000</timeStart>
      <degree>1</degree>
      <xp>2.35130000000000e-02</xp>
      <xp>7.61410000000000e-03</xp>
      <yp>3.58891000000000e-01</yp>
      <yp>-6.28700000000000e-04</yp>
    </interval>
  </meanPolarMotion>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_MEANPOLARMOTION_TYPE = "meanPolarMotion";

/***** CLASS ***********************************/

/** @brief MeanPolarMotions.
* The mean polar motion is described as polynomials in time intervals. */
class MeanPolarMotion
{
public:
  std::vector<Time>   timeStart;
  std::vector<UInt>   degree;
  std::vector<Vector> xp, yp;

  void compute(const Time &time, Double &xBar, Double &yBar) const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const MeanPolarMotion &x);
template<> void load(InArchive  &ar, MeanPolarMotion &x);

/** @brief Write into a MeanPolarMotion file. */
void writeFileMeanPolarMotion(const FileName &fileName, const MeanPolarMotion &x);

/** @brief Read from a MeanPolarMotion file. */
void readFileMeanPolarMotion(const FileName &fileName, MeanPolarMotion &x);

/// @}

/***********************************************/

#endif
