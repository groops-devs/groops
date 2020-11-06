/***********************************************/
/**
* @file tideGeneratingPotential.h
*
* @brief Tide generating potential.
*
* @author Torsten Mayer-Guerr
* @date 2008-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_TIDEGENERATINGPOTENTIAL__
#define __GROOPS_TIDEGENERATINGPOTENTIAL__

#include "base/importStd.h"
#include "base/time.h"
#include "base/matrix.h"
#include "base/doodson.h"

/***** CLASS ***********************************/

/**
* @brief Single constituent of the tide generatin potential.
* @ingroup base
* The potential can be computed with
* @f[ V_{tide}(\lambda,\vartheta,r) = \left(\frac{r}{R}\right)^2 ((c\cos(\theta_f)+s\sin(\theta_f))C_{2m} + (s\cos(\theta_f)-c\sin(\theta_f))S_{2m}) @f]
* R is @a DEFAULT_R, m is @a d[0] of the Doodson number (long periodic, diurnal, semidiurnal).
* c and s with unit m^2/s^2.
* @see TideGeneratingPotential
* @relates TideGeneratingPotential */
class TideGeneratingConstituent : public Doodson
{
public:
  Double c;         //!< cos part of potential [m^2/s^2].
  Double s;         //!< sin part of potential [m^2/s^2].

  /// Default Constructor.
  TideGeneratingConstituent() : c(0), s(0) {}

  /// Constructor.
  TideGeneratingConstituent(const Doodson &doodson, Double _c=0, Double _s=0) : Doodson(doodson), c(_c), s(_s) {}

  /// Amplitude [m^2/s^2].
  Double amplitude() const {return sqrt(c*c+s*s);}

  /// Admitance factor [m^2/s^2].
  Double admit() const {return c+s;}

  /** @brief Doodson-Wartburg phase correction.
  * Possible values: 0, PI, +-PI/2.
  * see IERS conventions 2003, table 6.4. */
  Double xi() const;
};

/***** CLASS ***********************************/

/**
 @brief Tide generating potential of sun and moon.
* @ingroup base
*/
class TideGeneratingPotential : public std::vector<TideGeneratingConstituent>
{
public:
/** @brief Doodson-Wartburg phase correction.
* Possible values: 0, PI, +-PI/2.
* see IERS conventions 2003, table 6.4. */
Double xi(const Doodson &doodson) const;
};

/***********************************************/

inline Double TideGeneratingConstituent::xi() const
{
  if(d[0]==0)
    return ((admit()>0) ? PI : 0);
  else if(d[0]==1)
    return ((admit()>0) ? -PI/2 : PI/2);
  else if(d[0]==2)
    return ((admit()>0) ? 0 : PI);
  return 0;
}

/***********************************************/

inline Double TideGeneratingPotential::xi(const Doodson &doodson) const
{
  try
  {
    if((doodson.d[0] == 0) && doodson == Doodson(std::vector<Int>{0, 0, 1, 0, 0, 0})) // special case: SA
      return 0;
    if((doodson.d[0] == 1) && doodson == Doodson(std::vector<Int>{1, 1,-1, 0, 0, 0})) // special case: S1
      return PI;
    if((doodson.d[0] == 3) && doodson == Doodson(std::vector<Int>{3, 0, 0, 0, 0, 0})) // special case: M3
      return PI;
    if(doodson.d[0] > 2)  // non-linear tides
      return 0;
    for(UInt i=0; i<size(); i++)
      if(doodson == at(i))
        return at(i).xi();
    throw(Exception("constituent "+doodson.code()+" not in list"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
