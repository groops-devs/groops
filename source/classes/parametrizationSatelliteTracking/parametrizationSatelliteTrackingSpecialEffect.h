/***********************************************/
/**
* @file parametrizationSatelliteTrackingSpecialEffect.h
*
* @brief Model GRACE special effects as time variable polynomials.
*
* @author Torsten Mayer-Guerr
* @author Saniya Behzadpour
* @date 2018-06-29
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKINGSPECIALEFFECT__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKINGSPECIALEFFECT__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingSpecialEffect = R"(
\subsection{SpecialEffect}\label{parametrizationSatelliteTrackingType:specialEffect}
Estimate deterministic signals in the GRACE K-Band measurements caused by Sun intrusions
into the star camera baffles of GRACE-A and eclipse transits of the satellites.
These events can be time-indexed beforehand using satellite position and orientation,
see \program{GraceSstSpecialEvents}. The shape of this short-period waveform is nearly
constant within one month and can be approximated by a polynomial.
The amplitude variation of the waveform can also be taken into account
by \configClass{parametrizationTemporal}{parametrizationTemporalType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationSatelliteTracking.h"
#include "files/fileInstrument.h"
#include "base/legendrePolynomial.h"

/***** CLASS ***********************************/

/** @brief Model GRACE special effects as time variable polynomials.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingSpecialEffect : public ParametrizationSatelliteTrackingBase
{
  Bool                       use;
  MiscValueArc               event;
  ParametrizationTemporalPtr temporal;
  UInt                       type;
  UInt                       degree;
  Double                     marginBefore, marginAfter;
  UInt                       idx;

public:
  ParametrizationSatelliteTrackingSpecialEffect(Config &config);

  Bool isPerArc() const {return FALSE;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount() const {return use * (degree+1) * temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/

inline ParametrizationSatelliteTrackingSpecialEffect::ParametrizationSatelliteTrackingSpecialEffect(Config &config)
{
  try
  {
    FileName    fileName;
    UInt        minCount;
    std::string choice;

    readConfig(config, "inputfileEvents", fileName, Config::MUSTSET, "", "instrument with GRACE events");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "eclipse1",       choice, "")) type = 1;
      if(readConfigChoiceElement(config, "eclipse2",       choice, "")) type = 2;
      if(readConfigChoiceElement(config, "starCameraBox1", choice, "")) type = 3;
      if(readConfigChoiceElement(config, "starCameraBox2", choice, "")) type = 4;
      if(readConfigChoiceElement(config, "starCameraBox3", choice, "")) type = 5;
      if(readConfigChoiceElement(config, "starCameraBox4", choice, "")) type = 6;
      if(readConfigChoiceElement(config, "starCameraBox5", choice, "")) type = 7;
      if(readConfigChoiceElement(config, "starCameraBox6", choice, "")) type = 8;
      endChoice(config);
    }
    readConfig(config, "marginLeft",           marginBefore,   Config::MUSTSET,  "20", "margin size (on both sides) [seconds]");
    readConfig(config, "marginRight",           marginAfter,    Config::MUSTSET,  "20", "margin size (on both sides) [seconds]");
    readConfig(config, "minNumberOfEvents", minCount, Config::DEFAULT,  "1",  "min. number of events to setup parameters");
    readConfig(config, "polynomialDegree",  degree,   Config::MUSTSET,  "10", "polynomial degree");
    readConfig(config, "temporal",          temporal, Config::MUSTSET,  "",   "");
    if(isCreateSchema(config)) return;

    event = InstrumentFile::read(fileName);

    UInt count = 0;
    for(UInt i=0; i<event.size(); i++)
      if(event.at(i).value == type)
        count++;
    use = (count >= minCount);

    idx = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationSatelliteTrackingSpecialEffect::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    return temporal->setInterval(timeStart, timeEnd, FALSE/*perArc*/);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingSpecialEffect::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    if(!use)
      return;

    std::string effectName;
    switch(type)
    {
      case 1: effectName = "eclipse1";       break;
      case 2: effectName = "eclipse2";       break;
      case 3: effectName = "starCameraBox1"; break;
      case 4: effectName = "starCameraBox2"; break;
      case 5: effectName = "starCameraBox3"; break;
      case 6: effectName = "starCameraBox4"; break;
      case 7: effectName = "starCameraBox5"; break;
      case 8: effectName = "starCameraBox6"; break;
    };

    std::vector<ParameterName> baseName;
    for(UInt n=0; n<=degree; n++)
      baseName.push_back(ParameterName("satellite1.satellite2", effectName+".legendrePolynomial.n"+n%"%i"s));
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingSpecialEffect::compute(UInt /*sstType*/, const std::vector<Time> &time, const Vector &/*sst0*/,
                                                          const Vector &/*position1*/, const Vector &/*position2*/, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                          const std::vector<Rotary3d> &/*rotSat1*/, const std::vector<Rotary3d> &/*rotSat2*/, MatrixSliceRef A)
{
  try
  {
    if(!use)
      return;

    UInt idx = 0;
    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      // find event
      while((idx < event.size()) && ((event.at(idx).value != type) || (event.at(idx).time+seconds2time(marginAfter) < time.at(idEpoch))))
        idx++;
      if(idx >= event.size())
        break;
      if(time.at(idEpoch) < event.at(idx).time-seconds2time(marginBefore))
        continue;

      std::vector<UInt>   index;
      std::vector<Double> factor;
      temporal->factors(event.at(idx).time, index, factor);

      const Double t  = 2*(time.at(idEpoch)-event.at(idx).time).seconds()/(marginBefore+marginAfter); // t = [-1,1]
      const Vector Pn = LegendrePolynomial::compute(t, degree);
      for(UInt i=0; i<index.size(); i++)
        axpy(factor.at(i), Pn.trans(), A.slice(idEpoch, (degree+1)*index.at(i), 1, degree+1));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
