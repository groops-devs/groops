/***********************************************/
/**
* @file parametrizationAccelerationPerRevolution.h
*
* @brief Oscillation per revoultion.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-18
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATIONPERREVOLUTION__
#define __GROOPS_PARAMETRIZATIONACCELERATIONPERREVOLUTION__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationPerRevolution = R"(
\subsection{PerRevolution}\label{parametrizationAccelerationType:perRevolution}
Oscillation per revolution.
)";
#endif


/***********************************************/

#include "base/import.h"
#include "classes/timeSeries/timeSeries.h"
#include "parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Oscillation per revoultion.
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationPerRevolution : public ParametrizationAccelerationBase
{
  std::vector<Time>  times;
  UInt               idxStart, idxEnd;
  UInt               order, countAxis;
  Bool               estimateX, estimateY, estimateZ;
  Bool               perArc;

public:
  ParametrizationAccelerationPerRevolution(Config &config);

  Bool isPerArc() const override {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount() const override {return 2*order*countAxis*(idxEnd-idxStart);}
  void parameterName(std::vector<ParameterName> &name) const override;
  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationAccelerationPerRevolution::ParametrizationAccelerationPerRevolution(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeries;
    readConfig(config, "order",     order,      Config::MUSTSET,  "1", "once, twice, ...");
    readConfig(config, "estimateX", estimateX,  Config::DEFAULT,  "1", "along");
    readConfig(config, "estimateY", estimateY,  Config::DEFAULT,  "1", "cross");
    readConfig(config, "estimateZ", estimateZ,  Config::DEFAULT,  "1", "radial");
    readConfig(config, "interval",  timeSeries, Config::DEFAULT,  "",  "setup new parameters each interval");
    readConfig(config, "perArc",    perArc,     Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    countAxis = estimateX+estimateY+estimateZ;

    times = timeSeries->times();
    if(times.size()==0)
    {
      times.push_back( Time() );
      times.push_back( date2time(2500,1,1) );
    }
    idxStart = 0;
    idxEnd   = times.size()-1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationAccelerationPerRevolution::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    const UInt idxStartOld = idxStart;
    const UInt idxEndOld   = idxEnd;

    idxStart = 0;
    while((idxStart+1<times.size()) && (timeStart>=times.at(idxStart+1)))
      idxStart++;
    idxEnd = idxStart;
    while((idxEnd<times.size()-1) && (timeEnd>times.at(idxEnd)))
      idxEnd++;

    return (idxStartOld != idxStart) || (idxEndOld != idxEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationPerRevolution::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=idxStart; i<idxEnd; i++)
    {
      std::string dateStr;
      if(times.front() != Time())
        dateStr = times.at(i).dateTimeStr()+"_"+times.at(i+1).dateTimeStr();

      for(UInt k=0; k<order; k++)
      {
        if(estimateX) name.push_back(ParameterName("satellite", "perRevolution.cos("+(k+1)%"%i"s+"*u).x", dateStr));
        if(estimateY) name.push_back(ParameterName("satellite", "perRevolution.cos("+(k+1)%"%i"s+"*u).y", dateStr));
        if(estimateZ) name.push_back(ParameterName("satellite", "perRevolution.cos("+(k+1)%"%i"s+"*u).z", dateStr));
        if(estimateX) name.push_back(ParameterName("satellite", "perRevolution.sin("+(k+1)%"%i"s+"*u).x", dateStr));
        if(estimateY) name.push_back(ParameterName("satellite", "perRevolution.sin("+(k+1)%"%i"s+"*u).y", dateStr));
        if(estimateZ) name.push_back(ParameterName("satellite", "perRevolution.sin("+(k+1)%"%i"s+"*u).z", dateStr));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationPerRevolution::compute(SatelliteModelPtr /*satellite*/, const Time &time, const Vector3d &position, const Vector3d &velocity,
                                                     const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/, MatrixSliceRef A)
{
  try
  {
    if((time<times.at(idxStart))||(time>=times.at(idxEnd)))
      return;

    // findex index (interval)
    UInt idx = idxStart;
    while(time>=times.at(idx+1))
      idx++;

    // Argument of latitude of satellite
    Vector3d z = normalize(crossProduct(position, velocity));
    Vector3d x = normalize(crossProduct(Vector3d(0,0,1), z));
    Vector3d y = crossProduct(z, x);
    Double   u = atan2(inner(position,y), inner(position,x));

    const Matrix rotary = (rotEarth*rotSat).matrix();
    Matrix R(3, countAxis);
    UInt idxAxis = 0;
    if(estimateX) copy(rotary.column(0),  R.column(idxAxis++));
    if(estimateY) copy(rotary.column(1),  R.column(idxAxis++));
    if(estimateZ) copy(rotary.column(2),  R.column(idxAxis++));

    for(UInt k=0; k<order; k++)
    {
      axpy(1e-9*cos((k+1)*u), R, A.column(2*order*countAxis*(idx-idxStart)+countAxis*(2*k+0), countAxis));
      axpy(1e-9*sin((k+1)*u), R, A.column(2*order*countAxis*(idx-idxStart)+countAxis*(2*k+1), countAxis));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

#endif
