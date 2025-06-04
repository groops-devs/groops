/***********************************************/
/**
* @file fileInstrument.cpp
*
* @brief Satellites instrument data organized in arcs.
*
* @author Torsten Mayer-Guerr
* @date 2004-11-29
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Instrument

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "inputOutput/logging.h"
#include "files/fileFormatRegister.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

GROOPS_REGISTER_FILEFORMAT(Instrument, FILE_INSTRUMENT_TYPE)

/***********************************************/

class MiscValuesOldEpoch : public MiscValuesEpoch
{
public:
  MiscValuesOldEpoch() : MiscValuesEpoch(0) {}
  virtual void load(InArchive  &ia) override;
};

/***********************************************/

Epoch *Epoch::create(Type type)
{
  try
  {
    if(type > 0)
      return new MiscValuesEpoch(type);
    switch(type)
    {
      case Epoch::INSTRUMENTTIME:        return new InstrumentTimeEpoch();
      case Epoch::MISCVALUE:             return new MiscValueEpoch();
      case Epoch::VECTOR3D:              return new Vector3dEpoch();
      case Epoch::COVARIANCE3D:          return new Covariance3dEpoch();
      case Epoch::ORBIT:                 return new OrbitEpoch();
      case Epoch::STARCAMERA:            return new StarCameraEpoch();
      case Epoch::ACCELEROMETER:         return new AccelerometerEpoch();
      case Epoch::SATELLITETRACKING:     return new SatelliteTrackingEpoch();
      case Epoch::GRADIOMETER:           return new GradiometerEpoch();
      case Epoch::GNSSRECEIVER:          return new GnssReceiverEpoch();
      case Epoch::OBSERVATIONSIGMA:      return new ObservationSigmaEpoch();
      case Epoch::MASS:                  return new MassEpoch();
      case Epoch::THRUSTER:              return new ThrusterEpoch();
      case Epoch::MAGNETOMETER:          return new MagnetometerEpoch();
      case Epoch::ACCHOUSEKEEPING:       return new AccHousekeepingEpoch();
      case Epoch::CLOCK:                 return new ClockEpoch();
      case Epoch::STARCAMERA1A:          return new StarCamera1AEpoch();
      case Epoch::ACCELEROMETER1A:       return new Accelerometer1AEpoch();
      case Epoch::SATELLITELASERRANGING: return new SatelliteLaserRangingEpoch();
      case Epoch::METEOROLOGICAL:        return new MeteorologicalEpoch();
      case Epoch::MISCVALUESOLD:         return new MiscValuesOldEpoch();
      case Epoch::MISCVALUES:            break;
      case Epoch::EMPTY:                 break;
    }

    throw(Exception("unknown instrument type ("+static_cast<Int>(type)%"%i)"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Epoch::getTypeName(Type type)
{
  try
  {
    if(type > 0)
      return "MISCVALUES("+static_cast<Double>(type)%"%i)"s;
    switch(type)
    {
      case Epoch::MISCVALUE:             return "MISCVALUE";
      case Epoch::INSTRUMENTTIME:        return "INSTRUMENTTIME";
      case Epoch::VECTOR3D:              return "VECTOR3D";
      case Epoch::COVARIANCE3D:          return "COVARIANCE3D";
      case Epoch::ORBIT:                 return "ORBIT";
      case Epoch::STARCAMERA:            return "STARCAMERA";
      case Epoch::ACCELEROMETER:         return "ACCELEROMETER";
      case Epoch::SATELLITETRACKING:     return "SATELLITETRACKING";
      case Epoch::GRADIOMETER:           return "GRADIOMETER";
      case Epoch::GNSSRECEIVER:          return "GNSSRECEIVER";
      case Epoch::OBSERVATIONSIGMA:      return "OBSERVATIONSIGMA";
      case Epoch::MASS:                  return "MASS";
      case Epoch::THRUSTER:              return "THRUSTER";
      case Epoch::MAGNETOMETER:          return "MAGNETOMETER";
      case Epoch::ACCHOUSEKEEPING:       return "ACCHOUSEKEEPING";
      case Epoch::CLOCK:                 return "CLOCK";
      case Epoch::STARCAMERA1A:          return "STARCAMERA1A";
      case Epoch::ACCELEROMETER1A:       return "ACCELEROMETER1A";
      case Epoch::SATELLITELASERRANGING: return "SATELLITELASERRANGING";
      case Epoch::METEOROLOGICAL:        return "METEOROLOGICAL";
      case Epoch::MISCVALUESOLD:         return "MISCVALUES";
      case Epoch::MISCVALUES:            return "MISCVALUES(0)";
      case Epoch::EMPTY:                 return "EMPTY";
    }
    return "unknown";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Epoch::fileFormatString(Type type)
{
  try
  {
    if(type > 0) // MISCVALUES
    {
      std::string str = "Time [MJD]             ";
      for(Int i=0; i<type; i++)
      {
        std::string str2 = "  data"+i%"%i"s;
        str += str2 + std::string(26-str2.size(), ' ');
      }
      return str;
    }

    switch(type)
    {                                       // "Time [MJD]             |                         |                         |                         |                         |                         |                         |                         |                         |                         |
      case Epoch::MISCVALUE:             return "Time [MJD]               data0                   ";
      case Epoch::INSTRUMENTTIME:        return "Time [MJD]             ";
      case Epoch::VECTOR3D:              return "Time [MJD]               data0: x                  data1: y                  data2: z                ";
      case Epoch::COVARIANCE3D:          return "Time [MJD]               data0: xx                 data1: yy                 data2: zz                 data3: xy                 data4: xz                 data5: yz               ";
      case Epoch::ORBIT:                 return "Time [MJD]               data0: pos x [m]          data1: pos y [m]          data2: pos z [m]          data3: vel x [m/s]        data4: vel y [m/s]        data5: vel z [m/s]        data6: acc x [m/s^2]      data7: acc y [m/s^2]      data8: acc z [m/s^2]    ";
      case Epoch::STARCAMERA:            return "Time [MJD]               data0: quaternion 0       data1: quaternion x       data2: quaternion y       data3: quaternion z     ";
      case Epoch::ACCELEROMETER:         return "Time [MJD]               data0: x [m/s^2]          data1: y [m/s^2]          data2: z [m/s^2]        ";
      case Epoch::SATELLITETRACKING:     return "Time [MJD]               data0: range [m]          data1: range-rate [m/s]   data2: range-acc [m/s^2]";
      case Epoch::GRADIOMETER:           return "Time [MJD]               data0: xx [1/s^2]         data1: yy [1/s^2]         data2: zz [1/s^2]         data3: xy [1/s^2]         data4: xz [1/s^2]         data5: yz [1/s^2]       ";
      case Epoch::GNSSRECEIVER:          return "";
      case Epoch::OBSERVATIONSIGMA:      return "Time [MJD]               data0: sigma            ";
      case Epoch::SATELLITELASERRANGING: return "Time [MJD]               data0: range [m]          data1: accuracy [m]       data2: redundancy         data3: window [s]         data4: wavelength [m]     data5: azmiuth [rad]      data6: elevation [rad]  ";
      case Epoch::METEOROLOGICAL:        return "Time [MJD]               data0: temperature [K]    data1: pressure [Pa]      data2: humidity [%]       data3: windSpeed [m/s]    data4: radiation [W/m^2]  data5: precip. [mm/d]   ";
      default:                           return "";
//       case Epoch::MASS:                 return "MASS";
//       case Epoch::THRUSTER:             return "THRUSTER";
//       case Epoch::MAGNETOMETER:         return "MAGNETOMETER";
//       case Epoch::ACCHOUSEKEEPING:      return "ACCHOUSEKEEPING";
//       case Epoch::CLOCK:                return "CLOCK";
//       case Epoch::STARCAMERA1A:         return "STARCAMERA1A";
//       case Epoch::ACCELEROMETER1A:      return "ACCELEROMETER1A";
//       case Epoch::MISCVALUESOLD:        return "MISCVALUES";
//       case Epoch::EMPTY:                return "EMPTY";
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt Epoch::dataCount(Type type, Bool mustDefined)
{
  try
  {
    if(type > 0)
      return static_cast<UInt>(type);
    switch(type)
    {
      case Epoch::INSTRUMENTTIME:        return 0;
      case Epoch::MISCVALUE:             return 1;
      case Epoch::VECTOR3D:              return 3;
      case Epoch::COVARIANCE3D:          return 6;
      case Epoch::ORBIT:                 return 9;
      case Epoch::STARCAMERA:            return 4;
      case Epoch::ACCELEROMETER:         return 3;
      case Epoch::SATELLITETRACKING:     return 3;
      case Epoch::GRADIOMETER:           return 6;
      case Epoch::OBSERVATIONSIGMA:      return 1;
      case Epoch::MASS:                  return 2;
      case Epoch::THRUSTER:              return 14;
      case Epoch::MAGNETOMETER:          return 13;
      case Epoch::ACCHOUSEKEEPING:       return 13;
      case Epoch::CLOCK:                 return 6;
      case Epoch::STARCAMERA1A:          return 9;
      case Epoch::ACCELEROMETER1A:       return 5;
      case Epoch::SATELLITELASERRANGING: return 7;
      case Epoch::METEOROLOGICAL:        return 6;
      case Epoch::GNSSRECEIVER:          if(mustDefined) throw(Exception("GNSSRECEIVER: Data columns not defined.")); return NULLINDEX;
      case Epoch::MISCVALUESOLD:         if(mustDefined) throw(Exception("Cannot determine the number of data columns of old file format (use FileConvert)")); return NULLINDEX;
      case Epoch::MISCVALUES:            if(mustDefined) throw(Exception("EMPTY: Cannot determine the number of data columns.")); return NULLINDEX;
      case Epoch::EMPTY:                 if(mustDefined) throw(Exception("EMPTY: Cannot determine the number of data columns.")); return NULLINDEX;
    }
    throw(Exception("unknown instrument type ("+static_cast<Int>(type)%"%i)"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Arc::Arc(const Arc &x) : epoch(x.epoch.size())
{
  try
  {
    for(UInt i=0; i<epoch.size(); i++)
      epoch.at(i) = std::unique_ptr<Epoch>(x.epoch.at(i)->clone());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc::Arc(const_MatrixSliceRef A, Epoch::Type type)
{
  try
  {
    if(type == Epoch::EMPTY)
    {
      if(A.columns() == 1) type = Epoch::INSTRUMENTTIME;
      if(A.columns() == 2) type = Epoch::MISCVALUE;
      if(A.columns() >= 3) type = static_cast<Epoch::Type>(A.columns()-1);
    }

    epoch.resize(A.rows());
    for(UInt i=0; i<A.rows(); i++)
    {
      epoch.at(i) = std::unique_ptr<Epoch>(Epoch::create(type));
      epoch.at(i)->time = mjd2time(A(i,0));
      if(A.columns() > 1)
        epoch.at(i)->setData(A.slice(i,1,1,A.columns()-1).trans());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc::Arc(const std::vector<Time> &times, const_MatrixSliceRef A, Epoch::Type type) : Arc(A, type)
{
  try
  {
    if(times.size() != A.rows())
      throw(Exception("Dimension error: times.size = "+times.size()%"%i, A("s+A.rows()%"%i x "s+A.columns()%"%i)"s));
    for(UInt i=0; i<size(); i++)
      at(i).time = times.at(i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc &Arc::operator=(const Arc &x)
{
  try
  {
    epoch.resize(x.epoch.size());
    for(UInt i=0; i<x.epoch.size(); i++)
      epoch.at(i) = std::unique_ptr<Epoch>(x.epoch.at(i)->clone());
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::sort()
{
  try
  {
    std::stable_sort(epoch.begin(), epoch.end(), [](const std::unique_ptr<Epoch> &x, const std::unique_ptr<Epoch> &y) {return x->time < y->time;});
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::removeDuplicateEpochs(Bool keepFirst, Double margin)
{
  try
  {
    if(keepFirst)
    {
      auto new_end = std::unique(epoch.begin(), epoch.end(), [&](const std::unique_ptr<Epoch> &x, const std::unique_ptr<Epoch> &y) {return std::fabs((x->time-y->time).seconds()) < margin;});
      epoch.erase(new_end, epoch.end());
    }
    else
    {
      auto new_end = std::unique(epoch.rbegin(), epoch.rend(), [&](const std::unique_ptr<Epoch> &x, const std::unique_ptr<Epoch> &y) {return std::fabs((x->time-y->time).seconds()) < margin;});
      epoch.erase(epoch.begin(), new_end.base());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc Arc::subArc(UInt start, UInt len) const
{
  try
  {
    if(len==0)
      return Arc();
    if((start+len)>size())
      throw(Exception("Length of sub-arc ("+len%"%i) starting at ("s+start%"%i) exceeds arc length ("s+size()%"%i)."s));

    Arc arc;
    arc.epoch.resize(len);
    for(UInt i=0; i<len; i++)
      arc.epoch.at(i) = std::unique_ptr<Epoch>(epoch.at(start+i)->clone());
    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::remove(UInt start, UInt len)
{
  try
  {
    if(len==0)
      return;
    if(start+len>size())
      throw(Exception("Length of sub-arc ("+len%"%i) starting at ("s+start%"%i) exceeds arc length ("s+size()%"%i)."s));
    epoch.erase(epoch.begin()+start, epoch.begin()+(start+len));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::append(const Arc &arc)
{
  try
  {
    if(&arc==this)
      throw(Exception("append same arc: not implemented"));
    if(arc.size() == 0)
      return;
    checkType(arc.getType());

    const UInt index = size();
    epoch.resize(index+arc.size());
    for(UInt i=0; i<arc.size(); i++)
      epoch.at(i+index) = std::unique_ptr<Epoch>(arc.epoch.at(i)->clone());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::synchronize(const std::vector<Time> &time, Double margin)
{
  try
  {
    if(size() == 0)
      return;
    UInt idxT=0, idxE=0, idxW = 0;
    for(;;)
    {
      if(idxE>=epoch.size())
        break;
      while((idxT<time.size())  && ((time.at(idxT)-epoch.at(idxE)->time).seconds() < -margin))
        idxT++;
      if(idxT>=time.size())
        break;
      while((idxE<epoch.size()) && ((time.at(idxT)-epoch.at(idxE)->time).seconds() > +margin))
        idxE++;
      if(idxE>=epoch.size())
        break;
      if(std::fabs((time.at(idxT)-epoch.at(idxE)->time).seconds()) > margin)
        break;
      epoch.at(idxW++) = std::move(epoch.at(idxE++));
    }
    epoch.resize(idxW);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::divide(const Time &minGap, UInt minArcLen, std::vector<Arc> &arcList) const
{
  try
  {
    if(size() == 0)
      return;

    // quick return possible?
    if((minGap == Time()) && (size() >= minArcLen))
    {
      arcList.push_back(*this);
      return;
    }

    UInt index = 0;
    while(index<size())
    {
      // find time gap
      UInt len = 1;
      while((index+len<size()) && ((epoch.at(index+len)->time-epoch.at(index+len-1)->time) < minGap))
        len++;
      if(len>=minArcLen)
        arcList.push_back(subArc(index,len));
      index += len;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::insert(UInt index, const Epoch &e)
{
  try
  {
    checkType(e.getType());
    epoch.insert(epoch.begin()+index, std::unique_ptr<Epoch>(e.clone()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::push_back(const Epoch &e)
{
  try
  {
    checkType(e.getType());
    epoch.push_back(std::unique_ptr<Epoch>(e.clone()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Time> Arc::times() const
{
  std::vector<Time> time(epoch.size());
  for(UInt i=0; i<epoch.size(); i++)
    time.at(i) = epoch.at(i)->time;
  return time;
}

/***********************************************/

Matrix Arc::matrix() const
{
  try
  {
    if(!size())
      return Matrix();
    Matrix A(size(), 1+at(0).data().rows(), Matrix::NOFILL);
    for(UInt i=0; i<size(); i++)
      A(i,0) = at(i).time.mjd();

    if(A.columns()>1)
      for(UInt i=0; i<size(); i++)
        copy(at(i).data().trans(), A.slice(i,1,1,A.columns()-1));

    if(getType() == Epoch::STARCAMERA)
      for(UInt i=1; i<size(); i++)
        if(inner(A.slice(i-1,1,1,4), A.slice(i,1,1,4))<0)
         A.slice(i,1,1,4) *= -1.;

    return A;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::checkSynchronized(const std::vector<std::reference_wrapper<const Arc>> &arcList)
{
  std::vector<Time> timesOld;
  for(const Arc &arc : arcList)
  {
    const std::vector<Time> times = arc.times();
    if(!timesOld.size())
      timesOld = times;
    if(times.size() && timesOld.size() && (times != timesOld))
      throw(Exception("instrument arc "+arc.getTypeName()+" is not synchronous with the other arcs"));
  }
}

/***********************************************/

void Arc::printStatistics(const Arc &arc)
{
  printStatistics(std::vector<Arc>(1, arc));
}

/***********************************************/

void Arc::printStatistics(const std::vector<Arc> &arcList)
{
  try
  {
    // number of epochs
    // ----------------
    UInt epochCount = 0;
    for(UInt arcNo=0; arcNo<arcList.size(); arcNo++)
      epochCount += arcList.at(arcNo).size();

    if(epochCount == 0)
    {
      logInfo<<"  arc count: "<<arcList.size()<<Log::endl;
      logInfo<<"  epochs:    "<<epochCount<<Log::endl;
      return;
    }

    // test times
    // ----------
    Bool notSorted = FALSE;
    UInt duplicateCount = 0;
    Time timeStart = date2time(9999,1,1), timeEnd, timeLast = date2time(-9999,1,1);
    for(UInt arcNo=0; arcNo<arcList.size(); arcNo++)
      for(UInt i=0; i<arcList.at(arcNo).size(); i++)
      {
        if(arcList.at(arcNo).at(i).time == timeLast)
          duplicateCount++;
        if(arcList.at(arcNo).at(i).time < timeLast)
          notSorted = TRUE;
        timeLast  = arcList.at(arcNo).at(i).time;
        timeStart = std::min(timeStart, timeLast);
        timeEnd   = std::max(timeEnd,   timeLast);
      }

    logInfo<<"  time start:      "<<timeStart.dateTimeStr()<<Log::endl;
    logInfo<<"  time end:        "<<timeEnd.dateTimeStr()<<Log::endl;
    logInfo<<"  epochs:          "<<epochCount<<Log::endl;
    if(notSorted)
    {
      logWarning<<"  epochs are not sorted!"<<Log::endl;
      return;
    }
    if(duplicateCount)
      logInfo<<"  duplicates:      "<<duplicateCount<<Log::endl;

    // median sampling
    // ---------------
    std::vector<Time> times;
    for(const Arc &arc : arcList)
    {
      auto arcTimes = arc.times();
      times.insert(times.end(), arcTimes.begin(), arcTimes.end());
    }
    const Double sampling = medianSampling(times).seconds();
    logInfo<<"  median sampling: "<<sampling<<" seconds"<<Log::endl;

    // arc statistics
    // --------------
    if(arcList.size() == 1)
    {
      UInt countGaps = 0;
      for(UInt i=1; i<arcList.at(0).size(); i++)
        if((arcList.at(0).at(i).time-arcList.at(0).at(i-1).time).seconds() > 1.5*sampling)
          countGaps++;
      logInfo<<"  gaps:            "<<countGaps<<Log::endl;
    }
    else
    {
      Double  meanLen  = 0;
      UInt    maxLen   = 0;
      UInt    minLen   = MAX_UINT;
      Time    maxTime;
      Time    minTime = seconds2time(100*365*86400.);
      Time    meanTime;

      for(UInt arcNo=0; arcNo<arcList.size(); arcNo++)
      {
        UInt size = arcList.at(arcNo).size();
        if(size==0)
          continue;
        Time time = arcList.at(arcNo).at(size-1).time - arcList.at(arcNo).at(0).time;

        maxLen   = std::max(maxLen, size);
        minLen   = std::min(minLen, size);
        meanLen += size;

        maxTime   = std::max(maxTime, time);
        minTime   = std::min(minTime, time);
        meanTime += time;
      }

      UInt minCount = 0;
      UInt maxCount = 0;
      for(UInt arcNo=0; arcNo<arcList.size(); arcNo++)
      {
        if(maxLen == arcList.at(arcNo).size()) maxCount++;
        if(minLen == arcList.at(arcNo).size()) minCount++;
      }

      meanLen  *= 1./arcList.size();
      meanTime *= 1./arcList.size();

      logInfo<<"  arc count:       "<<arcList.size()<<Log::endl;
      logInfo<<"  max. arc length: "<<maxTime.str() <<" with "<<maxLen <<" epochs\t ("<<maxCount<<" arcs)"<<Log::endl;
      logInfo<<"  min. arc length: "<<minTime.str() <<" with "<<minLen <<" epochs\t ("<<minCount<<" arcs)"<<Log::endl;
      logInfo<<"  mean arc length: "<<meanTime.str()<<" with "<<meanLen<<" epochs"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Arc::load(InArchive  &ia)
{
  try
  {
    UInt typeInt, count;
    ia>>nameValue("type",       typeInt);
    ia>>nameValue("pointCount", count);
    epoch.resize(count);
    for(UInt i=0; i<count; i++)
    {
      ia>>beginGroup("epoch");
      epoch.at(i) = std::unique_ptr<Epoch>(Epoch::create(static_cast<Epoch::Type>(typeInt)));
      epoch.at(i)->load(ia);
      ia>>endGroup("epoch");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Arc::save(OutArchive &oa) const
{
  oa<<nameValue("type",       static_cast<Int>(getType()));
  oa<<nameValue("pointCount", size());
  std::string comment = Epoch::fileFormatString(getType());
  oa.comment(comment);
  oa.comment(std::string(comment.size(), '='));
  for(UInt i=0; i<size(); i++)
  {
    oa<<beginGroup("epoch");
    epoch.at(i)->save(oa);
    oa<<endGroup("epoch");
  }
}

/***********************************************/
/***********************************************/

void InstrumentFile::open(const FileName &name)
{
  try
  {
    close();
    if(!name.empty())
    {
      file.open(name, ""/*arbitrary type*/, std::max(FILE_INSTRUMENT_VERSION, FILE_MATRIX_VERSION));
      fileName = name;
      index = 0;
      if(file.type() == FILE_INSTRUMENT_TYPE)
      {
        if(file.version() < 20200123)
        {
          UInt typeInt;
          file>>nameValue("satelliteType", typeInt);
          switch(typeInt)
          {
            case  0: type = Epoch::EMPTY;                break;
            case 23: type = Epoch::MISCVALUESOLD;        break;
            case 15: type = Epoch::INSTRUMENTTIME;       break;
            case 22: type = Epoch::MISCVALUE;            break;
            case 24: type = Epoch::VECTOR3D;             break;
            case  8: type = Epoch::COVARIANCE3D;         break;
            case  1: type = Epoch::ORBIT;                break;
            case  2: type = Epoch::STARCAMERA;           break;
            case  3: type = Epoch::ACCELEROMETER;        break;
            case  5: type = Epoch::SATELLITETRACKING;    break;
            case  4: type = Epoch::GRADIOMETER;          break;
            case 13: type = Epoch::GNSSRECEIVER;         break;
            case 20: type = Epoch::OBSERVATIONSIGMA;     break;
            case 19: type = Epoch::MASS;                 break;
            case 17: type = Epoch::THRUSTER;             break;
            case 18: type = Epoch::MAGNETOMETER;         break;
            case 21: type = Epoch::ACCHOUSEKEEPING;      break;
            case 26: type = Epoch::CLOCK;                break;
            case 25: type = Epoch::STARCAMERA1A;         break;
            case 27: type = Epoch::ACCELEROMETER1A;      break;
            default: throw(Exception("unsupported old instrument type"));
          }
        }
        else
        {
          Int typeInt;
          file>>nameValue("satelliteType", typeInt);
          type = static_cast<Epoch::Type>(typeInt);
        }
        file>>nameValue("arcCount", arcCount_);
      }
      else if(file.type().empty() || (file.type() == FILE_MATRIX_TYPE))
      {
        arcCount_ = 1;
        file>>nameValue("matrix", A);
        if(A.columns() == 0) type = Epoch::EMPTY;
        if(A.columns() == 1) type = Epoch::INSTRUMENTTIME;
        if(A.columns() == 2) type = Epoch::MISCVALUE;
        if(A.columns() >= 3) type = static_cast<Epoch::Type>(A.columns()-1);
      }
      else
        throw(Exception("file type is '"+file.type()+"' but must be '"+FILE_INSTRUMENT_TYPE+"' or '"+FILE_MATRIX_TYPE+"'"));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentFile::close()
{
  if(!fileName.empty())
    file.close();
  fileName  = FileName();
  arcCount_ = 0;
  type      = Epoch::EMPTY;
}

/***********************************************/

Arc InstrumentFile::readArc(UInt i)
{
  try
  {
    if(fileName.empty())
      return Arc();

    if(i>=arcCount_)
      throw(Exception("index >= arcCount"));

    // behind arc in file -> restart at beginning
    if(i<index)
      open(FileName(fileName));

    // special case: convert matrix to instrument arc
    if(file.type().empty() || (file.type() == FILE_MATRIX_TYPE))
    {
      Matrix B;
      std::swap(A, B);
      index++;
      return Arc(B, type);
    }

    Arc arc;
    std::unique_ptr<Epoch> epoch;
    while(index <= i)
    {
      arc = Arc();
      UInt count;
      file>>beginGroup("arc");
      file>>nameValue("pointCount", count);
      for(UInt i=0; i<count; i++)
      {
        if(!epoch)
          epoch = std::unique_ptr<Epoch>(Epoch::create(type));
        file>>nameValue("epoch", *epoch);
        arc.push_back(*epoch);
      }
      file>>endGroup("arc");
      index++;
    }
    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc InstrumentFile::read(const FileName &name)
{
  try
  {
    InstrumentFile file(name);
    Arc arc;
    for(UInt arcNo=0; arcNo<file.arcCount(); arcNo++)
      arc.append(file.readArc(arcNo));
    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentFile::checkArcCount(const std::vector<std::reference_wrapper<const InstrumentFile>> &fileList)
{
  UInt arcCountOld = 0;
  for(const InstrumentFile &file : fileList)
  {
    if(!arcCountOld)
      arcCountOld = file.arcCount();
    if(arcCountOld && file.arcCount() && (file.arcCount() != arcCountOld))
      throw(Exception("number of arcs ("+file.arcCount()%"%i) of instrument file <"s+file.fileName.str()+"> is different compared to the other arcs ("+arcCountOld%"%i)"s));
  }
}

/***********************************************/
/***********************************************/

void MiscValuesOldEpoch::load(InArchive &ia)
{
  if(ia.version() < 20170920)
  {
    ia >> nameValue ("time",   time);
    ia >> nameValue ("values", values);
    return;
  }
  UInt count;
  ia >> nameValue ("time",  time);
  ia >> nameValue ("count", count);
  values = Vector(count);
  for(UInt i=0; i<count; i++)
    ia >> nameValue ("value", values(i));
}

/***********************************************/
/***********************************************/

void MiscValuesEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time", time);
  for(UInt i=0; i<values.rows(); i++)
    oa << nameValue ("value", values(i));
}

/***********************************************/

void MiscValuesEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time", time);
  for(UInt i=0; i<values.size(); i++)
    ia >> nameValue ("value", values(i));
}

/***********************************************/

Vector MiscValuesEpoch::data() const
{
  return values;
}

/***********************************************/

void MiscValuesEpoch::setData(const Vector &x)
{
  values = x;
}

/***********************************************/
/***********************************************/

void InstrumentTimeEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time", time);
}

/***********************************************/

void InstrumentTimeEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time", time);
}

/***********************************************/

Vector InstrumentTimeEpoch::data() const
{
  return Vector(0);
}

/***********************************************/

void InstrumentTimeEpoch::setData(const Vector &/*x*/)
{
}

/***********************************************/
/***********************************************/

void MiscValueEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",  time);
  oa << nameValue ("value", value);
}

/***********************************************/

void MiscValueEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",  time);
  ia >> nameValue ("value", value);
}

/***********************************************/

Vector MiscValueEpoch::data() const
{
  Vector x(1);
  x(0) = value;
  return x;
}

/***********************************************/

void MiscValueEpoch::setData(const Vector &x)
{
  value = x(0);
}

/***********************************************/
/***********************************************/

void Vector3dEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",     time);
  oa << nameValue ("vector3d", vector3d);
}

/***********************************************/

void Vector3dEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",     time);
  ia >> nameValue ("vector3d", vector3d);
}

/***********************************************/

Vector Vector3dEpoch::data() const
{
  return vector3d.vector();
}

/***********************************************/

void Vector3dEpoch::setData(const Vector &x)
{
  vector3d.x() = x(0);
  vector3d.y() = x(1);
  vector3d.z() = x(2);
}

/***********************************************/
/***********************************************/

void Covariance3dEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",       time);
  oa << nameValue ("covariance", covariance);
}

/***********************************************/

void Covariance3dEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",       time);
  ia >> nameValue ("covariance", covariance);
}

/***********************************************/

Vector Covariance3dEpoch::data() const
{
  Vector x(6);
  x(0) = covariance.xx();
  x(1) = covariance.yy();
  x(2) = covariance.zz();
  x(3) = covariance.xy();
  x(4) = covariance.xz();
  x(5) = covariance.yz();
  return x;
}

/***********************************************/

void Covariance3dEpoch::setData(const Vector &x)
{
  covariance.xx() = x(0);
  covariance.yy() = x(1);
  covariance.zz() = x(2);
  covariance.xy() = x(3);
  covariance.xz() = x(4);
  covariance.yz() = x(5);
}

/***********************************************/
/***********************************************/

void OrbitEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",         time);
  oa << nameValue ("position",     position);
  oa << nameValue ("velocity",     velocity);
  oa << nameValue ("acceleration", acceleration);
}

/***********************************************/

void OrbitEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",         time);
  ia >> nameValue ("position",     position);
  ia >> nameValue ("velocity",     velocity);
  ia >> nameValue ("acceleration", acceleration);
}

/***********************************************/

Vector OrbitEpoch::data() const
{
  Vector x(9);
  x(0) = position.x();
  x(1) = position.y();
  x(2) = position.z();
  x(3) = velocity.x();
  x(4) = velocity.y();
  x(5) = velocity.z();
  x(6) = acceleration.x();
  x(7) = acceleration.y();
  x(8) = acceleration.z();
  return x;
}

/***********************************************/

void OrbitEpoch::setData(const Vector &x)
{
  position.x()     = x(0);
  position.y()     = x(1);
  position.z()     = x(2);
  velocity.x()     = x(3);
  velocity.y()     = x(4);
  velocity.z()     = x(5);
  acceleration.x() = x(6);
  acceleration.y() = x(7);
  acceleration.z() = x(8);
}

/***********************************************/
/***********************************************/

void StarCameraEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",   time);
  oa << nameValue ("rotary", rotary);
}

/***********************************************/

void StarCameraEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",   time);
  ia >> nameValue ("rotary", rotary);
}

/***********************************************/

Vector StarCameraEpoch::data() const
{
  return rotary.quaternion();
}

/***********************************************/

void StarCameraEpoch::setData(const Vector &x)
{
  rotary = Rotary3d(1./norm(x)*x);
}

/***********************************************/
/***********************************************/

void AccelerometerEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",         time);
  oa << nameValue ("acceleration", acceleration);
}

/***********************************************/

void AccelerometerEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",         time);
  ia >> nameValue ("acceleration", acceleration);
}

/***********************************************/

Vector AccelerometerEpoch::data() const
{
  Vector x(3);
  x(0) = acceleration.x();
  x(1) = acceleration.y();
  x(2) = acceleration.z();
  return x;
}

/***********************************************/

void AccelerometerEpoch::setData(const Vector &x)
{
  acceleration.x() = x(0);
  acceleration.y() = x(1);
  acceleration.z() = x(2);
}

/***********************************************/
/***********************************************/

void SatelliteTrackingEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",       time);
  oa << nameValue ("range",      range);
  oa << nameValue ("rangeRate",  rangeRate);
  oa << nameValue ("rangeAcceleration", rangeAcceleration);
}

/***********************************************/

void SatelliteTrackingEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",       time);
  ia >> nameValue ("range",      range);
  ia >> nameValue ("rangeRate",  rangeRate);
  ia >> nameValue ("rangeAcceleration", rangeAcceleration);
  if(ia.version() < 20200123)
  {
    UInt phaseIndex;
    ia >> nameValue ("phaseIndex", phaseIndex);
  }
}

/***********************************************/

Vector SatelliteTrackingEpoch::data() const
{
  Vector x(3);
  x(0) = range;
  x(1) = rangeRate;
  x(2) = rangeAcceleration;
  return x;
}

/***********************************************/

void SatelliteTrackingEpoch::setData(const Vector &x)
{
  range             = x(0);
  rangeRate         = x(1);
  rangeAcceleration = x(2);
}

/***********************************************/
/***********************************************/

void GradiometerEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",            time);
  oa << nameValue ("gravityGradient", gravityGradient);
}

/***********************************************/

void GradiometerEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",            time);
  ia >> nameValue ("gravityGradient", gravityGradient);
}

/***********************************************/

Vector GradiometerEpoch::data() const
{
  Vector x(6);
  x(0) = gravityGradient.xx();
  x(1) = gravityGradient.yy();
  x(2) = gravityGradient.zz();
  x(3) = gravityGradient.xy();
  x(4) = gravityGradient.xz();
  x(5) = gravityGradient.yz();
  return x;
}

/***********************************************/

void GradiometerEpoch::setData(const Vector &x)
{
  gravityGradient.xx() = x(0);
  gravityGradient.yy() = x(1);
  gravityGradient.zz() = x(2);
  gravityGradient.xy() = x(3);
  gravityGradient.xz() = x(4);
  gravityGradient.yz() = x(5);
}

/***********************************************/
/***********************************************/

void GnssReceiverEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",        time);
  oa << nameValue ("type",        obsType);
  oa << nameValue ("satellite",   satellite);
  oa << nameValue ("observation", observation);
  oa << nameValue ("clockError",  clockError);
}

/***********************************************/

void GnssReceiverEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",        time);
  ia >> nameValue ("type",        obsType);
  ia >> nameValue ("satellite",   satellite);
  ia >> nameValue ("observation", observation);
  ia >> nameValue ("clockError",  clockError);
}

/***********************************************/

Vector GnssReceiverEpoch::data() const
{
  try
  {
    throw(Exception("not implemented yet"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverEpoch::setData(const Vector& /*data*/)
{
  try
  {
    throw(Exception("not implemented yet"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void ObservationSigmaEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",  time);
  oa << nameValue ("sigma", sigma);
}

/***********************************************/

void ObservationSigmaEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",  time);
  ia >> nameValue ("sigma", sigma);
}

/***********************************************/

Vector ObservationSigmaEpoch::data() const
{
  Vector x(1);
  x(0) = sigma;
  return x;
}

/***********************************************/

void ObservationSigmaEpoch::setData(const Vector &x)
{
  sigma = x(0);
}

/***********************************************/
/***********************************************/

void MassEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",         time);
  oa << nameValue ("massThruster", massThr);
  oa << nameValue ("massTank",     massTank);
}

/***********************************************/

void MassEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",         time);
  ia >> nameValue ("massThruster", massThr);
  ia >> nameValue ("massTank",     massTank);
}

/***********************************************/

Vector MassEpoch::data() const
{
  Vector x(2);
  x(0) = massThr;
  x(1) = massTank;
  return x;
}

/***********************************************/

void MassEpoch::setData(const Vector &x)
{
  massThr  = x(0);
  massTank = x(1);
}

/***********************************************/
/***********************************************/

void ThrusterEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",     time);
  oa << nameValue ("onTime1",  onTime1);
  oa << nameValue ("onTime2",  onTime2);
  oa << nameValue ("onTime3",  onTime3);
  oa << nameValue ("onTime4",  onTime4);
  oa << nameValue ("onTime5",  onTime5);
  oa << nameValue ("onTime6",  onTime6);
  oa << nameValue ("onTime7",  onTime7);
  oa << nameValue ("onTime8",  onTime8);
  oa << nameValue ("onTime9",  onTime9);
  oa << nameValue ("onTime10", onTime10);
  oa << nameValue ("onTime11", onTime11);
  oa << nameValue ("onTime12", onTime12);
  oa << nameValue ("onTime13", onTime13);
  oa << nameValue ("onTime14", onTime14);
}

/***********************************************/

void ThrusterEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",     time);
  ia >> nameValue ("onTime1",  onTime1);
  ia >> nameValue ("onTime2",  onTime2);
  ia >> nameValue ("onTime3",  onTime3);
  ia >> nameValue ("onTime4",  onTime4);
  ia >> nameValue ("onTime5",  onTime5);
  ia >> nameValue ("onTime6",  onTime6);
  ia >> nameValue ("onTime7",  onTime7);
  ia >> nameValue ("onTime8",  onTime8);
  ia >> nameValue ("onTime9",  onTime9);
  ia >> nameValue ("onTime10", onTime10);
  ia >> nameValue ("onTime11", onTime11);
  ia >> nameValue ("onTime12", onTime12);
  ia >> nameValue ("onTime13", onTime13);
  ia >> nameValue ("onTime14", onTime14);
}

/***********************************************/

Vector ThrusterEpoch::data() const
{
  Vector x(14);
  x(0)  = onTime1;
  x(1)  = onTime2;
  x(2)  = onTime3;
  x(3)  = onTime4;
  x(4)  = onTime5;
  x(5)  = onTime6;
  x(6)  = onTime7;
  x(7)  = onTime8;
  x(8)  = onTime9;
  x(9)  = onTime10;
  x(10) = onTime11;
  x(11) = onTime12;
  x(12) = onTime13;
  x(13) = onTime14;
  return x;
}

/***********************************************/

void ThrusterEpoch::setData(const Vector &x)
{
  onTime1  = x(0);
  onTime2  = x(1);
  onTime3  = x(2);
  onTime4  = x(3);
  onTime5  = x(4);
  onTime6  = x(5);
  onTime7  = x(6);
  onTime8  = x(7);
  onTime9  = x(8);
  onTime10 = x(9);
  onTime11 = x(10);
  onTime12 = x(11);
  onTime13 = x(12);
  onTime14 = x(13);
}

/***********************************************/
/***********************************************/

void MagnetometerEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",                     time);
  oa << nameValue ("magneticField",            magneticField);
  oa << nameValue ("torquerA",                 torquerA);
  oa << nameValue ("torquerB",                 torquerB);
  oa << nameValue ("magneticFieldCalibration", magneticFieldCalibration);
  oa << nameValue ("torquerCalibration",       torquerCalibration);
}

/***********************************************/

void MagnetometerEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",                     time);
  ia >> nameValue ("magneticField",            magneticField);
  ia >> nameValue ("torquerA",                 torquerA);
  ia >> nameValue ("torquerB",                 torquerB);
  ia >> nameValue ("magneticFieldCalibration", magneticFieldCalibration);
  ia >> nameValue ("torquerCalibration",       torquerCalibration);
}

/***********************************************/

Vector MagnetometerEpoch::data() const
{
  Vector x(13);
  x(0)   = magneticField.x();
  x(1)   = magneticField.y();
  x(2)   = magneticField.z();
  x(3)   = torquerA.x();
  x(4)   = torquerA.y();
  x(5)   = torquerA.z();
  x(6)   = torquerB.x();
  x(7)   = torquerB.y();
  x(8)   = torquerB.z();
  x(9)   = magneticFieldCalibration.x();
  x(10)  = magneticFieldCalibration.y();
  x(11)  = magneticFieldCalibration.z();
  x(12)  = torquerCalibration;
  return x;
}

/***********************************************/

void MagnetometerEpoch::setData(const Vector &x)
{
  magneticField.x()            = x(0);
  magneticField.y()            = x(1);
  magneticField.z()            = x(2);
  torquerA.x()                 = x(3);
  torquerA.y()                 = x(4);
  torquerA.z()                 = x(5);
  torquerB.x()                 = x(6);
  torquerB.y()                 = x(7);
  torquerB.z()                 = x(8);
  magneticFieldCalibration.x() = x(9);
  magneticFieldCalibration.y() = x(10);
  magneticFieldCalibration.z() = x(11);
  torquerCalibration           = x(12);
}

/***********************************************/
/***********************************************/

void AccHousekeepingEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",               time);
  oa << nameValue ("biasVoltage",        biasVoltage);
  oa << nameValue ("vd",                 vd);
  oa << nameValue ("xOut",               xOut);
  oa << nameValue ("yOut",               yOut);
  oa << nameValue ("temperatureSU",      tempSU);
  oa << nameValue ("temperatureICU",     tempICU);
  oa << nameValue ("temperatureCore",    tempCore);
  oa << nameValue ("temperatureICUConv", tempICUConv);
  oa << nameValue ("blockNumber",        blkNrICU);
}

/***********************************************/

void AccHousekeepingEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",               time);
  ia >> nameValue ("biasVoltage",        biasVoltage);
  ia >> nameValue ("vd",                 vd);
  ia >> nameValue ("xOut",               xOut);
  ia >> nameValue ("yOut",               yOut);
  ia >> nameValue ("temperatureSU",      tempSU);
  ia >> nameValue ("temperatureICU",     tempICU);
  ia >> nameValue ("temperatureCore",    tempCore);
  ia >> nameValue ("temperatureICUConv", tempICUConv);
  ia >> nameValue ("blockNumber",        blkNrICU);
}

/***********************************************/

Vector AccHousekeepingEpoch::data() const
{
  Vector x(13);
  x(0)  = biasVoltage;
  x(1)  = vd;
  x(2)  = xOut.x();
  x(3)  = xOut.y();
  x(4)  = xOut.z();
  x(5)  = yOut.x();
  x(6)  = yOut.y();
  x(7)  = yOut.z();
  x(8)  = tempSU;
  x(9)  = tempICU;
  x(10) = tempCore;
  x(11) = tempICUConv;
  x(12) = blkNrICU;
  return x;
}

/***********************************************/

void AccHousekeepingEpoch::setData(const Vector &x)
{
  biasVoltage = x(0);
  vd          = x(1);
  xOut.x()    = x(2);
  xOut.y()    = x(3);
  xOut.z()    = x(4);
  yOut.x()    = x(5);
  yOut.y()    = x(6);
  yOut.z()    = x(7);
  tempSU      = x(8);
  tempICU     = x(9);
  tempCore    = x(10);
  tempICUConv = x(11);
  blkNrICU    = x(12);
}

/***********************************************/
/***********************************************/

void ClockEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",        time);
  oa << nameValue ("rcvTime",     rcvTime);
  oa << nameValue ("epsTime",     epsTime);
  oa << nameValue ("epsError",    epsError);
  oa << nameValue ("epsDrift",    epsDrift);
  oa << nameValue ("driftError",  driftError);
  oa << nameValue ("qualityFlag", qualityFlag);
}

/***********************************************/

void ClockEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",        time);
  ia >> nameValue ("rcvTime",     rcvTime);
  ia >> nameValue ("epsTime",     epsTime);
  ia >> nameValue ("epsError",    epsError);
  ia >> nameValue ("epsDrift",    epsDrift);
  ia >> nameValue ("driftError",  driftError);
  ia >> nameValue ("qualityFlag", qualityFlag);
}

/***********************************************/

Vector ClockEpoch::data() const
{
  Vector x(6);
  x(0)  = rcvTime;
  x(1)  = epsTime;
  x(2)  = epsError;
  x(3)  = epsDrift;
  x(4)  = driftError;
  x(5)  = qualityFlag;
  return x;
}

/***********************************************/

void ClockEpoch::setData(const Vector &x)
{
  rcvTime     = x(0);
  epsTime     = x(1);
  epsError    = x(2);
  epsDrift    = x(3);
  driftError  = x(4);
  qualityFlag = x(5);
}

/***********************************************/
/***********************************************/

void StarCamera1AEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",      time);
  oa << nameValue ("rcvTime",   rcvTime);
  oa << nameValue ("epsTime",   epsTime);
  oa << nameValue ("scaDesign", scaDesign);
  oa << nameValue ("q0",        q0);
  oa << nameValue ("q1",        q1);
  oa << nameValue ("q2",        q2);
  oa << nameValue ("q3",        q3);
  oa << nameValue ("nLocks",    nLocks);
  oa << nameValue ("nStars",    nStars);
}

/***********************************************/

void StarCamera1AEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",      time);
  ia >> nameValue ("rcvTime",   rcvTime);
  ia >> nameValue ("epsTime",   epsTime);
  ia >> nameValue ("scaDesgin", scaDesign);
  ia >> nameValue ("q0",        q0);
  ia >> nameValue ("q1",        q1);
  ia >> nameValue ("q2",        q2);
  ia >> nameValue ("q3",        q3);
  ia >> nameValue ("nLocks",    nLocks);
  ia >> nameValue ("nStars",    nStars);
}

/***********************************************/

Vector StarCamera1AEpoch::data() const
{
  Vector x(9);
  x(0) = rcvTime;
  x(1) = epsTime;
  x(2) = scaDesign;
  x(3) = q0;
  x(4) = q1;
  x(5) = q2;
  x(6) = q3;
  x(7) = nLocks;
  x(8) = nStars;
  return x;
}

/***********************************************/

void StarCamera1AEpoch::setData(const Vector &x)
{
  rcvTime   = x(0);
  epsTime   = x(1);
  scaDesign = x(2);
  q0        = x(3);
  q1        = x(4);
  q2        = x(5);
  q3        = x(6);
  nLocks    = x(7);
  nStars    = x(8);
}

/***********************************************/
/***********************************************/

void Accelerometer1AEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",         time);
  oa << nameValue ("rcvTimeInt",   rcvTimeInt);
  oa << nameValue ("rcvTimeFrac",  rcvTimeFrac);
  oa << nameValue ("acceleration", acceleration);
}

/***********************************************/

void Accelerometer1AEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",         time);
  ia >> nameValue ("rcvTimeInt",   rcvTimeInt);
  ia >> nameValue ("rcvTimeFrac",  rcvTimeFrac);
  ia >> nameValue ("acceleration", acceleration);
}

/***********************************************/

Vector Accelerometer1AEpoch::data() const
{
  Vector x(5);
  x(0) = rcvTimeInt;
  x(1) = rcvTimeFrac;
  x(2) = acceleration.x();
  x(3) = acceleration.y();
  x(4) = acceleration.z();
  return x;
}

/***********************************************/

void Accelerometer1AEpoch::setData(const Vector &x)
{
  rcvTimeInt       = x(0);
  rcvTimeFrac      = x(1);
  acceleration.x() = x(2);
  acceleration.y() = x(3);
  acceleration.z() = x(4);
}

/***********************************************/
/***********************************************/

void SatelliteLaserRangingEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",        time);
  oa << nameValue ("range",       range);
  oa << nameValue ("accuracy",    accuracy);
  oa << nameValue ("redundancy",  redundancy);
  oa << nameValue ("window",      window);
  oa << nameValue ("wavelength",  wavelength);
  oa << nameValue ("azmiuth",     azmiuth);
  oa << nameValue ("elevation",   elevation);
}

/***********************************************/

void SatelliteLaserRangingEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",       time);
  ia >> nameValue ("range",      range);
  ia >> nameValue ("accuracy",   accuracy);
  ia >> nameValue ("redundancy", redundancy);
  ia >> nameValue ("window",     window);
  ia >> nameValue ("wavelength", wavelength);
  ia >> nameValue ("azmiuth",    azmiuth);
  ia >> nameValue ("elevation",  elevation);
}

/***********************************************/

Vector SatelliteLaserRangingEpoch::data() const
{
  Vector x(7);
  x(0) = range;
  x(1) = accuracy;
  x(2) = redundancy;
  x(3) = window;
  x(4) = wavelength;
  x(5) = azmiuth;
  x(6) = elevation;
  return x;
}

/***********************************************/

void SatelliteLaserRangingEpoch::setData(const Vector &x)
{
  range      = x(0);
  accuracy   = x(1);
  redundancy = x(2);
  window     = x(3);
  wavelength = x(4);
  azmiuth    = x(5);
  elevation  = x(6);
}

/***********************************************/
/***********************************************/

void MeteorologicalEpoch::save(OutArchive &oa) const
{
  oa << nameValue ("time",           time);
  oa << nameValue ("temperature",    temperature);
  oa << nameValue ("pressure",       pressure);
  oa << nameValue ("humidity",       humidity);
  oa << nameValue ("windSpeed",      windSpeed);
  oa << nameValue ("solarRadiation", solarRadiation);
  oa << nameValue ("precipitation",  precipitation);
}

/***********************************************/

void MeteorologicalEpoch::load(InArchive &ia)
{
  ia >> nameValue ("time",           time);
  ia >> nameValue ("temperature",    temperature);
  ia >> nameValue ("pressure",       pressure);
  ia >> nameValue ("humidity",       humidity);
  ia >> nameValue ("windSpeed",      windSpeed);
  ia >> nameValue ("solarRadiation", solarRadiation);
  ia >> nameValue ("precipitation",  precipitation);
}

/***********************************************/

Vector MeteorologicalEpoch::data() const
{
  Vector x(6);
  x(0) = temperature;
  x(1) = pressure;
  x(2) = humidity;
  x(3) = windSpeed;
  x(4) = solarRadiation;
  x(5) = precipitation;
  return x;
}

/***********************************************/

void MeteorologicalEpoch::setData(const Vector &x)
{
  temperature    = x(0);
  pressure       = x(1);
  humidity       = x(2);
  windSpeed      = x(3);
  solarRadiation = x(4);
  precipitation  = x(5);
}

/***********************************************/
/***********************************************/
