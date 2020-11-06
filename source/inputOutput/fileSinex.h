/***********************************************/
/**
* @file fileSinex.h
*
* @brief SINEX file representation.
*
* @author Sebastian Strasser
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
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief SINEX file representation. */
class Sinex
{
  enum BlockType
  {
    UNKNOWN,
    FILE_REFERENCE,
    FILE_COMMENT,
    INPUT_HISTORY,
    INPUT_FILES,
    INPUT_ACKNOWLEDGEMENTS,
    NUTATION_DATA,
    PRECESSION_DATA,
    SOURCE_ID,
    SITE_ID,
    SITE_DATA,
    SITE_RECEIVER,
    SITE_ANTENNA,
    SITE_GPS_PHASE_CENTER,
    SITE_GAL_PHASE_CENTER,
    SITE_ECCENTRICITY,
    SATELLITE_ID,
    SATELLITE_PHASE_CENTER,
    BIAS_EPOCHS,
    BIAS_DESCRIPTION,
    BIAS_SOLUTION,
    SOLUTION_EPOCHS,
    SOLUTION_DISCONTINUITY,
    SOLUTION_STATISTICS,
    SOLUTION_ESTIMATE,
    SOLUTION_APRIORI,
    SOLUTION_MATRIX_ESTIMATE,
    SOLUTION_MATRIX_APRIORI,
    SOLUTION_NORMAL_EQUATION_VECTOR,
    SOLUTION_NORMAL_EQUATION_MATRIX,
    SATELLITE_IDENTIFIER,
    SATELLITE_PRN,
    SATELLITE_FREQUENCY_CHANNEL,
    SATELLITE_MASS,
    SATELLITE_COM,
    SATELLITE_ECCENTRICITY,
    SATELLITE_TX_POWER
  };

  class SinexBlock;
  typedef std::shared_ptr<SinexBlock> SinexBlockPtr;

  std::string _header;                        /// SINEX file header line
  std::map<BlockType, SinexBlockPtr> _blocks; /// SINEX blocks

  static BlockType blockType(const std::string &label);
  static std::string resize(std::string str, UInt length) { str.resize(length, ' '); return str; }
  static std::string format(Double value, UInt length=6, UInt precision=4);

public:
  class Epoch;
  class Parameter;
  class SinexText;
  class SinexSolutionEpochs;
  class SinexSolutionVector;
  class SinexSolutionMatrix;
  class SinexSolutionStatistics;
  typedef std::shared_ptr<SinexText>               SinexTextPtr;
  typedef std::shared_ptr<SinexSolutionEpochs>     SinexSolutionEpochsPtr;
  typedef std::shared_ptr<SinexSolutionVector>     SinexSolutionVectorPtr;
  typedef std::shared_ptr<SinexSolutionMatrix>     SinexSolutionMatrixPtr;
  typedef std::shared_ptr<SinexSolutionStatistics> SinexSolutionStatisticsPtr;

  Sinex() {}
  Sinex(const FileName &fileName) { readFile(fileName); }

  // read and write related
  void        readConfigHeader(Config &config);
  std::string header() const;
  void        readFile(const FileName &fileName);
  void        writeFile(const FileName &fileName) const;

  // block related
  Bool hasBlock(const std::string &label) const { return (_blocks.find(blockType(label)) != _blocks.end()); }
  template<typename T> std::shared_ptr<T> addBlock(const std::string &label);
  template<typename T> std::shared_ptr<T> getBlock(const std::string &label) const;
  void removeBlock(const std::string &label) { _blocks.erase(blockType(label)); }
  SinexTextPtr addSiteIdBlock(const std::vector<GnssStationInfo> &stationInfos);
  SinexTextPtr addSiteReceiverBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef);
  SinexTextPtr addSiteAntennaBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef);
  SinexTextPtr addSiteGpsPhaseCenterBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, const std::string &antennaModel);
  SinexTextPtr addSiteGalileoPhaseCenterBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, const std::string &antennaModel);
  SinexTextPtr addSiteEccentricityBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef);
  SinexTextPtr addSatelliteIdBlock(const std::vector<GnssStationInfo> &transmitterInfos, const Time &timeRef);
  SinexTextPtr addSatellitePhaseCenter(const std::vector<GnssStationInfo> &transmitterInfos, const Time &timeRef, const std::string &antennaModel);

  // conversion
  static std::string time2str(Time time, Bool fourDigitYear=FALSE);
  static Time        str2time(const std::string &line, std::size_t pos, Bool fourDigitYear=FALSE);
};

/***********************************************/

/** @brief SINEX vector block line/entry. */
class Sinex::Parameter
{
public:
  UInt        parameterIndex;
  std::string parameterType;
  std::string siteCode;
  std::string pointCode;
  std::string solutionId;
  Time        time;
  std::string unit;
  std::string constraintCode;
  Double      value;
  Double      sigma;

  Parameter() : parameterIndex(0), value(0), sigma(0) {}
};

/***********************************************/

/** @brief SINEX epochs block line/entry. */
class Sinex::Epoch
{
public:
  std::string siteCode;
  std::string pointCode;
  std::string solutionId;
  std::string observationCode;
  Time        timeStart;
  Time        timeEnd;
  Time        timeMean;


  Epoch() {}
  Epoch(const std::string &siteCode_, const std::string &pointCode_, const std::string &solutionId_,
        const std::string &observationCode_, const Time &timeStart_, const Time &timeEnd_) :
    siteCode(siteCode_), pointCode(pointCode_), solutionId(solutionId_), observationCode(observationCode_),
    timeStart(timeStart_), timeEnd(timeEnd_), timeMean(0.5*(timeStart_+timeEnd_)) {}
};

/***********************************************/

/** @brief SINEX base block representation. */
class Sinex::SinexBlock
{
  std::string _label;
  BlockType   _blockType;

protected:
  virtual BlockType blockType() const { return _blockType; }

public:
  SinexBlock(const std::string &label) : _label(label), _blockType(Sinex::blockType(label)) {}
  virtual ~SinexBlock() {}
  virtual std::string header()    const = 0;
  virtual std::string label()     const { return _label; }
  virtual void        readLine(const std::string &line) = 0;
  virtual Bool        writeBlock(std::ostream &file) const = 0;
};

/***********************************************/

/** @brief SINEX text block representation. */
class Sinex::SinexText : public Sinex::SinexBlock
{
  std::vector<std::string> _lines;

public:
  SinexText(const std::string &label) : SinexBlock(label) {}
  virtual std::string header() const;
  virtual void        readLine(const std::string &line) { _lines.push_back(line); }
  virtual Bool        writeBlock(std::ostream &file) const;
  virtual UInt        size()   const { return _lines.size(); }
  virtual std::vector<std::string> lines() const { return _lines; }
};

/***********************************************/

/** @brief SINEX epochs block representation. */
class Sinex::SinexSolutionEpochs : public Sinex::SinexBlock
{
  std::vector<Epoch> _epochs;

public:
  SinexSolutionEpochs(const std::string &label) : SinexBlock(label) {}
  virtual std::string header() const { return "*SITE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_"; }
  virtual void        readLine(const std::string &line);
  virtual Bool        writeBlock(std::ostream &file) const;
  virtual UInt        size()   const { return _epochs.size(); }
  virtual std::vector<Epoch> epochs() const { return _epochs; }
  virtual void        addEpoch(Epoch epoch) { _epochs.push_back(epoch); }
};

/***********************************************/

/** @brief SINEX vector block representation. */
class Sinex::SinexSolutionVector : public Sinex::SinexBlock
{
  std::vector<Parameter> _parameters;

public:
  SinexSolutionVector(const std::string &label) : SinexBlock(label) {}
  virtual std::string header() const;
  virtual void        readLine(const std::string &line);
  virtual Bool        writeBlock(std::ostream &file) const;
  virtual UInt        size()   const { return _parameters.size(); }
  virtual Vector      vector() const;
  virtual void        addParameter(Parameter parameter) { _parameters.push_back(parameter); }
  virtual std::vector<ParameterName> parameterNames() const;
  virtual std::vector<Parameter>     parameters()     const { return _parameters; }
};

/***********************************************/

/** @brief SINEX matrix block representation. */
class Sinex::SinexSolutionMatrix : public Sinex::SinexBlock
{
public:
  enum Type { CORRELATION, COVARIANCE, INFORMATION };

  Matrix _matrix;
  Type   _type;

  SinexSolutionMatrix(const std::string &label, const UInt size);
  virtual std::string label()  const;
  virtual std::string header() const { return "*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______"; }
  virtual void        readLine(const std::string &line);
  virtual Bool        writeBlock(std::ostream &file) const;
  virtual Type        type()   const { return _type; }
  virtual Matrix      matrix() const { return _matrix; }
  virtual void        setMatrix(const_MatrixSliceRef matrix, Type type = INFORMATION) { _matrix = matrix; _type = type; }
};

/***********************************************/

/** @brief SINEX solution statistics block representation. */
class Sinex::SinexSolutionStatistics : public Sinex::SinexBlock
{
  std::map<std::string, Double> _values;

public:
  SinexSolutionStatistics(const std::string &label) : SinexBlock(label) {}
  virtual std::string header() const { return "*____STATISTICAL_PARAMETER_____ _______VALUE(S)_______"; }
  virtual void        readLine(const std::string &line);
  virtual Bool        writeBlock(std::ostream &file) const;
  virtual Double      value(const std::string &name)  const;
  virtual void        addValue(const std::string &name, Double value) { _values[name] = value; }
};

/***********************************************/

template<typename T>
std::shared_ptr<T> Sinex::addBlock(const std::string &label)
{
  try
  {
    BlockType type = blockType(label);
    if(_blocks.find(type) != _blocks.end())
      throw(Exception("SINEX block already exists: " + label));

    if(     type == FILE_REFERENCE)                   _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == FILE_COMMENT)                     _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == INPUT_HISTORY)                    _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == INPUT_FILES)                      _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == INPUT_ACKNOWLEDGEMENTS)           _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == NUTATION_DATA)                    logWarning << "SINEX block not implemented yet: " << label << Log::endl;
    else if(type == PRECESSION_DATA)                  logWarning << "SINEX block not implemented yet: " << label << Log::endl;
    else if(type == SOURCE_ID)                        logWarning << "SINEX block not implemented yet: " << label << Log::endl;
    else if(type == SITE_ID)                          _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SITE_DATA)                        logWarning << "SINEX block not implemented yet: " << label << Log::endl;
    else if(type == SITE_RECEIVER)                    _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SITE_ANTENNA)                     _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SITE_GPS_PHASE_CENTER)            _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SITE_GAL_PHASE_CENTER)            _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SITE_ECCENTRICITY)                _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_ID)                     _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_PHASE_CENTER)           _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == BIAS_EPOCHS)                      logWarning << "SINEX block not implemented yet: " << label << Log::endl;
    else if(type == BIAS_DESCRIPTION)                 _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == BIAS_SOLUTION)                    _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SOLUTION_EPOCHS)                  _blocks[type] = SinexBlockPtr(new SinexSolutionEpochs(label));
    else if(type == SOLUTION_DISCONTINUITY)           _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SOLUTION_STATISTICS)              _blocks[type] = SinexBlockPtr(new SinexSolutionStatistics(label));
    else if(type == SOLUTION_ESTIMATE)                _blocks[type] = SinexBlockPtr(new SinexSolutionVector(label));
    else if(type == SOLUTION_APRIORI)                 _blocks[type] = SinexBlockPtr(new SinexSolutionVector(label));
    else if(type == SOLUTION_MATRIX_ESTIMATE)         _blocks[type] = SinexBlockPtr(new SinexSolutionMatrix(label, getBlock<SinexSolutionVector>("SOLUTION/ESTIMATE")->size()));
    else if(type == SOLUTION_MATRIX_APRIORI)          _blocks[type] = SinexBlockPtr(new SinexSolutionMatrix(label, getBlock<SinexSolutionVector>("SOLUTION/APRIORI")->size()));
    else if(type == SOLUTION_NORMAL_EQUATION_VECTOR)  _blocks[type] = SinexBlockPtr(new SinexSolutionVector(label));
    else if(type == SOLUTION_NORMAL_EQUATION_MATRIX)
    {
      UInt size = 0;
      try
      {
        size = getBlock<SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR")->size();
      }
      catch(std::exception &/*e*/)
      {
        size = getBlock<SinexSolutionVector>("SOLUTION/ESTIMATE")->size();
      }
      _blocks[type] = SinexBlockPtr(new SinexSolutionMatrix(label, size));
    }
    else if(type == SATELLITE_IDENTIFIER)             _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_PRN)                    _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_FREQUENCY_CHANNEL)      _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_MASS)                   _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_COM)                    _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_ECCENTRICITY)           _blocks[type] = SinexBlockPtr(new SinexText(label));
    else if(type == SATELLITE_TX_POWER)               _blocks[type] = SinexBlockPtr(new SinexText(label));
    else
      logWarning << "Cannot add unknown SINEX block: " << label << Log::endl;

    if(_blocks.find(type) != _blocks.end())
      return std::dynamic_pointer_cast<T>(_blocks.at(type));
    else
      return nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename T>
std::shared_ptr<T> Sinex::getBlock(const std::string &label) const
{
  try
  {
    auto iter = _blocks.find(blockType(label));
    if(iter == _blocks.end())
      throw(Exception("SINEX block not found: "+label));

    std::shared_ptr<T> t = std::dynamic_pointer_cast<T>(iter->second);
    if(t == nullptr)
      throw(Exception("Cannot cast to requested block type: "+label));
    return t;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
