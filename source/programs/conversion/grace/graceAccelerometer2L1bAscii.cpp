/***********************************************/
/**
* @file graceAccelerometer2L1bAscii.cpp
*
* @brief Convert GROOPS accelerometer files to the GRACE SDS L1B ASCII format.
*
* @author Andreas Kvas
* @date 2020-06-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert GROOPS accelerometer files to the GRACE SDS L1B ASCII format.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert GROOPS accelerometer files to the GRACE SDS L1B ASCII format.
* Convert GROOPS accelerometer files to the GRACE SDS L1B ASCII format.
* @ingroup programsConversionGroup */
class GraceAccelerometer2L1bAscii
{
private:
  std::string headerVariables() const;
  std::string headerAttributes(UInt epochCount, const Time &t0, const Time &timeStart, const Time &timeEnd,
                               const std::vector<std::string> &globalAttributes) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceAccelerometer2L1bAscii, SINGLEPROCESS, "Convert GROOPS accelerometer files to the GRACE SDS L1B ASCII format.", Conversion)

/***********************************************/

void GraceAccelerometer2L1bAscii::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameAcc;
    std::string graceId;
    std::vector<std::string> globalAttributes;

    readConfig(config, "outputfileAscii",        fileNameOut,      Config::MUSTSET,  "", "ASCII outputfile");
    readConfig(config, "inputfileAccelerometer", fileNameAcc,      Config::MUSTSET,  "", "GROOPS acceleromter file");
    readConfig(config, "satelliteIdentifier",    graceId,          Config::MUSTSET,  "", "satellite identifier (A or B for GRACE, C or D for GRACE-FO)");
    readConfig(config, "globalAttributes",       globalAttributes, Config::OPTIONAL, "", "additional attributes as 'key: value' pairs");
    if(isCreateSchema(config)) return;

    AccelerometerArc accelerometer = InstrumentFile::read(fileNameAcc);

    const UInt epochCount = accelerometer.size();
    Time timeStart = accelerometer.at(0).time;
    Time timeEnd = accelerometer.at(epochCount-1).time;

    Time t0 = date2time(2000, 1, 1, 12);

    OutFile file(fileNameOut);
    file<<headerAttributes(epochCount, t0, timeStart, timeEnd, globalAttributes);
    file<<headerVariables()<<std::endl;
    file<<"# End of YAML header"<<std::endl;;
    for(UInt k = 0; k < epochCount; k++)
    {
      Int t = static_cast<Int>(std::round((accelerometer.at(k).time - t0).seconds()));
      Vector3d acc = accelerometer.at(k).acceleration;
      file<<t<<" "<<graceId<<" "<<acc.x()%"%17.15e"s<<" "<<acc.y()%"%17.15e"s<<" "<<acc.z()%"%17.15e"s;
      file<<" 0 0 0 "<<0.0%"%17.15e"s<<" "<<0.0%"%17.15e"s<<" "<<0.0%"%17.15e"s<<" 00000000"<<std::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GraceAccelerometer2L1bAscii::headerAttributes(UInt epochCount, const Time &t0, const Time &timeStart, const Time &timeEnd,
                                                          const std::vector<std::string> &globalAttributes) const
{
  std::stringstream ss;
  ss<<"header:"<<std::endl;
  ss<<"  dimensions:"<<std::endl;
  ss<<"    num_records: "<<epochCount<<std::endl;
  ss<<"  global_attributes:"<<std::endl;
  ss<<"    time_coverage_start: "<<timeStart%"%DT%H:%M:%SZ"s<<std::endl;
  ss<<"    time_coverage_stop: "<<timeEnd%"%DT%H:%M:%SZ"s<<std::endl;
  for(auto s : globalAttributes)
    ss<<"    "<<s<<std::endl;
  ss<<"non-standard_attributes:"<<std::endl;
  ss<<"  epoch_time: "<<t0%"%DT%H:%M:%SZ"s<<std::endl;
  ss<<"  start_time_epoch_secs: "<<(timeStart-t0).seconds()%"%i"s<<std::endl;
  ss<<"  stop_time_epoch_secs: "<<(timeEnd-t0).seconds()%"%i"s<<std::endl;

  return ss.str();
}

/***********************************************/

std::string GraceAccelerometer2L1bAscii::headerVariables() const
{
  std::string header_string =
R"(variables:
  - gps_time:
      comment: 1st column
      coverage_content_type: referenceInformation
      long_name: Continuous seconds past 01-Jan-2000 11:59:47 UTC
      units: second
  - GRACEFO_id:
      comment: 2nd column
      coverage_content_type: referenceInformation
      long_name: satellite identifier
      units: char
      valid_range: A, B, C,D
      value_meanings:
        - A = GRACE 1 (GRACE A)
        - B = GRACE 1 (GRACE B)
        - C = GRACE-FO 1 (GRACE C)
        - D = GRACE-FO 2 (GRACE D)
  - lin_accl_x:
      comment: 3rd column
      coverage_content_type: physicalMeasurement
      long_name: Linear acceleration along X-axis
      units: m/s2
  - lin_accl_y:
      comment: 4th column
      coverage_content_type: physicalMeasurement
      long_name: Linear acceleration along Y-axis
      units: m/s2
  - lin_accl_z:
      comment: 5th column
      coverage_content_type: physicalMeasurement
      long_name: Linear acceleration along Z-axis
      units: m/s2
  - ang_accl_x:
      comment: 6th column
      coverage_content_type: physicalMeasurement
      long_name: Angular acceleration about X-axis, 0 for ACT1B
      units: rad/s2
  - ang_accl_y:
      comment: 7th column
      coverage_content_type: physicalMeasurement
      long_name: Angular acceleration about Y-axis, 0 for ACT1B
      units: rad/s2
  - ang_accl_z:
      comment: 8th column
      coverage_content_type: physicalMeasurement
      long_name: Angular acceleration about Z-axis, 0 for ACT1B
      units: rad/s2
  - acl_x_res:
      comment: 9th column
      coverage_content_type: modelResult
      long_name: Residual of lin_accl_x from non-CRN-filtered value, always 0
      units: m/s2
  - acl_y_res:
      comment: 10th column
      coverage_content_type: modelResult
      long_name: Residual of lin_accl_y from non-CRN-filtered value, always 0
      units: m/s2
  - acl_z_res:
      comment: 11th column
      coverage_content_type: modelResult
      long_name: Residual of lin_accl_z from non-CRN-filtered value, always 0
      units: m/s2
  - qualflg:
      comment: 12th column
      coverage_content_type: qualityInformation
      flag_masks: 1b, 2b, 4b, 8b, 16b, 32b, 64b, 128b
      flag_meanings:
        - bit 0 = Not defined
        - bit 1 = Not defined
        - bit 2 = Not defined
        - bit 3 = Not defined
        - bit 4 = Not defined)";

    return header_string;
}

/***********************************************/
