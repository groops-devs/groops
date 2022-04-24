/***********************************************/
/**
* @file fileConvert.cpp
*
* @brief Converts GROOPS file between different file formats (ASCII, XML, binary).
*
* @author Sebastian Strasser
* @date 2016-02-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts GROOPS file between different file formats (ASCII, XML, binary),
see \reference{file formats}{general.fileFormat} for details.
It prints also some information about the content.
Therefore it can be used to get an idea about the content of binary files.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileArchive.h"
#include "files/fileAdmittance.h"
#include "files/fileArcList.h"
#include "files/fileDoodsonEarthOrientationParameter.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileEarthOrientationParameter.h"
#include "files/fileEarthTide.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGriddedData.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "files/fileMatrix.h"
#include "files/fileMeanPolarMotion.h"
#include "files/fileNormalEquation.h"
#include "files/fileOceanPoleTide.h"
#include "files/fileParameterName.h"
#include "files/filePolygon.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileStringTable.h"
#include "files/fileTideGeneratingPotential.h"
#include "files/fileTimeSplinesGravityfield.h"
#include "files/fileVariationalEquation.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Converts GROOPS file between different file formats (ASCII, XML, binary).
* @ingroup programsGroup */
class FileConvert
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(FileConvert, SINGLEPROCESS, "Converts GROOPS file between different file formats (ASCII, XML, binary).",
                        System, Instrument, VariationalEquation, DoodsonHarmonics, Grid, Matrix, PotentialCoefficients, TimeSplines)

/***********************************************/

void FileConvert::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameInput, fileNameOutput;

    renameDeprecatedConfig(config, "outputFile", "outputfile", date2time(2020, 8, 20));
    renameDeprecatedConfig(config, "inputFile",  "inputfile",  date2time(2020, 8, 20));

    readConfig(config, "outputfile", fileNameOutput, Config::OPTIONAL, "", "GROOPS formats: .xml, .txt, .dat");
    readConfig(config, "inputfile",  fileNameInput,  Config::MUSTSET,  "", "GROOPS formats: .xml, .txt, .dat");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"Converting from <"<<fileNameInput<<"> to <"<< fileNameOutput<<">"<<Log::endl;
    std::string type;
    {
      InFileArchive infile(fileNameInput, ""/*arbitrary type*/);
      logInfo<<"  type    = "<<infile.type()<<Log::endl;
      logInfo<<"  version = "<<infile.version()<<Log::endl;
      type = infile.type();
    }

    if(type.empty() && !fileNameOutput.empty())
      throw(Exception("Not a GROOPS file or with very old format without type information"));

    // =============================================

    if(type == FILE_INSTRUMENT_TYPE)
    {
      InstrumentFile file(fileNameInput);
      const Epoch::Type type = file.getType();
      logInfo<<"  instr   = "<<Epoch::getTypeName(type)<<" ("<<static_cast<Int>(type)<<")"<<Log::endl;
      std::vector<Arc> arcList(file.arcCount(), type);
      for(UInt arcNo=0; arcNo<file.arcCount(); arcNo++)
        arcList.at(arcNo) = file.readArc(arcNo);
      Arc::printStatistics(arcList);
      if(!fileNameOutput.empty())
        InstrumentFile::write(fileNameOutput, arcList);
      return;
    }

    // =============================================

    if(type == FILE_VARIATIONALEQUATION_TYPE)
    {
      FileVariationalEquation file(fileNameInput);
      if(file.satellite())
        logInfo<<"  name    = "<<file.satellite()->satelliteName<<Log::endl;
      logInfo<<"  arcs    = "<<file.arcCount()<<Log::endl;
      if(!fileNameOutput.empty())
      {
        std::vector<VariationalEquationArc> arcs(file.arcCount());
        for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
          arcs.at(arcNo) = file.readArc(arcNo);
        writeFileVariationalEquation(fileNameOutput, file.satellite(), arcs);
      }
      return;
    }

    // =============================================

    if(type == FILE_DOODSONHARMONIC_TYPE)
    {
      DoodsonHarmonic x;
      readFileDoodsonHarmonic(fileNameInput,  x);
      for(UInt i=0; i<x.doodson.size(); i++)
        logInfo<<"  "<<x.doodson.at(i).code()<<" ("<<x.doodson.at(i).name()<<")"<<Log::endl;
      if(!fileNameOutput.empty())
        writeFileDoodsonHarmonic(fileNameOutput, x);
      return;
    }

    // =============================================

    if(type == FILE_MATRIX_TYPE)
    {
      Matrix x;
      readFileMatrix(fileNameInput,  x);
      std::string type = "matrix";
      if(x.getType() == Matrix::SYMMETRIC)
        type += "Symmetric"s + (x.isUpper() ? "Upper"s : "Lower"s);
      else if(x.getType() == Matrix::TRIANGULAR)
        type += "Triangular"s + (x.isUpper() ? "Upper"s : "Lower"s);
      logInfo<<"  "<<type<<"("<<x.rows()<<" x "<<x.columns()<<")"<<Log::endl;
      if(!fileNameOutput.empty())
        writeFileMatrix(fileNameOutput, x);
      return;
    }

    // =============================================

    if(type == FILE_GRIDDEDDATA_TYPE || type == "gridRectangular")
    {
      GriddedData x;
      readFileGriddedData(fileNameInput,  x);
      MiscGriddedData::printStatistics(x);
      if(!fileNameOutput.empty())
        writeFileGriddedData(fileNameOutput, x);
      return;
    }

    // =============================================

    if(type == FILE_GRIDDEDDATATIMESERIES_TYPE)
    {
      InFileGriddedDataTimeSeries file(fileNameInput);
      logInfo<<"  spline degree = "<<file.splineDegree()<<Log::endl;
      logInfo<<"  node count    = "<<file.nodeCount()<<Log::endl;
      logInfo<<"  data columns  = "<<file.dataCount()<<Log::endl;
      logInfo<<"  sampling      = "<<medianSampling(file.times()).mjd()<<" days"<<Log::endl;
      logInfo<<"  interval      = ["<<file.times().front().dateTimeStr()<<", "<<file.times().back().dateTimeStr()<<"]"<<Log::endl;
      MiscGriddedData::printStatistics(file.grid());

      if(!fileNameOutput.empty())
      {
        std::vector<Matrix> data(file.nodeCount());
        for(UInt i=0; i<file.nodeCount(); i++)
          data.at(i) = file.data(i);
        writeFileGriddedDataTimeSeries(fileNameOutput, file.splineDegree(), file.times(), file.grid(), data);
      }
      return;
    }

    // =============================================

    if(type == FILE_TIMESPLINESGRAVITYFIELD_TYPE)
    {
      InFileTimeSplinesGravityfield file(fileNameInput);
      logInfo<<"  spline degree = "<<file.splineDegree()<<Log::endl;
      logInfo<<"  node count    = "<<file.nodeCount()<<Log::endl;
      logInfo<<"  sampling      = "<<medianSampling(file.times()).mjd()<<" days"<<Log::endl;
      logInfo<<"  interval      = ["<<file.times().front().dateTimeStr()<<", "<<file.times().back().dateTimeStr()<<"]"<<Log::endl;

      if(!fileNameOutput.empty())
      {
        Double GM = DEFAULT_GM, R = DEFAULT_R;
        std::vector<Matrix> cnm(file.nodeCount());
        std::vector<Matrix> snm(file.nodeCount());
        for(UInt i=0; i<file.nodeCount(); i++)
        {
          SphericalHarmonics harm = file.sphericalHarmonics(i);
          GM        = harm.GM();
          R         = harm.R();
          cnm.at(i) = harm.cnm();
          snm.at(i) = harm.snm();
        }
        writeFileTimeSplinesGravityfield(fileNameOutput, GM, R, file.splineDegree(), file.times(), cnm, snm);
      }
      return;
    }

    // =============================================

    if(type == FILE_TIMESPLINESCOVARIANCE_TYPE)
    {
      InFileTimeSplinesCovariance file(fileNameInput);
      logInfo<<"  GM            = "<<file.GM()<<Log::endl;
      logInfo<<"  R             = "<<file.R()<<Log::endl;
      logInfo<<"  minDegree     = "<<file.minDegree()<<Log::endl;
      logInfo<<"  maxDegree     = "<<file.maxDegree()<<Log::endl;
      logInfo<<"  spline degree = "<<file.splineDegree()<<Log::endl;
      logInfo<<"  node count    = "<<file.nodeCount()<<Log::endl;
      logInfo<<"  sampling      = "<<medianSampling(file.times()).mjd()<<" days"<<Log::endl;
      logInfo<<"  interval      = ["<<file.times().front().dateTimeStr()<<", "<<file.times().back().dateTimeStr()<<"]"<<Log::endl;

      if(!fileNameOutput.empty())
      {
        std::vector<Matrix> C(file.nodeCount());
        for(UInt i=0; i<file.nodeCount(); i++)
          C.at(i) = file.covariance(i);
        writeFileTimeSplinesCovariance(fileNameOutput, file.GM(), file.R(), file.minDegree(), file.maxDegree(), file.splineDegree(), file.times(), C);
      }
      return;
    }

    // =============================================

    if(fileNameOutput.empty())
      return;

    // =============================================

    // Standard file types
    // -------------------
    if(type == FILE_ADMITTANCE_TYPE)
    {
      Admittance x;
      readFileAdmittance (fileNameInput,  x);
      writeFileAdmittance(fileNameOutput, x);
    }
    else if(type == FILE_ARCLIST_TYPE)
    {
      std::vector<UInt> arcsInterval;
      std::vector<Time> timesInterval;
      readFileArcList (fileNameInput,  arcsInterval, timesInterval);
      writeFileArcList(fileNameOutput, arcsInterval, timesInterval);
    }
    else if(type == FILE_DOODSONEARTHORIENTATIONPARAMETER_TYPE)
    {
      DoodsonEop x;
      readFileDoodsonEarthOrientationParameter (fileNameInput,  x);
      writeFileDoodsonEarthOrientationParameter(fileNameOutput, x);
    }
    // else if(type == FILE_DOODSONHARMONIC_TYPE) see above
    else if(type == FILE_EARTHORIENTATIONPARAMETER_TYPE)
    {
      Matrix x;
      readFileEarthOrientationParameter (fileNameInput,  x);
      writeFileEarthOrientationParameter(fileNameOutput, x);
    }
    else if(type == FILE_EARTHTIDE_TYPE)
    {
      Matrix kReal, kImag, kPlus;
      Matrix doodson20, doodson21, doodson22;
      Vector ampIp20, ampOp20, ampIp21, ampOp21, amp22;
      Double h2_0, h2_2, l2_0, l2_2, l21_1, l22_1, h21_imag, l21_imag, h22_imag, l22_imag, h3, l3;
      Matrix deformationArg21, deformationArg20;
      Vector dR21_ip, dR21_op, dR20_ip, dR20_op, dT21_ip, dT21_op, dT20_ip, dT20_op;
      readFileEarthTide(fileNameInput, kReal, kImag, kPlus, doodson20, doodson21, doodson22,
                        ampIp20, ampOp20, ampIp21, ampOp21, amp22, h2_0, h2_2, l2_0, l2_2, l21_1, l22_1, h21_imag, l21_imag, h22_imag, l22_imag, h3, l3,
                        deformationArg21, deformationArg20, dR21_ip, dR21_op, dR20_ip, dR20_op, dT21_ip, dT21_op, dT20_ip, dT20_op);
      writeFileEarthTide(fileNameOutput, kReal, kImag, kPlus, doodson20, doodson21, doodson22,
                        ampIp20, ampOp20, ampIp21, ampOp21, amp22, h2_0, h2_2, l2_0, l2_2, l21_1, l22_1, h21_imag, l21_imag, h22_imag, l22_imag, h3, l3,
                        deformationArg21, deformationArg20, dR21_ip, dR21_op, dR20_ip, dR20_op, dT21_ip, dT21_op, dT20_ip, dT20_op);
    }
    else if(type == FILE_GNSSANTENNADEFINITION_TYPE)
    {
      std::vector<GnssAntennaDefinitionPtr> x;
      readFileGnssAntennaDefinition (fileNameInput,  x);
      writeFileGnssAntennaDefinition(fileNameOutput, x);
    }
    else if(type == FILE_GNSSSIGNALBIAS_TYPE)
    {
      GnssSignalBias x;
      readFileGnssSignalBias (fileNameInput,  x);
      writeFileGnssSignalBias(fileNameOutput, x);
    }
    else if(type == FILE_GNSSSTATIONINFO_TYPE)
    {
      std::vector<GnssStationInfo> x;
      readFileGnssStationInfo (fileNameInput,  x);
      writeFileGnssStationInfo(fileNameOutput, x);
    }
    // else if(type == FILE_GRIDRECTANGULAR_TYPE)  see above
    // else if(type == FILE_MATRIX_TYPE)  see above
    else if(type == FILE_MEANPOLARMOTION_TYPE)
    {
      MeanPolarMotion x;
      readFileMeanPolarMotion (fileNameInput,  x);
      writeFileMeanPolarMotion(fileNameOutput, x);
    }
    else if(type == FILE_NORMALEQUATION_TYPE)
    {
      throw(Exception("Conversion of FILE_NORMALEQUATION_TYPE not implemented."));
    }
    else if(type == FILE_OCEANPOLETIDE_TYPE)
    {
      SphericalHarmonics harmReal, harmImag;
      readFileOceanPoleTide (fileNameInput,  harmReal, harmImag);
      writeFileOceanPoleTide(fileNameOutput, harmReal, harmImag);
    }
    else if(type == FILE_PARAMETERNAME_TYPE)
    {
      std::vector<ParameterName> x;
      readFileParameterName (fileNameInput,  x);
      writeFileParameterName(fileNameOutput, x);
    }
    else if(type == FILE_POLYGON_TYPE)
    {
      std::vector<Polygon> x;
      readFilePolygon (fileNameInput,  x);
      writeFilePolygon(fileNameOutput, x);
    }
    // else if(type == FILE_INSTRUMENT_TYPE)  see above
    else if(type == FILE_SATELLITEMODEL_TYPE)
    {
      std::vector<SatelliteModelPtr> x;
      readFileSatelliteModel (fileNameInput,  x);
      writeFileSatelliteModel(fileNameOutput, x);
    }
    else if(type == FILE_POTENTIALCOEFFICIENTS_TYPE)
    {
      SphericalHarmonics x;
      readFileSphericalHarmonics (fileNameInput,  x);
      writeFileSphericalHarmonics(fileNameOutput, x);
    }
    // else if(type == FILE_STRINGLIST_TYPE)  not really GROOPS format
    else if(type == FILE_TIDEGENERATINGPOTENTIAL_TYPE)
    {
      TideGeneratingPotential x;
      readFileTideGeneratingPotential (fileNameInput,  x);
      writeFileTideGeneratingPotential(fileNameOutput, x);
    }
    // else if(type == FILE_TIMESPLINESGRAVITYFIELD_TYPE) see above
    // else if(type == FILE_TIMESPLINESCOVARIANCE_TYPE)   see above
    // else if(type == FILE_VALUEGRID_TYPE)               see above
    // else if(type == FILE_VARIATIONALEQUATION_TYPE)     see above
    else
      throw(Exception("Conversion of file type not implemented."));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

