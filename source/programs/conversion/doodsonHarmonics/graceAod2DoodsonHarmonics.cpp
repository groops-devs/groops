/***********************************************/
/**
* @file graceAod2DoodsonHarmonics.cpp
*
* @brief Convert AOD1B tides into DoodsonHarmonics.
*
* @author Torsten Mayer-Guerr
* @date 2017-05-04
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the atmospheric and ocean tidal products (AOD1B)
from the GRACE SDS format into \configFile{outputfileDoodsonHarmonics}{doodsonHarmonic}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileTideGeneratingPotential.h"

/***** CLASS ***********************************/

/** @brief Convert AOD1B into DoodsonHarmonics.
* @ingroup programsConversionGroup */
class GraceAod2DoodsonHarmonics
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceAod2DoodsonHarmonics, SINGLEPROCESS, "Convert AOD1B tides into DoodsonHarmonics.", Conversion, DoodsonHarmonics)
GROOPS_RENAMED_PROGRAM(GraceCsrAod2DoodsonHarmonics, GraceAod2DoodsonHarmonics, date2time(2020, 6, 14))

/***********************************************/

void GraceAod2DoodsonHarmonics::run(Config &config)
{
  try
  {
    FileName fileNameDoodson, fileNameTGP;
    std::vector<FileName> fileNameInput;

    readConfig(config, "outputfileDoodsonHarmonics",       fileNameDoodson, Config::MUSTSET, "", "" );
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP,     Config::MUSTSET, "{groopsDataDir}/tides/generatingTide_HW95.txt", "to compute Xi phase correction");
    readConfig(config, "inputfile",                        fileNameInput,   Config::MUSTSET, "", "" );
    if(isCreateSchema(config)) return;

    // ======================================================

    // read tide generating potential (TGP)
    // ------------------------------------
    logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
    TideGeneratingPotential tgp;
    readFileTideGeneratingPotential(fileNameTGP, tgp);

    // ======================================================

    Double GM = DEFAULT_GM;
    Double R  = DEFAULT_R;
    std::vector<Doodson>  doodson;
    std::vector<Matrix>   cnmCos, snmCos, cnmSin, snmSin;

    for(UInt i=0; i<fileNameInput.size(); i++)
    {
      logStatus<<"read file <"<<fileNameInput.at(i)<<">"<<Log::endl;
      InFile file(fileNameInput.at(i));

      // Header
      // ------
      std::string line;
      std::string name;
      UInt dataCount = 0;
      UInt version = 9999;
      UInt degree  = 180;
      for(;;)
      {
        std::getline(file, line);
        if(file.eof())
          break;

        if(line.find("PARTIAL TIDE")==0)
          name = String::trim(line.substr(31,4));

        if(line.find("CONSTANT GM")==0)
          GM = String::toDouble(line.substr(31, 22));

        if(line.find("CONSTANT A")==0)
          R = String::toDouble(line.substr(31, 22));

        if(line.find("MAXIMUM DEGREE")==0)
          degree = String::toInt(line.substr(32, 3));

        if(line.find("NUMBER OF DATA SETS")==0)
          dataCount = String::toInt(line.substr(31, 22));

        if(line.find("SOFTWARE VERSION")==0)
          version = String::toInt(line.substr(50, 2));

        // Header fertig
        if(line.find("END OF HEADER")==0)
          break;
      }

      if(version==9999)
        logWarning<<"found not SOFTWARE VERSION"<<Log::endl;

      if(name.empty())
        logWarning<<"found not PARTIAL TIDE"<<Log::endl;

      if(dataCount!=2)
        logWarning<<"NUMBER OF DATA SETS not found or not 2: "<<dataCount<<Log::endl;

      doodson.push_back(Doodson(name));

      for(UInt k=0; k<dataCount; k++)
      {
        // Data Header
        std::getline(file, line);
        std::string type = line.substr(41,3);

        Matrix cnm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        Matrix snm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        for(UInt n=0; n<=degree; n++)
          for(UInt m=0; m<=n; m++)
          {
            UInt n2,m2;
            file>>n2>>m2;
            file>>cnm(n2,m2)>>snm(n2,m2);
            // rest of the line
            std::getline(file, line);
          }

        if(type=="cos")
        {
          cnmCos.push_back(cnm);
          snmCos.push_back(snm);
        }
        else if(type=="sin")
        {
          cnmSin.push_back(cnm);
          snmSin.push_back(snm);
        }
        else
        {
          logError<<"Unknown type: "<<type<<Log::endl;
        }
      }
    }

    // ======================================================

    auto cnmCosOld = cnmCos;
    auto snmCosOld = snmCos;
    auto cnmSinOld = cnmSin;
    auto snmSinOld = snmSin;
    for(UInt i=0; i<doodson.size(); i++)
    {
      const Double xi = tgp.xi(doodson.at(i));
      cnmCos.at(i)  = cnmCosOld.at(i) * cos(xi) - cnmSinOld.at(i) * sin(xi);
      snmCos.at(i)  = snmCosOld.at(i) * cos(xi) - snmSinOld.at(i) * sin(xi);
      cnmSin.at(i)  = cnmCosOld.at(i) * sin(xi) + cnmSinOld.at(i) * cos(xi);
      snmSin.at(i)  = snmCosOld.at(i) * sin(xi) + snmSinOld.at(i) * cos(xi);
    }

    // ======================================================

    // Write tides
    // -----------
    logStatus<<"write tides to <"<<fileNameDoodson<<">"<<Log::endl;
    writeFileDoodsonHarmonic(fileNameDoodson, DoodsonHarmonic(GM, R, doodson, cnmCos, snmCos, cnmSin, snmSin));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
