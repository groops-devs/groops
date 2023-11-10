/***********************************************/
/**
* @file grace2PotentialCoefficients.cpp
*
* @brief read GRACE data.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts potential coefficients from the GRACE SDS format
into \file{potential coefficients file}{potentialCoefficients}.
The program supports file formats for RL04 to RL06.

Within the program, the variables \verb|epochStart|, \verb|epochEnd| and \verb|epochMid|
are populated with the corresponding time-stamps in the file.
These can be used in to \configFile{outputfilePotentialCoefficients}{potentialCoefficients}
to auto-generate the file name.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileSphericalHarmonics.h"

/***** CLASS ***********************************/

/** @brief Read CSR GRACE data.
* @ingroup programsConversionGroup */
class Grace2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Grace2PotentialCoefficients, SINGLEPROCESS, "read GRACE data", Conversion, Grace, PotentialCoefficients)

/***********************************************/

void Grace2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameIn;
    readConfig(config, "outputfilePotentialCoefficients", fileNameOut, Config::MUSTSET, "", "variables: epochStart, epochEnd, epochMid");
    readConfig(config, "inputfile",                       fileNameIn,  Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);

    Double  GM = 0.3986004415e15;
    Double  R  = 0.6378136460e07;
    Matrix  cnm, snm, sigma2cnm, sigma2snm;
    Time    timeStart = date2time(9999,1,1), timeEnd;

    Bool seekR = FALSE;
    Bool seekGM = FALSE;
    Bool seekDegree = FALSE;

    std::string line;
    while(std::getline(file, line))
    {
      if(line.find("dimension") != std::string::npos)
        seekDegree=TRUE;
      if( (line.find("degree") != std::string::npos) && seekDegree )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        UInt degree = INFINITYDEGREE;
        ss>>degree;
        cnm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        snm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2cnm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2snm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        seekDegree = FALSE;
      }

      if(line.find("SHM ")==0)
      {
        UInt degree = String::toInt(line.substr(6, 5));
        cnm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        snm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2cnm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2snm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      }
      if(line.find("EARTH ")==0)
      {
        GM = String::toDouble(line.substr(6, 16));
        R  = String::toDouble(line.substr(23, 16));
      }

      if(line.find("mean_equator_radius") != std::string::npos)
        seekR = TRUE;
      if( (line.find("value") != std::string::npos) && seekR )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        ss>>R;
        seekR = FALSE;
      }

      if(line.find("earth_gravity_param") != std::string::npos)
        seekGM = TRUE;
      if( (line.find("value") != std::string::npos) && seekGM )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        ss>>GM;
        seekGM = FALSE;
      }

      if((line.find("GRCOEF")==0)||(line.find("GRCOF2")==0))
      {
        const UInt n    = String::toInt(line.substr(6, 5));
        const UInt m    = String::toInt(line.substr(11, 5));
        cnm(n,m)        = String::toDouble(line.substr(17, 18));
        snm(n,m)        = String::toDouble(line.substr(36, 18));
        sigma2cnm(n,m)  = String::toDouble(line.substr(55, 10));
        sigma2snm(n,m)  = String::toDouble(line.substr(66, 10));
        sigma2cnm(n,m) *= sigma2cnm(n,m);
        sigma2snm(n,m) *= sigma2snm(n,m);

        const UInt yearStart  = String::toInt(line.substr(77, 4));
        const UInt monthStart = String::toInt(line.substr(81, 2));
        const UInt dayStart   = String::toInt(line.substr(83, 2));
        timeStart = std::min(timeStart, date2time(yearStart, monthStart, dayStart));

        const UInt yearEnd    = String::toInt(line.substr(91, 4));
        const UInt monthEnd   = String::toInt(line.substr(95, 2));
        const UInt dayEnd     = String::toInt(line.substr(97, 2));
        timeEnd = std::max(timeEnd, date2time(yearEnd, monthEnd, dayEnd));
      }
    }

    VariableList fileNameVariableList;
    fileNameVariableList.setVariable("epochStart", timeStart.mjd());
    fileNameVariableList.setVariable("epochEnd",   timeEnd.mjd());
    fileNameVariableList.setVariable("epochMid",   (0.5*timeStart+0.5*timeEnd).mjd());
    logStatus<<"writing potential coefficients to file <"<<fileNameOut(fileNameVariableList)<<">"<<Log::endl;
    writeFileSphericalHarmonics(fileNameOut(fileNameVariableList), SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
