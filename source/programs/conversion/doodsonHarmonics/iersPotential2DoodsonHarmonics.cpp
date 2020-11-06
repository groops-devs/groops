/***********************************************/
/**
* @file iersPotential2DoodsonHarmonics.cpp
*
* @brief Read ocean tide file in IERS format.
*
* @author Torsten Mayer-Guerr
* @date 2019-11-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read ocean tide file in IERS format.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileDoodsonHarmonic.h"

/***** CLASS ***********************************/

/** @brief Read ocean tide file in IERS format.
* @ingroup programsConversionGroup */
class IersPotential2DoodsonHarmonics
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(IersPotential2DoodsonHarmonics, SINGLEPROCESS, "Read ocean tide file in IERS format", Conversion, DoodsonHarmonics)

/***********************************************/

void IersPotential2DoodsonHarmonics::run(Config &config)
{
  try
  {
    FileName  fileNameOut;
    FileName  fileNameIn, fileNameTGP;
    UInt      countHeader;
    Double    GM, R;
    UInt      minDegree, maxDegree = INFINITYDEGREE;

    readConfig(config, "outputfileDoodsonHarmoncis", fileNameOut, Config::MUSTSET, "",  "");
    readConfig(config, "inputfile",                  fileNameIn,  Config::MUSTSET, "",  "");
    readConfig(config, "headerLines",                countHeader, Config::MUSTSET, "4", "skip number of header lines");
    readConfig(config, "minDegree",                  minDegree,   Config::DEFAULT, "0", "");
    readConfig(config, "maxDegree",                  maxDegree,   Config::MUSTSET, "",  "");
    readConfig(config, "GM",                         GM,          Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                          R,           Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    // ==============================

    logStatus<<"read file from <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);

    // skip headerlines
    std::string line;
    for(UInt i=0; i<countHeader; i++)
      std::getline(file, line);

    std::vector<Doodson> doodson;
    std::vector<Matrix>  cnmCos, cnmSin, snmCos, snmSin;
    for(;;)
    {
      std::getline(file, line);
      if(file.eof())
        break;
      if(line.empty())
        continue;

      // Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-11)
      // Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100)
      // (long period from FES2002 up to (50,50) + equilibrium Om1/Om2, atmospheric tide NOT included)
      // Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
      // 255.555 M2    2   1 -12.07164  -4.38919    -3.09008   1.50139
      std::stringstream ss(line);
      std::string doodstring, name;
      UInt   n, m;
      Double cPlus, sPlus, cMinus, sMinus;
      ss>>doodstring>>name>>n>>m>>cPlus>>sPlus>>cMinus>>sMinus;
      if(doodstring.size() == 6)
        doodstring = '0'+doodstring;

      // Constituent already exist?
      UInt idx = std::distance(doodson.begin(), std::find(doodson.begin(), doodson.end(), Doodson(doodstring)));
      if(idx >= doodson.size())
      {
        doodson.push_back(Doodson(doodstring));
        logStatus<<doodstring<<" "<<name<<Log::endl;
        cnmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        cnmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
      }

      // iers conventions 2010, eq. (6.15)
      if((n>=minDegree) && (n<=maxDegree))
      {
        cnmCos.at(idx)(n,m) =  (cPlus + cMinus) * 1e-11;
        snmCos.at(idx)(n,m) =  (sPlus - sMinus) * 1e-11;
        cnmSin.at(idx)(n,m) =  (sPlus + sMinus) * 1e-11;
        snmSin.at(idx)(n,m) = -(cPlus - cMinus) * 1e-11;
      }
    }

    // write results
    // -------------
    logStatus<<"write tides to <"<<fileNameOut<<">"<<Log::endl;
    writeFileDoodsonHarmonic(fileNameOut, DoodsonHarmonic(GM, R, doodson, cnmCos, snmCos, cnmSin, snmSin));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
