/***********************************************/
/**
* @file iersWaterHeight2DoodsonHarmonics.cpp
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
#include "files/fileTideGeneratingPotential.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Read ocean tide file in IERS format.
* @ingroup programsConversionGroup */
class IersWaterHeight2DoodsonHarmonics
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(IersWaterHeight2DoodsonHarmonics, SINGLEPROCESS, "Read ocean tide file in IERS format", Conversion, DoodsonHarmonics)

/***********************************************/

void IersWaterHeight2DoodsonHarmonics::run(Config &config)
{
  try
  {
    FileName  fileNameOut;
    FileName  fileNameIn, fileNameTGP;
    UInt      countHeader;
    KernelPtr kernel;
    Double    factor;
    Double    GM, R;
    UInt      minDegree, maxDegree = INFINITYDEGREE;

    readConfig(config, "outputfileDoodsonHarmoncis",       fileNameOut, Config::MUSTSET, "",  "");
    readConfig(config, "inputfile",                        fileNameIn,  Config::MUSTSET, "",  "");
    readConfig(config, "headerLines",                      countHeader, Config::MUSTSET, "3", "skip number of header lines");
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP, Config::MUSTSET, "{groopsDataDir}/tides/generatingTide_HW95.txt", "to compute Xi phase correction");
    readConfig(config, "kernel",                           kernel,      Config::MUSTSET, "waterHeight", "data type of input values");
    readConfig(config, "factor",                           factor,      Config::DEFAULT, "0.01", "to convert in SI units");
    readConfig(config, "minDegree",                        minDegree,   Config::DEFAULT, "0", "");
    readConfig(config, "maxDegree",                        maxDegree,   Config::MUSTSET, "",  "");
    readConfig(config, "GM",                               GM,          Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                R,           Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    // ==============================

    // read tide generating potential (TGP)
    // ------------------------------------
    logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
    TideGeneratingPotential tgp;
    readFileTideGeneratingPotential(fileNameTGP, tgp);

    // Koeffizienten einlesen
    // ----------------------
    logStatus<<"read file from <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);

    // skip headerlines
    std::string line;
    for(UInt i=0; i<countHeader; i++)
      std::getline(file, line);

    std::vector<Doodson> doodson;
    std::vector<Matrix>  cnmCos, cnmSin, snmCos, snmSin;
    std::vector<Double>  xi;
    const Vector kn = factor * R/GM * kernel->coefficients(Vector3d(R,0,0), maxDegree);
    for(;;)
    {
      std::getline(file, line);
      if(file.eof())
        break;
      if(line.empty())
        continue;

      // Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100) in cm
      // (long period from FES2002 up to (50,50) + equilibrium Om1/Om2, atmospheric tide NOT included)
      // Doodson Darw  n   m    Csin+     Ccos+       Csin-     Ccos-       C+   eps+      C-   eps-
      // 255.555 M2    2   1 -0.991591 -0.360545   -0.253804  0.123315   1.0551 250.019 0.2822 295.914
      std::stringstream ss(line);
      std::string doodstring, name;
      UInt   n, m;
      Double cSinPlus, cCosPlus, cSinMinus, cCosMinus;
      Double cBarPlus, epsPlus, cBarMinus, epsMinus;
      ss>>doodstring>>name>>n>>m>>cCosPlus>>cSinPlus>>cCosMinus>>cSinMinus>>cBarPlus>>epsPlus>>cBarMinus>>epsMinus;
      if(doodstring.size() == 6)
        doodstring = '0'+doodstring;

      // Constituent already exist?
      UInt idx = std::distance(doodson.begin(), std::find(doodson.begin(), doodson.end(), Doodson(doodstring)));
      if(idx >= doodson.size())
      {
        doodson.push_back(Doodson(doodstring));
        cnmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        cnmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        xi.push_back(tgp.xi(doodson.at(idx)));
        logStatus<<doodstring<<" "<<name<<std::string(6ul-std::min(static_cast<std::size_t>(6), name.size()), ' ')<<"xi = "<<RAD2DEG*xi.back()<<Log::endl;
      }

      if((n>=minDegree) && (n<=maxDegree))
      {
        // variant 1: iers conventions 2010, eq. (6.20) and (6.21)
        const Double cPlus  = cCosPlus  * std::cos(xi.at(idx)) + cSinPlus  * std::sin(xi.at(idx));
        const Double cMinus = cCosMinus * std::cos(xi.at(idx)) + cSinMinus * std::sin(xi.at(idx));
        const Double sPlus  = cSinPlus  * std::cos(xi.at(idx)) - cCosPlus  * std::sin(xi.at(idx));
        const Double sMinus = cSinMinus * std::cos(xi.at(idx)) - cCosMinus * std::sin(xi.at(idx));

        // variant 2: iers conventions 2010, eq. (6.21)
//         const Double cPlus  = cBarPlus  * std::sin(DEG2RAD*epsPlus  + xi.at(idx));
//         const Double cMinus = cBarMinus * std::sin(DEG2RAD*epsMinus + xi.at(idx));
//         const Double sPlus  = cBarPlus  * std::cos(DEG2RAD*epsPlus  + xi.at(idx));
//         const Double sMinus = cBarMinus * std::cos(DEG2RAD*epsMinus + xi.at(idx));

        // iers conventions 2010, eq. (6.15)
        cnmCos.at(idx)(n,m) =  kn(n) * (cPlus + cMinus);
        snmCos.at(idx)(n,m) =  kn(n) * (sPlus - sMinus);
        cnmSin.at(idx)(n,m) =  kn(n) * (sPlus + sMinus);
        snmSin.at(idx)(n,m) = -kn(n) * (cPlus - cMinus);
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
