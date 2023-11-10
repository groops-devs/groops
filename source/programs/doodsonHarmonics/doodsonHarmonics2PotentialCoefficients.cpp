/***********************************************/
/**
* @file doodsonHarmonics2PotentialCoefficients.cpp
*
* @brief Write a lot of files with potential coefficients from harmonic tide model.
*
* @author Torsten Mayer-Guerr
* @date 2007-10-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
The \configFile{inputfileDoodsonHarmonics}{doodsonHarmonic} contains a Fourier series of a time variable
gravitational potential at specific tidal frequencies (tides)
\begin{equation}
V(\M x,t) = \sum_{f} V_f^c(\M x)\cos(\theta_f(t)) + V_f^s(\M x)\sin(\theta_f(t)),
\end{equation}
where $V_f^c(\M x)$ and $V_f^s(\M x)$ are spherical harmonics expansions.
If set the expansions are limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly. The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.

The \configFile{outputfilePotentialCoefficients}{potentialCoefficients} is not a single file but a series of files.
For each spherical harmonics expansion $V_f^c(\M x)$ and $V_f^s(\M x)$ a separate file is created
where the variables \config{variableLoopName}, \config{variableLoopDoodson}, \config{variableLoopCosSin} are set accordingly.
The file name should contain these variables, e.g. \verb|coeff.{name}.{doodson}.{cossin}.gfc|.

If \config{applyXi} the Doodson-Warburg phase correction (see IERS conventions) is applied to the cos/sin
potentialCoefficients before.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileTideGeneratingPotential.h"

/***** CLASS ***********************************/

/** @brief Write a lot of files with PotentialCoefficients from harmonic tide model.
* @ingroup programsGroup */
class DoodsonHarmonics2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonics2PotentialCoefficients, SINGLEPROCESS, "write a lot of files with potential coefficients from harmonic tide model.", DoodsonHarmonics, PotentialCoefficients)

/***********************************************/

void DoodsonHarmonics2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut, fileNameIn;
    FileName    fileNameTGP;
    std::string nameName, nameDoodson, nameCosSin, nameIndex, nameCount;
    UInt        minDegree, maxDegree = INFINITYDEGREE;
    Double      GM, R;
    Bool        applyXi;

    readConfig(config, "outputfilePotentialCoefficients",  fileNameOut, Config::MUSTSET,  "coeff.{name}.{doodson}.{cossin}.gfc", "");
    readConfig(config, "variableLoopName",                 nameName,    Config::OPTIONAL, "name",    "variable with darwins's name of each constituent");
    readConfig(config, "variableLoopDoodson",              nameDoodson, Config::OPTIONAL, "doodson", "variable with doodson code of each constituent");
    readConfig(config, "variableLoopCosSin",               nameCosSin,  Config::OPTIONAL, "cossin",  "variable with 'cos' or 'sin' of each constituent");
    readConfig(config, "variableLoopIndex",                nameIndex,   Config::OPTIONAL, "",        "variable with index of each constituent (starts with zero)");
    readConfig(config, "variableLoopCount",                nameCount,   Config::OPTIONAL, "",        "variable with total number of constituents");
    readConfig(config, "inputfileDoodsonHarmonics",        fileNameIn,  Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP, Config::OPTIONAL, "{groopsDataDir}/tides/generatingTide_HW95.txt", "to compute Xi phase correction");
    readConfig(config, "minDegree",                        minDegree,   Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",                        maxDegree,   Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                               GM,          Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                R,           Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "applyXi",                          applyXi,     Config::DEFAULT,  "0", "apply Doodson-Warburg phase correction (see IERS conventions)");
    if(isCreateSchema(config)) return;

    // read tide generating potential (TGP)
    // ------------------------------------
    TideGeneratingPotential tgp;
    if(applyXi)
    {
      if(fileNameTGP.empty())
        throw(Exception("Need TideGeneratingPotential to compute xi phase correction"));
      logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
      readFileTideGeneratingPotential(fileNameTGP, tgp);
    }

    // read ocean tide file
    // --------------------
    logStatus<<"read doodson harmonics file <"<<fileNameIn<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameIn, d);

    // write tides
    // -----------
    VariableList varList;
    if(!nameName.empty())    varList.undefineVariable(nameName);
    if(!nameDoodson.empty()) varList.undefineVariable(nameDoodson);
    if(!nameCosSin.empty())  varList.undefineVariable(nameCosSin);
    if(!nameIndex.empty())   varList.undefineVariable(nameIndex);
    if(!nameCount.empty())   varList.setVariable(nameCount, static_cast<Double>(d.doodson.size()));
    logStatus<<"writing potential coefficients to files <"<<fileNameOut(varList)<<">"<<Log::endl;

    for(UInt i=0; i<d.doodson.size(); i++)
    {
      logStatus<<"  "<<d.doodson.at(i).code()<<" "<<d.doodson.at(i).name()<<Log::endl;
      Double xi = 0.;
      if(applyXi)
        xi = tgp.xi(d.doodson.at(i));

      const Matrix cnmCos = d.cnmCos.at(i) * cos(xi) - d.cnmSin.at(i) * sin(xi);
      const Matrix snmCos = d.snmCos.at(i) * cos(xi) - d.snmSin.at(i) * sin(xi);
      const Matrix cnmSin = d.cnmCos.at(i) * sin(xi) + d.cnmSin.at(i) * cos(xi);
      const Matrix snmSin = d.snmCos.at(i) * sin(xi) + d.snmSin.at(i) * cos(xi);

      const SphericalHarmonics harmCos(d.GM, d.R, cnmCos, snmCos);
      const SphericalHarmonics harmSin(d.GM, d.R, cnmSin, snmSin);

      if(!nameName.empty())    varList.setVariable(nameName, d.doodson.at(i).name());
      if(!nameDoodson.empty()) varList.setVariable(nameDoodson, d.doodson.at(i).code());
      if(!nameIndex.empty())   varList.setVariable(nameIndex, static_cast<Double>(i));
      if(!nameCosSin.empty())  varList.setVariable(nameCosSin, "cos");
      writeFileSphericalHarmonics(fileNameOut(varList), harmCos.get(maxDegree, minDegree, GM, R));
      if(!nameCosSin.empty())  varList.setVariable(nameCosSin, "sin");
      writeFileSphericalHarmonics(fileNameOut(varList), harmSin.get(maxDegree, minDegree, GM, R));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
