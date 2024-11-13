/***********************************************/
/**
* @file doodsonAdmittance2SupplementaryFiles.cpp
*
* @brief Write ascii doodson and admittance matrix.
*
* @author Torsten Mayer-Guerr
* @date 2022-07-22
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
The publication of an ocean tide model includes not only the atlas
in the form of spherical harmonics coefficients,
but also the matrix of Doodson multipliers (\config{outputfileDoodsonMatrix})
and the \config{outputfileAdmittanceMatrix}.

The \config{outputfileMajorTideList} contains the \config{fileNames}
for each contituent.
The required information is taken from the
\configFile{inputfileAdmittance}{admittance}.

See also \program{DoodsonHarmonics2PotentialCoefficients}.
)";

/***********************************************/

#include "base/import.h"
#include "base/doodson.h"
#include "inputOutput/file.h"
#include "files/fileMatrix.h"
#include "files/fileAdmittance.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Write ascii doodson and admittance matrix.
* @ingroup programsConversionGroup */
class DoodsonAdmittance2SupplementaryFiles
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonAdmittance2SupplementaryFiles, SINGLEPROCESS, "write ascii doodson and admittance matrix", Conversion, DoodsonHarmonics)
GROOPS_RENAMED_PROGRAM(Admittance2Iers, DoodsonAdmittance2SupplementaryFiles, date2time(2023, 7, 1))

/***********************************************/

void DoodsonAdmittance2SupplementaryFiles::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameTideList, fileNameDoodson, fileNameMatrix;
    FileName fileNames;
    FileName fileNameAdmittance;

    readConfig(config, "outputfileMajorTideList",    fileNameTideList,   Config::OPTIONAL, "fileList.txt",   "");
    readConfig(config, "fileNames",                  fileNames,          Config::OPTIONAL, "model_{doodson}_{name}_{cossin}.gfc", "template for fileList, variables: doodson, name, cossin");
    readConfig(config, "outputfileDoodsonMatrix",    fileNameDoodson,    Config::OPTIONAL, "doodson.txt",    "");
    readConfig(config, "outputfileAdmittanceMatrix", fileNameMatrix,     Config::OPTIONAL, "admittance.txt", "");
    readConfig(config, "inputfileAdmittance",        fileNameAdmittance, Config::MUSTSET,  "", "interpolation of minor constituents");
    if(isCreateSchema(config)) return;

    logStatus<<"read admittance file <"<<fileNameAdmittance<<">"<<Log::endl;
    Admittance admit;
    readFileAdmittance(fileNameAdmittance, admit);

    if(!fileNameTideList.empty())
    {
      logStatus<<"write major tide list to file <"<<fileNameTideList<<">"<<Log::endl;
      OutFile file(fileNameTideList);
      for(const auto &dood : admit.doodsonMajor)
      {
        VariableList varList;
        varList.setVariable("doodson", dood.code());
        varList.setVariable("name",    dood.name());
        varList.setVariable("cossin", "cos");
        file<<fileNames(varList).str()<<" ";
        varList.setVariable("cossin", "sin");
        file<<fileNames(varList).str()<<std::endl;
      }
    }

    if(!fileNameDoodson.empty())
    {
      logStatus<<"write doodson coefficients to file <"<<fileNameDoodson<<">"<<Log::endl;
      OutFile file(fileNameDoodson);
      for(const auto &dood : admit.doodsonMinor)
      {
        for(UInt i=0; i<6; i++)
          file<<" "<<dood.d[i]%"%3i"s;
        file<<std::endl;
      }
    }

    if(!fileNameMatrix.empty())
    {
      logStatus<<"write admittance matrix to file <"<<fileNameMatrix<<">"<<Log::endl;
      OutFile file(fileNameMatrix);
      for(UInt i=0; i<admit.admittance.rows(); i++)
      {
        for(UInt k=0; k<admit.admittance.columns(); k++)
          file<<admit.admittance(i,k)%"%16.8e"s;
        file<<std::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
