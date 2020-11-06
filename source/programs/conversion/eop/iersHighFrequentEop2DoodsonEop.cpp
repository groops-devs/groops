/***********************************************/
/**
* @file iersHighFrequentEop2DoodsonEop.cpp
*
* @brief Read Diurnal and Subdiurnal Earth Orientation variations.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-15
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read Diurnal and Subdiurnal Earth Orientation variations according to updated IERS 2010 conventions
and write them as \configFile{outputfileDoodsonEOP}{doodsonEarthOrientationParameter}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileDoodsonEarthOrientationParameter.h"

/***** CLASS ***********************************/

/** @brief Read Diurnal and Subdiurnal Earth Orientation variations.
* @ingroup programsConversionGroup */
class IersHighFrequentEop2DoodsonEop
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(IersHighFrequentEop2DoodsonEop, SINGLEPROCESS, "Read Diurnal and Subdiurnal Earth Orientation variations.", Conversion)

/***********************************************/

static Bool stripComments(std::istream &stream)
{
  try
  {
    char c;
    if(!(stream>>c))
      return FALSE;
    stream.putback(c);
    if(c != '#')
      return TRUE;
    std::string line;
    std::getline(stream, line); // skip rest of line
    return stripComments(stream);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void IersHighFrequentEop2DoodsonEop::run(Config &config)
{
  try
  {
    FileName fileNameIn, fileNameDoodsonEop;

    readConfig(config, "outputfileDoodsonEOP", fileNameDoodsonEop, Config::MUSTSET, "", "");
    readConfig(config, "inputfile",            fileNameIn,         Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    std::vector<Doodson> doodson;
    std::vector<Double>  xpCos, xpSin, ypCos, ypSin, ut1Cos, ut1Sin, lodCos, lodSin;

    while(stripComments(file))
    {
      std::string name, doodName;
      Int         arg;
      Double      period;
      Double      xs, xc, ys, yc, us, uc, ls, lc;
      file>>name>>arg>>arg>>arg>>arg>>arg>>arg>>doodName>>period>>xs>>xc>>ys>>yc>>us>>uc>>ls>>lc;

      doodson.push_back(Doodson(doodName));
      xpCos.push_back(1e-6*xc); // micro arc sec -> arc sec
      xpSin.push_back(1e-6*xs);
      ypCos.push_back(1e-6*yc);
      ypSin.push_back(1e-6*ys);
      ut1Cos.push_back(1e-6*uc); // micro seconds -> seconds
      ut1Sin.push_back(1e-6*us);
      lodCos.push_back(1e-6*lc);
      lodSin.push_back(1e-6*ls);
    }

    DoodsonEop doodsonEop;
    doodsonEop.doodson = doodson;
    doodsonEop.coeff   = Matrix(doodson.size(), 8);
    for(UInt i=0; i<doodson.size(); i++)
    {
      doodsonEop.coeff(i, 0) = xpCos.at(i);
      doodsonEop.coeff(i, 1) = xpSin.at(i);
      doodsonEop.coeff(i, 2) = ypCos.at(i);
      doodsonEop.coeff(i, 3) = ypSin.at(i);
      doodsonEop.coeff(i, 4) = ut1Cos.at(i);
      doodsonEop.coeff(i, 5) = ut1Sin.at(i);
      doodsonEop.coeff(i, 6) = lodCos.at(i);
      doodsonEop.coeff(i, 7) = lodSin.at(i);
    }

    // Save to file
    // ------------
    logStatus<<"write doodson EOP file <"<<fileNameDoodsonEop<<">"<<Log::endl;
    writeFileDoodsonEarthOrientationParameter(fileNameDoodsonEop, doodsonEop);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
