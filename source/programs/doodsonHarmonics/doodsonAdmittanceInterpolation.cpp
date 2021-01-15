/***********************************************/
/**
* @file doodsonAdmittanceInterpolation.cpp
*
* @brief Visualize the interpolation of the minor tides.
*
* @author Torsten Mayer-Guerr
* @date 2010-05-13
* update 2012-08-02 (DR): write TGP amplitudes in 2nd column
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
To visualize the interpolation of the minor tides.
The output is a \file{matrix}{matrix} with the first column containing the tidal frequency,
the second column is the tide generating amplitude (from \configFile{inputfileTideGeneratingPotential}{tideGeneratingPotential}), and the following
columns the contribution of the major tides to the this tidal frequency as defined in in \configFile{inputfileAdmittance}{admittance}.

\fig{!hb}{0.8}{doodsonAdmittanceInterpolation}{fig:doodsonAdmittanceInterpolation}{Linear interpolation of minor tides in the diurnal band.}
)";

/***********************************************/

#include "base/import.h"
#include "base/doodson.h"
#include "files/fileMatrix.h"
#include "files/fileAdmittance.h"
#include "files/fileTideGeneratingPotential.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Visualize the interpolation of the minor tides.
* @ingroup programsGroup */
class DoodsonAdmittanceInterpolation
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonAdmittanceInterpolation, SINGLEPROCESS, "visualize the interpolation of the minor tides.", DoodsonHarmonics)

/***********************************************/

void DoodsonAdmittanceInterpolation::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      outName, admittanceName, tgpName;

    readConfig(config, "outputfile",                       outName,        Config::MUSTSET,  "", "");
    readConfig(config, "inputfileAdmittance",              admittanceName, Config::MUSTSET,  "", "interpolation of minor constituents");
    readConfig(config, "inputfileTideGeneratingPotential", tgpName,        Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    // read admittace file
    // -------------------
    logStatus<<"read admittance file <"<<admittanceName<<">"<<Log::endl;
    Admittance admit;
    readFileAdmittance(admittanceName, admit);

    // read tide generating potential
    // ------------------------------
    std::vector<Double> majorAmpl(admit.doodsonMajor.size(), 1.0);
    std::vector<Double> ampl(admit.doodsonMinor.size(), 1.0);
    if(!tgpName.empty())
    {
      logStatus<<"read tide generating potential <"<<tgpName<<">"<<Log::endl;
      TideGeneratingPotential tgp;
      readFileTideGeneratingPotential(tgpName, tgp);

      for(UInt i=0; i<tgp.size(); i++)
      {
        for(UInt k=0; k<admit.doodsonMajor.size(); k++)
          if(tgp.at(i) == admit.doodsonMajor.at(k))
            majorAmpl.at(k) = tgp.at(i).admit();
        for(UInt k=0; k<admit.doodsonMinor.size(); k++)
          if(tgp.at(i) == admit.doodsonMinor.at(k))
            ampl.at(k) = tgp.at(i).admit();
      }
    }

    // Computing
    // ---------
    logStatus<<"computing interpolation"<<Log::endl;
    Matrix A(admit.doodsonMinor.size(), 2+admit.doodsonMajor.size());
    for(UInt i=0; i<admit.doodsonMinor.size(); i++)
    {
      A(i,0) = admit.doodsonMinor.at(i).frequency()/(2*PI);
      A(i,1) = ampl.at(i);
      for(UInt k=0; k<admit.doodsonMajor.size(); k++)
        A(i,2+k) = admit.admittance(k,i)/ampl.at(i)*majorAmpl.at(k);
    }

    // save results
    // ------------
    logStatus<<"write interpolation to file <"<<outName<<">"<<Log::endl;
    writeFileMatrix(outName, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
