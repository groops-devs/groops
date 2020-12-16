/***********************************************/
/**
* @file gnssEop2IgsErp.cpp
*
* @brief Write GNSS Earth orientation parameters to IGS ERP file format.
**
* @author Sebastian Strasser
* @date 2019-05-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write GNSS Earth orientation parameters to \href{https://files.igs.org/pub/data/format/erp.txt}{IGS ERP file format}.

Requires polar motion, polar motion rate, dUT1 and LOD parameters in the solution
vector \configFile{inputfileSolution}{matrix} and their sigmas in \configFile{inputfileSigmax}{matrix}.
Solution usually comes out of \program{GnssProcessing}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Write GNSS Earth orientation parameters to IGS ERP file format.
* @ingroup programsConversionGroup */
class GnssEop2IgsErp
{
public:
  class Epoch
  {
  public:
    FileName fileNameSolution, fileNameSigmax, fileNameParameterNames, fileNameStationList;
    std::vector<FileName> fileNameTransmitterList;
    Time time;
  };

  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssEop2IgsErp, SINGLEPROCESS, "Write GNSS Earth orientation parameters to IGS ERP file format.", Conversion, Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssEop2IgsErp::Epoch &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileSolution",          var.fileNameSolution,        Config::MUSTSET,  "", "parameter vector");
  readConfig(config, "inputfileSigmax",            var.fileNameSigmax,          Config::MUSTSET,  "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
  readConfig(config, "inputfileParameterNames",    var.fileNameParameterNames,  Config::MUSTSET,  "", "parameter names");
  readConfig(config, "inputfileTransmitterList",   var.fileNameTransmitterList, Config::OPTIONAL, "", "transmitter PRNs used in solution (used for transmitter count)");
  readConfig(config, "inputfileStationList",       var.fileNameStationList,     Config::OPTIONAL, "", "stations used in solution (used for station count)");
  readConfig(config, "time",                       var.time,                    Config::MUSTSET,  "", "reference time for epoch");
  endSequence(config);
  return TRUE;
}

/***********************************************/


void GnssEop2IgsErp::run(Config &config)
{
  try
  {
    FileName fileNameOut;
    std::vector<Epoch> epochs;
    std::vector<std::string> comments;

    readConfig(config, "outputfileIgsErp", fileNameOut, Config::MUSTSET,  "", "IGS ERP file");
    readConfig(config, "epoch",            epochs,      Config::MUSTSET,  "", "e.g. daily solution");
    readConfig(config, "comment",          comments,    Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    std::vector<std::string> lines;
    for(const auto &epoch : epochs)
    {
      logStatus<<"reading solution from <"<<epoch.fileNameSolution<<">"<<Log::endl;
      Vector x;
      readFileMatrix(epoch.fileNameSolution, x);

      logStatus<<"reading standard deviations from <"<<epoch.fileNameSigmax<<">"<<Log::endl;
      Vector sigmax;
      readFileMatrix(epoch.fileNameSigmax, sigmax);

      logStatus<<"reading parameter names from <"<<epoch.fileNameParameterNames<<">"<<Log::endl;
      std::vector<ParameterName> parameterNames;
      readFileParameterName(epoch.fileNameParameterNames, parameterNames);

      std::vector<std::string> transmitterList;
      for(const auto &fileName : epoch.fileNameTransmitterList)
      {
        logStatus<<"reading transmitter list from <"<<fileName<<">"<<Log::endl;
        readFileStringList(fileName, transmitterList);
      }

      logStatus<<"reading station list from <"<<epoch.fileNameStationList<<">"<<Log::endl;
      std::vector<std::string> stationList;
      readFileStringList(epoch.fileNameStationList, stationList);

      Double xp = 0, yp = 0, xpRate = 0, ypRate = 0, ut1 = 0, lod = 0;
      Double sigmaXp = 0, sigmaYp = 0, sigmaXpRate = 0, sigmaYpRate = 0, sigmaUt1 = 0, sigmaLod = 0;
      auto extract = [&x, &sigmax](UInt index, Double factor, Double &value, Double &sigma)
      {
        value = factor*x(index);
        sigma = factor*sigmax(index);
      };
      for(UInt i = 0; i < parameterNames.size(); i++)
        if(parameterNames.at(i).object == "earth")
        {
          if(parameterNames.at(i).type == "polarMotion.xp")
          {
            if(!parameterNames.at(i).temporal.size())
              extract(i, 1e3, xp, sigmaXp); // mas => E-6"
            else if(parameterNames.at(i).temporal.substr(0,6) == "trend.")
              extract(i, 1e3, xpRate, sigmaXpRate); // mas/d => E-6"
          }
          else if(parameterNames.at(i).type == "polarMotion.yp")
          {
            if(!parameterNames.at(i).temporal.size())
              extract(i, 1e3, yp, sigmaYp); // mas => E-6"
            else if(parameterNames.at(i).temporal.substr(0,6) == "trend.")
              extract(i, 1e3, ypRate, sigmaYpRate); // mas/d => E-6"
          }
          else if(parameterNames.at(i).type == "UT1")
          {
            if(!parameterNames.at(i).temporal.size())
              extract(i, 1e4, ut1, sigmaUt1); // ms => .1us
            else if(parameterNames.at(i).temporal.substr(0,6) == "trend.")
              extract(i, 1e4, lod, sigmaLod); // ms/d => .1us/d
          }
        }

      std::stringstream ss;
      ss << epoch.time.mjd()%"%9.2f"s << xp%"%9i"s << yp%"%9i"s << ut1%"%9i"s << lod%"%8i"s << sigmaXp%"%7i"s << sigmaYp%"%7i"s << sigmaUt1%"%7i"s << sigmaLod%"%7i"s
         << stationList.size()%"%5i"s << "    0" << transmitterList.size()%"%5i"s << xpRate%"%7i"s << ypRate%"%7i"s << sigmaXpRate%"%7i"s << sigmaYpRate%"%7i"s;
      lines.push_back(ss.str());
    }

    logStatus<<"writing ERPs to <"<<fileNameOut<<">"<<Log::endl;
    OutFile outfile(fileNameOut);
    outfile << "version 2" << std::endl;
    for(const auto &comment : comments)
      outfile << comment << std::endl;
    outfile << std::string(115, '-') << std::endl;
    outfile << "   MJD       Xpole    Ypole  UT1-UTC     LOD   Xsig   Ysig  UTsig LODsig   Nr   Nf   Nt    Xrt    Yrt Xrtsig Yrtsig" << std::endl;
    outfile << "              e-6\"     e-6\"     e-7s  e-7s/d   e-6\"   e-6\"   e-7s e-7s/d                e-6\"/d e-6\"/d e-6\"/d e-6\"/d" << std::endl;
    outfile << std::string(115, '-') << std::endl;
    for(const auto &line : lines)
      outfile << line << std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
