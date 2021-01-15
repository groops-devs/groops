/***********************************************/
/**
* @file goceXml2Gradiometer.cpp
*
* @brief Read ESA XML GOCE Data.
*
* @author Torsten Mayer-Guerr
* @date 2010-10-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read ESA XML GOCE Data.
The \config{outputfileGradiometer} is written as \file{instrument file (GRADIOMETER)}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/xml.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read ESA XML GOCE Data.
* @ingroup programsConversionGroup */
class GoceXml2Gradiometer
{
  void readFileGoceGradiometer(const FileName &fileName, GradiometerArc &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GoceXml2Gradiometer, SINGLEPROCESS, "read ESA XML GOCE Data", Conversion, Instrument)

/***********************************************/

void GoceXml2Gradiometer::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> fileName;

    readConfig(config, "outputfileGradiometer", outName,  Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",             fileName, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    GradiometerArc arc;
    for(UInt i=0; i<fileName.size(); i++)
    {
      logStatus<<"read file <"<<fileName.at(i)<<">"<<Log::endl;
      readFileGoceGradiometer(fileName.at(i), arc);
    }

    logStatus<<"write gradiometer file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GoceXml2Gradiometer::readFileGoceGradiometer(const FileName &fileName, GradiometerArc &arc)
{
  try
  {
    GradiometerEpoch epoch;

    InFile     file(fileName);
    XmlNodePtr rootNode = XmlNode::read(file);
    XmlNodePtr dataNode = getChild(rootNode, "Data_Block", TRUE);
    XmlNodePtr ggtNode  = getChild(dataNode, "EGG_GGT_DS", TRUE);

    UInt epochCount = childCount(ggtNode, "EGG_GGT_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(ggtNode, "EGG_GGT_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr tnsNode = getChild(epochNode, "Gravity_Grad_Tensor", TRUE);
      std::string tnsStr;
      readXml(tnsNode, "U_G", tnsStr, TRUE);

      std::stringstream ss(tnsStr);
      ss>>epoch.gravityGradient.xx();
      ss>>epoch.gravityGradient.yy();
      ss>>epoch.gravityGradient.zz();
      ss>>epoch.gravityGradient.xy();
      ss>>epoch.gravityGradient.xz();
      ss>>epoch.gravityGradient.yz();

      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
