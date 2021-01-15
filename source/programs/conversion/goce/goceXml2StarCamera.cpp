/***********************************************/
/**
* @file goceXml2StarCamera.cpp
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
)";

/***********************************************/

#include "programs/program.h"
#include "parser/xml.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read ESA XML GOCE Data.
* @ingroup programsConversionGroup */
class GoceXml2StarCamera
{
  void readFileGoceStarCamera(const FileName &fileName, StarCameraArc &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GoceXml2StarCamera, SINGLEPROCESS, "read ESA XML GOCE Data", Conversion, Instrument)

/***********************************************/

void GoceXml2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> fileName;

    readConfig(config, "outputfileStarCamera", outName, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile", fileName, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;
    for(UInt i=0; i<fileName.size(); i++)
    {
      logStatus<<"read file '"<<fileName.at(i)<<"'"<<Log::endl;
      readFileGoceStarCamera(fileName.at(i), arc);
    }

    // Daten speichern
    // ---------------
    logStatus<<"write gradiometer file <"<<outName<<">"<<Log::endl;
    std::list<Arc> arcList; arcList.push_back(arc);
    InstrumentFile::write(outName, arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GoceXml2StarCamera::readFileGoceStarCamera(const FileName &fileName, StarCameraArc &arc)
{
  try
  {
    StarCameraEpoch epoch;

    InFile     file(fileName);
    XmlNodePtr rootNode = XmlNode::read(file);
    XmlNodePtr dataNode = getChild(rootNode, "Data_Block", TRUE);
    XmlNodePtr ggtNode  = getChild(dataNode, "EGG_IAQ_DS", TRUE);

    UInt epochCount = childCount(ggtNode, "EGG_IAQ_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(ggtNode, "EGG_IAQ_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr tnsNode = getChild(epochNode, "Corr_Quat", TRUE);
      std::string tnsStr;
      readXml(tnsNode, "Q_Grad", tnsStr, TRUE);

      std::stringstream ss(tnsStr);
      Vector q(4);
      ss>>q(1)>>q(2)>>q(3)>>q(0);
      epoch.rotary = Rotary3d(q);

      if(arc.size() && (epoch.time <= arc.at(arc.size()-1).time))
        {logWarning<<"(epoch.time <= arc.at(arc.size()-1).time)"<<Log::endl; continue;}
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
