/***********************************************/
/**
* @file goceXmlEggNom1b.cpp
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
class GoceXmlEggNom1b
{
  void readGradiometer   (XmlNodePtr &dataNode, GradiometerArc    &arc);
  void readAccelerometer (XmlNodePtr &dataNode, AccelerometerArc  &arc);
  void readStarCamera    (XmlNodePtr &dataNode, StarCameraArc     &arc);
  void readAngularRate   (XmlNodePtr &dataNode, AccelerometerArc  &arc);
  void readAngularAcc    (XmlNodePtr &dataNode, AccelerometerArc  &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GoceXmlEggNom1b, SINGLEPROCESS, "Read ESA XML GOCE Data", Conversion, Instrument)

/***********************************************/

void GoceXmlEggNom1b::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName sggName, accName, scaName, angRateName, angAccName;
    std::vector<FileName> fileNames;

    readConfig(config, "outputfileGradiometer",    sggName,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAccelerometer",  accName,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileStarCamera",     scaName,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAngularRate",    angRateName,   Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAngularAcc",     angAccName,    Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",                fileNames,     Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    GradiometerArc     gradiometerArc;
    AccelerometerArc   accelerometerArc;
    StarCameraArc      starCameraArc;
    AccelerometerArc   angRateArc;
    AccelerometerArc   angAccArc;

    for(auto &fileName : fileNames)
    {
      logStatus<<"read file <"<<fileName<<">"<<Log::endl;
      InFile     file(fileName);
      XmlNodePtr rootNode = XmlNode::read(file);
      XmlNodePtr dataNode = getChild(rootNode, "Data_Block", TRUE);

      if(!sggName.empty())     readGradiometer   (dataNode, gradiometerArc);
      if(!accName.empty())     readAccelerometer (dataNode, accelerometerArc);
      if(!scaName.empty())     readStarCamera    (dataNode, starCameraArc);
      if(!angRateName.empty()) readAngularRate   (dataNode, angRateArc);
      if(!angAccName.empty())  readAngularAcc    (dataNode, angAccArc);
    }

    // Daten speichern
    // ---------------
    if(!sggName.empty())
    {
      logStatus<<"write gradiometer file <"<<sggName<<">"<<Log::endl;
      InstrumentFile::write(sggName, gradiometerArc);
      Arc::printStatistics(gradiometerArc);
    }

    if(!accName.empty())
    {
      logStatus<<"write accelerometer file <"<<accName<<">"<<Log::endl;
      InstrumentFile::write(accName, accelerometerArc);
      Arc::printStatistics(accelerometerArc);
    }

    if(!scaName.empty())
    {
      logStatus<<"write starCamera file <"<<scaName<<">"<<Log::endl;
      InstrumentFile::write(scaName, starCameraArc);
      Arc::printStatistics(starCameraArc);
    }

    if(!angRateName.empty())
    {
      logStatus<<"write angular rate file <"<<angRateName<<">"<<Log::endl;
      InstrumentFile::write(angRateName, angRateArc);
      Arc::printStatistics(angRateArc);
    }

    if(!angRateName.empty())
    {
      logStatus<<"write angular acceleration file <"<<angAccName<<">"<<Log::endl;
      InstrumentFile::write(angAccName, angAccArc);
      Arc::printStatistics(angAccArc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GoceXmlEggNom1b::readGradiometer(XmlNodePtr &dataNode, GradiometerArc &arc)
{
  try
  {
    GradiometerEpoch epoch;

    XmlNodePtr node  = getChild(dataNode, "EGG_GGT_DS", TRUE);

    UInt epochCount = childCount(node, "EGG_GGT_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(node, "EGG_GGT_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr node2 = getChild(epochNode, "Gravity_Grad_Tensor", TRUE);
      std::string tnsStr;
      readXml(node2, "U_G", tnsStr, TRUE);

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
    logWarning<<"In GoceXmlEggNom1b::readGradiometer:\n"<<e.what()<<" continue..."<<Log::endl;
  }
}

/***********************************************/

void GoceXmlEggNom1b::readAccelerometer(XmlNodePtr &dataNode, AccelerometerArc &arc)
{
  try
  {
    AccelerometerEpoch epoch;

    XmlNodePtr node  = getChild(dataNode, "EGG_CCD_DS", TRUE);

    UInt epochCount = childCount(node, "EGG_CCD_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(node, "EGG_CCD_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr node2 = getChild(epochNode, "Acc_Ccm", TRUE);
      std::string xStr, yStr, zStr;
      readXml(node2, "X", xStr, TRUE);
      readXml(node2, "Y", yStr, TRUE);
      readXml(node2, "Z", zStr, TRUE);

      Double tmp;
      std::stringstream ssx(xStr);
      std::stringstream ssy(yStr);
      std::stringstream ssz(zStr);
      ssx>>epoch.acceleration.x()>>tmp>>tmp;
      ssy>>tmp>>epoch.acceleration.y()>>tmp;
      ssz>>tmp>>tmp>>epoch.acceleration.z();

      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    logWarning<<"In GoceXmlEggNom1b::readAccelerometer:\n"<<e.what()<<" continue..."<<Log::endl;
  }
}

/***********************************************/

void GoceXmlEggNom1b::readStarCamera(XmlNodePtr &dataNode, StarCameraArc &arc)
{
  try
  {
    StarCameraEpoch epoch;

    XmlNodePtr node       = getChild(dataNode, "EGG_IAQ_DS", TRUE);
    UInt       epochCount = childCount(node, "EGG_IAQ_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(node, "EGG_IAQ_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr node2 = getChild(epochNode, "Corr_Quat", TRUE);
      std::string str;
      readXml(node2, "Q_Grad", str, TRUE);

      Vector q(4);
      std::stringstream ss(str);
      ss>>q(1)>>q(2)>>q(3)>>q(0);
      epoch.rotary = Rotary3d(q);

      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    logWarning<<"In GoceXmlEggNom1b::readStarCamera:\n"<<e.what()<<" continue..."<<Log::endl;
  }
}

/***********************************************/

void GoceXmlEggNom1b::readAngularRate(XmlNodePtr &dataNode,  AccelerometerArc &arc)
{
  try
  {
    AccelerometerEpoch epoch;

    XmlNodePtr node       = getChild(dataNode, "EGG_GAR_DS", TRUE);
    UInt       epochCount = childCount(node, "EGG_GAR_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(node, "EGG_GAR_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr node2 = getChild(epochNode, "Corr_Est_Ang_Rate", TRUE);
      std::string text;
      readXml(node2, "Iar_Est", text, TRUE);

      std::stringstream ss(text);
      ss>>epoch.acceleration.x();
      ss>>epoch.acceleration.y();
      ss>>epoch.acceleration.z();

      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    logWarning<<"In GoceXmlEggNom1b::readAngularRate:\n"<<e.what()<<" continue..."<<Log::endl;
  }
}

/***********************************************/

void GoceXmlEggNom1b::readAngularAcc(XmlNodePtr &dataNode,  AccelerometerArc &arc)
{
  try
  {
    AccelerometerEpoch epoch;

    XmlNodePtr node       = getChild(dataNode, "EGG_CGA_DS", TRUE);
    UInt       epochCount = childCount(node, "EGG_CGA_1i", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(node, "EGG_CGA_1i", TRUE);
      Double t = 0;
      readXml(epochNode, "Tt_GPS", t, TRUE);
      epoch.time = seconds2time(t) + date2time(1980,1,6);

      XmlNodePtr node2 = getChild(epochNode, "Cal_Grad_Ang_Acc", TRUE);
      std::string text;
      readXml(node2, "CGA", text, TRUE);

      std::stringstream ss(text);
      ss>>epoch.acceleration.x();
      ss>>epoch.acceleration.y();
      ss>>epoch.acceleration.z();

      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    logWarning<<"In GoceXmlEggNom1b::readAngularAcc:\n"<<e.what()<<" continue..."<<Log::endl;
  }
}

/***********************************************/
/***********************************************/
