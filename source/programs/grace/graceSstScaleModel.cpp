/***********************************************/
/**
* @file graceSstScaleModel.cpp
*
* @brief First order Ensemble Averaging
*
* @author Saniya Behzadpour
* @date 2018-05-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This programs estimate satellite-to-satellite-tracking (SST) deterministic signals
due to eclipse transits and low-SNR values from post-fit residuals.
The low-SNR effects are estimated by directly using the residual values.
The ensemble averaging method is used to characterize the average properties of eclipse transit signal shapes across all transit events.
Each shape is assigned to one arc of 3 hours (default). This can be modefied by enabling \config{averagingInterval}.
)";

/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief First order Ensemble Averaging
* @ingroup programsGroup */
class GraceSstScaleModel
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceSstScaleModel, SINGLEPROCESS, "Ensamble Averaging of eclipse transit signals", Grace, Instrument)

/***********************************************/

void GraceSstScaleModel::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outNameEclipse, outNameSnr;
    FileName inName1, inName2, inNameSst, inNameSnr;
    UInt     timeMargin, waveLength, neighborNumber;
    UInt     arcCount = 0;
    UInt     countSst = 0;
    Bool     perArc = TRUE;
    Bool     lowSnr = FALSE;
    Bool     eclipseTransit = FALSE;

    readConfig(config, "inputfileGraceResiduals", inNameSst,  Config::MUSTSET, "",   "SST Residuals");
    readConfig(config, "timeMargin",              timeMargin, Config::MUSTSET, "25", "epochs before instrumental events");
    readConfig(config, "waveLength",              waveLength, Config::MUSTSET, "60", "length of the sample wave");
    if(readConfigSequence(config, "estimateEclipseTransitScale", Config::OPTIONAL, "", ""))
    {
      eclipseTransit = TRUE;
      readConfig(config, "outputfileScaleModel" ,        outNameEclipse, Config::MUSTSET, "", "");
      readConfig(config, "inputfileGrace1EclipseFactor", inName1,        Config::MUSTSET, "", "GRACE-A eclipse factors computed with integrated orbit");
      readConfig(config, "inputfileGrace2EclipseFactor", inName2,        Config::MUSTSET, "", "GRACE-B eclipse factors computed with integrated orbit");
      if(readConfigSequence(config, "averagingInterval", Config::OPTIONAL, "", ""))
      {
        perArc = FALSE;
        readConfig(config, "nearestNeighborNumber", neighborNumber, Config::DEFAULT, "24", "");
        endSequence(config);
      }
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateLowSnrScale", Config::OPTIONAL, "", ""))
    {
      lowSnr = TRUE;
      readConfig(config, "outputfileScaleModel", outNameSnr, Config::MUSTSET, "", "");
      readConfig(config, "inputfileGraceSstSNR", inNameSnr,  Config::MUSTSET, "", "GRACE SNR values");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // =======================

    std::vector<MiscValueArc>         arcListEF1, arcListEF2;
    std::vector<SatelliteTrackingArc> arcListRes;
    SatelliteTrackingArc              sstArc_1, sstArc_2;

    logStatus<<"read satellite data"<< lowSnr <<Log::endl;
    InstrumentFile fileSst(inNameSst);
    arcCount = fileSst.arcCount();
    arcListRes.resize(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        arcListRes.at(arcNo)  = fileSst.readArc(arcNo);

    if(eclipseTransit)
    {
      InstrumentFile fileEF1(inName1);
      InstrumentFile fileEF2(inName2);
      InstrumentFile::checkArcCount({fileSst, fileEF1, fileEF2});
      arcListEF1.resize(arcCount);
      arcListEF2.resize(arcCount);
      std::vector<std::vector<Double>> sumWaveP (arcCount); //sum of the waves in each arc with positive criteria
      std::vector<std::vector<Double>> sumWaveN (arcCount); //sum of the waves in each arc with negative criteria
      std::vector<UInt> nP (arcCount); //number of the waves in each arc with positive criteria
      std::vector<UInt> nN (arcCount); //number of the waves in each arc with positive criteria

      logStatus<<"find candidates and compute the sum"<<Log::endl;
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        arcListEF1.at(arcNo)  = fileEF1.readArc(arcNo);
        arcListEF2.at(arcNo)  = fileEF2.readArc(arcNo);
        sumWaveP.at(arcNo).resize(waveLength);
        sumWaveN.at(arcNo).resize(waveLength);
        nP.at(arcNo) = 1;
        nN.at(arcNo) = 1;
        Double criteria = 0;

        countSst = arcListRes.at(arcNo).size();
        for(UInt i=0; i<countSst; i++)
        {
          if(i < countSst - timeMargin)
          {
            criteria = arcListEF2.at(arcNo).at(i+timeMargin).value - arcListEF1.at(arcNo).at(i+timeMargin).value;
            if(criteria >0)
            {
              UInt maxj = waveLength;
              if (countSst - i < maxj)
                  maxj = countSst - i;
              for(UInt j=i; j< i+maxj; j++)
                sumWaveP.at(arcNo).at(j-i) += arcListRes.at(arcNo).at(j).rangeRate;
              i = i + maxj-1;
              nP.at(arcNo) = nP.at(arcNo) + 1 ;
            }
            else if(criteria < 0)
            {
              UInt maxj = waveLength;
              if (countSst - i < maxj)
                  maxj = countSst - i;
              for(UInt j=i; j< i+maxj; j++)
                sumWaveN.at(arcNo).at(j-i) += arcListRes.at(arcNo).at(j).rangeRate;
              i = i + maxj-1;
              nN.at(arcNo) = nN.at(arcNo) + 1 ;
            }
          }
        }
      }

      UInt m = perArc ? arcCount: (UInt) arcCount/neighborNumber; //number of data segments
      Matrix sumTotalP (waveLength,m);
      Matrix sumTotalN (waveLength,m);
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        for(UInt i=0; i< m; i++)
          if ((arcNo >= i*arcCount/m)&&(arcNo < (i+1)*arcCount/m))
          {
            for(UInt j=0; j< waveLength; j++)
              sumTotalP(j,i) += sumWaveP.at(arcNo).at(j) / nP.at(arcNo);
            for(UInt j=0; j< waveLength; j++)
              sumTotalN(j,i) += sumWaveN.at(arcNo).at(j) / nN.at(arcNo);
          }
      // Compute the desired signal
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        UInt k =0;
        for(UInt i=0; i< m; i++)
          if ((arcNo >= i*arcCount/m)&&(arcNo < (i+1)*arcCount/m))
              k= i;
        Double criteria = 0;
        countSst = arcListRes.at(arcNo).size();
        for(UInt i=0; i<countSst; i++)
        {
          if(i < (countSst - timeMargin))
          {
            criteria = arcListEF2.at(arcNo).at(i+timeMargin).value - arcListEF1.at(arcNo).at(i+timeMargin).value;
            if(criteria > 0)
            {
              UInt maxj = waveLength;
              if (countSst - i < maxj)
                  maxj = countSst - i;
              for(UInt j=i; j< i+maxj; j++)
              {
                SatelliteTrackingEpoch epoch;
                epoch.time  = arcListRes.at(arcNo).at(j).time;
                epoch.range = epoch.rangeRate = epoch.rangeAcceleration = sumTotalP(j-i,k)/(arcCount/m);
                sstArc_1.push_back(epoch);
              }
              i = i + maxj-1;
            }
            else if(criteria < 0)
            {
              UInt maxj = waveLength;
              if (countSst - i < maxj)
                  maxj = countSst - i;
              for(UInt j=i; j< i+maxj; j++)
              {
                SatelliteTrackingEpoch epoch;
                epoch.time  = arcListRes.at(arcNo).at(j).time;
                epoch.range = epoch.rangeRate = epoch.rangeAcceleration = sumTotalN(j-i,k)/(arcCount/m);
                sstArc_1.push_back(epoch);
              }
              i = i + maxj-1;
            }
            else
            {
              SatelliteTrackingEpoch epoch;
              epoch.time  = arcListRes.at(arcNo).at(i).time;
              epoch.range = epoch.rangeRate = epoch.rangeAcceleration = 0;
              sstArc_1.push_back(epoch);
            }
          }
          else
          {
            SatelliteTrackingEpoch epoch;
            epoch.time  = arcListRes.at(arcNo).at(i).time;
            epoch.range = epoch.rangeRate = epoch.rangeAcceleration = 0;
            sstArc_1.push_back(epoch);
          }
        }
      }
      logStatus<<"write the output eclipse transit signal <"<<outNameEclipse<<">"<<Log::endl;
      InstrumentFile::write(outNameEclipse, sstArc_1);
    } //if (eclipseTransit)

    if(lowSnr)
    {
      InstrumentFile fileSnr(inNameSnr);
      InstrumentFile::checkArcCount({fileSst, fileSnr});
      std::vector<MiscValuesArc> arcListSnr;
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        arcListSnr.push_back( fileSnr.readArc(arcNo) );
        countSst = arcListRes.at(arcNo).size();
        for(UInt i=0; i<countSst; i++)
        {
          SatelliteTrackingEpoch epoch;
          epoch.time  = arcListRes.at(arcNo).at(i).time;
          epoch.range = epoch.rangeRate = epoch.rangeAcceleration = 0;
          if (arcListSnr.at(arcNo).at(i).values(2) < 590)
               epoch.range = epoch.rangeRate = epoch.rangeAcceleration = arcListRes.at(arcNo).at(i).rangeRate;
          sstArc_2.push_back(epoch);
        }
      }
      logStatus<<"write the output low SNR signal <"<<outNameSnr<<">"<<Log::endl;
      InstrumentFile::write(outNameSnr, sstArc_2);
    } //if (lowSnr)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***************************************/

