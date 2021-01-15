/***********************************************/
/**
* @file gnssNanu2Matrix.cpp
*
* @brief Convert NANU (Notice Advisory to NAVSTAR users) files to outage Matrix files per GPS PRN and day.
*
* NANU file source: https://celestrak.com/GPS/NANU/
*
* @author Sebastian Strasser
* @date 2016-11-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert NANU (Notice Advisory to NAVSTAR users) files to outage Matrix files per GPS PRN and day.

NANU file source: \url{https://celestrak.com/GPS/NANU/}

Extracts satellite outage information from \config{inputfileNanu} and writes a matrix
file \configFile{outputfileMatrix}{matrix} for each PRN and each day. Each line represents an outage,
with start time in the first column and end time in the second column of the matrix.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "parser/dataVariables.h"
#include "inputOutput/file.h"
#include "files/fileMatrix.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Convert NANU (Notice Advisory to NAVSTAR users) files to outage Matrix files per GPS PRN and day.
* @ingroup programsConversionGroup */
class GnssNanu2Matrix
{
  class Outage
  {
  public:
    Time        timeStart;
    Time        timeEnd;
    std::string nanuType;
    UInt        nanuNumber;
    UInt        referenceNumber;
    UInt        prn;

    Bool operator< (const Outage& outage) const
    {
      return (this->nanuNumber < outage.nanuNumber);
    }

    Bool operator== (const Outage& outage) const
    {
      if(this->nanuType  != outage.nanuType)  return FALSE;
      if(this->prn       != outage.prn)       return FALSE;
      if(this->timeStart != outage.timeStart) return FALSE;
      if(this->timeEnd   != outage.timeEnd)   return FALSE;
      return TRUE;
    }
  };

  // implemented NANU types
  std::vector<std::string> nanuTypes = {"FCSTSUMM", "FCSTUUFN", "FCSTCANC", "UNUSUFN", "UNUSABLE", "UNUNOREF", "USABINIT", "DECOM"};

  Bool readNanuFile(const FileName &nanuFile, Outage &outage) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssNanu2Matrix, SINGLEPROCESS, "Convert NANU files to outage Matrix file per GPS PRN.", Conversion, Gnss)

/***********************************************/

void GnssNanu2Matrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameMatrix;
    std::vector<FileName> fileNameNanu;
    TimeSeriesPtr         timesIntervalPtr;
    Bool                  initiallyUnusable = TRUE;

    readConfig(config, "outputfileMatrix",  fileNameMatrix,    Config::MUSTSET,  "outage_{loopTime}.txt", "PRN is appended to file name");
    readConfig(config, "inputfileNanu",     fileNameNanu,      Config::MUSTSET,  "",                      "NANU file");
    readConfig(config, "intervals",         timesIntervalPtr,  Config::MUSTSET,  "",                      "for {loopTime} variable in outputfile");
    readConfig(config, "initiallyUnusable", initiallyUnusable, Config::DEFAULT,  "1",                     "set time periods before initial interval to unusable");
    if(isCreateSchema(config)) return;

    std::vector<Time> timesInterval;
    if(timesIntervalPtr)
      timesInterval = timesIntervalPtr->times();

    VariableList fileNameVariableList;
    addTimeVariables(fileNameVariableList);

    // read NANU files
    logStatus << "read NANU files" << Log::endl;
    const UInt nSat = 32;
    std::vector<std::vector<Outage>> outages(nSat);
    UInt fileCount = 0;
    for(const auto &fileName : fileNameNanu)
    {
      Outage outage;
      try
      {
        if(readNanuFile(fileName, outage))
        {
          outages.at(outage.prn-1).push_back(outage);
          fileCount++;
        }
      }
      catch(std::exception &e)
      {
        logWarning << "error reading file <" << fileName << ">: " << e.what() << Log::endl;
      }
    }
    logInfo << "  files read: " << fileCount << Log::endl;

    // clean up outages (delete forecasts if there is a summary, etc.)
    for(UInt iSat = 0; iSat < nSat; iSat++)
    {
      // sort outages by NANU number
      std::sort(outages.at(iSat).begin(), outages.at(iSat).end());

      // remove duplicate outages
      for(UInt i = 1; i < outages.at(iSat).size(); i++)
      {
        Outage outage = outages.at(iSat).at(i);
        Outage outagePrev = outages.at(iSat).at(i-1);

        if(outage == outagePrev)
        {
          // correct reference number of subsequent NANU messages
          for(UInt j = i; j < outages.at(iSat).size(); j++)
            if(outages.at(iSat).at(j).referenceNumber == outage.nanuNumber)
              outages.at(iSat).at(j).referenceNumber = outagePrev.nanuNumber;

          outages.at(iSat).erase(outages.at(iSat).begin()+i);
          i--;
        }
      }

      // clean up outages per satellite
      for(UInt i = 0; i < outages.at(iSat).size(); i++)
      {
        Outage outage = outages.at(iSat).at(i);

        // remove reference NANU
        if(outage.referenceNumber != NULLINDEX)
          for(UInt j = i; j --> 0 ;)
            if(outages.at(iSat).at(j).nanuNumber == outage.referenceNumber)
            {
              outages.at(iSat).erase(outages.at(iSat).begin()+j);
              i--;
              break;
            }

        // remove forecast outage if it was canceled
        if(outage.nanuType == "FCSTCANC")
        {
          // remove forecast outage NANU
          for(UInt j = i; j --> 0 ;)
            if(outages.at(iSat).at(j).nanuNumber == outage.referenceNumber)
            {
              outages.at(iSat).erase(outages.at(iSat).begin()+j);
              i--;
              break;
            }

          // remove forecast canceled NANU
          outages.at(iSat).erase(outages.at(iSat).begin()+i);
          i--;
          continue;
        }

        // remove NANU if start time >= end time
        if(outage.timeStart >= outage.timeEnd)
        {
          outages.at(iSat).erase(outages.at(iSat).begin()+i);
          i--;
          continue;
        }

        // remove 'initially usable' NANU messages and set timeEnd of preceding outage accordingly
        if(outage.nanuType == "USABINIT")
        {
          if(i > 0 && outages.at(iSat).at(i-1).timeEnd == date2time(2500,1,1))
            outages.at(iSat).at(i-1).timeEnd = outage.timeStart;
          else
          {
            logWarning << outage.prn%"G%02i: "s << "no 'until further notice' NANU found before NANU " << outage.nanuNumber << " (" << outage.nanuType << ")" << Log::endl;
            for(UInt j = 0; j < i; j++)
              logInfo << outages.at(iSat).at(j).prn%"  preceding NANU for G%02i: "s  << outages.at(iSat).at(j).nanuNumber << " " <<
                         outages.at(iSat).at(j).nanuType << " " << outages.at(iSat).at(j).prn << Log::endl;

            if(i == 0)
            {
              logInfo << outage.prn%"G%02i: set to unusable until NANU "s << outage.nanuNumber << " (" << outage.nanuType << ") start time " << outage.timeStart.dateTimeStr() << Log::endl;
              Outage outageNew;
              outageNew.prn       = outage.prn;
              outageNew.nanuType  = "UNUSABLE";
              outageNew.timeStart = timesInterval.at(0);
              outageNew.timeEnd   = outage.timeStart;
              outages.at(iSat).insert(outages.at(iSat).begin(), outageNew);
              i++;
            }
          }

          outages.at(iSat).erase(outages.at(iSat).begin()+i);
          i--;
          continue;
        }
      }

      if(initiallyUnusable)
      {
        Outage outageNew;
        outageNew.nanuType  = "UNUSABLE";
        outageNew.timeStart = Time();
        outageNew.timeEnd   = timesInterval.at(0);
        outages.at(iSat).insert(outages.at(iSat).begin(), outageNew);
      }


      // check for time period overlaps
      for(UInt i = 1; i < outages.at(iSat).size(); i++)
      {
        Outage outage = outages.at(iSat).at(i);
        Outage outagePrev = outages.at(iSat).at(i-1);

        if(outagePrev.timeEnd > outage.timeStart)
          logWarning << (iSat+1)%"G%02i"s << " NANU overlap: " << outagePrev.nanuNumber << " (" << outagePrev.nanuType << ") timeEnd " << outagePrev.timeEnd.dateTimeStr() << " > "
                     << outage.nanuNumber << " (" << outage.nanuType << ") timeStart: " << outage.timeStart.dateTimeStr() << " ==> manual investigation required" << Log::endl;
      }
    }

    // write outages to matrix file per satellite and day
    logStatus << "write outage matrix files" << Log::endl;
    Single::forEach(nSat, [&](UInt iSat)
    {
      // generate map of daily outage periods
      std::map<Int,std::vector<Time>> timeStart, timeEnd;
      for(UInt i = 0; i < outages.at(iSat).size(); i++)
      {
        const Int mjdStart = std::max(outages.at(iSat).at(i).timeStart.mjdInt(), timesInterval.at(0).mjdInt());
        const Int mjdEnd   = std::min(static_cast<Int>(std::ceil(outages.at(iSat).at(i).timeEnd.mjd())), timesInterval.back().mjdInt());
        for(Int mjdInt = mjdStart; mjdInt < mjdEnd; mjdInt++)
        {
          timeStart[mjdInt].push_back(std::max(Time(mjdInt,   0), outages.at(iSat).at(i).timeStart));
          timeEnd  [mjdInt].push_back(std::min(Time(mjdInt+1, 0), outages.at(iSat).at(i).timeEnd));
        }
      }

      // write daily outage periods to files
      for(const auto &kv : timeStart)
      {
        const Int  mjdInt = kv.first;
        const UInt count  = kv.second.size();
        Matrix A(count, 2);
        for(UInt i = 0; i < count; i++)
        {
          A(i,0) = timeStart[mjdInt].at(i).mjd();
          A(i,1) = timeEnd  [mjdInt].at(i).mjd();
        }

        UInt idInterval;
        for(idInterval = 0; idInterval < timesInterval.size(); idInterval++)
          if(timesInterval.at(idInterval) >= Time(mjdInt, 0))
            break;

        evaluateTimeVariables(idInterval, timesInterval.at(idInterval), timesInterval.at(idInterval+1), fileNameVariableList);
        writeFileMatrix(fileNameMatrix(fileNameVariableList).appendBaseName((iSat+1)%".G%02i"s), A);
      }
    });
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssNanu2Matrix::readNanuFile(const FileName &nanuFile, Outage &outage) const
{
  try
  {
    InFile file(nanuFile);

    UInt iLine = 0;
    std::string line;
    while(std::getline(file, line))
    {
      iLine++;

      // skip header lines
      if(iLine <= 2)
        continue;

      // NANU type
      if(iLine == 3)
      {
        if(line.find("NANU TYPE:") == std::string::npos)
          return FALSE;

        // check if NANU type is implemented
        outage.nanuType = String::trim(line.substr(18));
        if(std::find(nanuTypes.begin(), nanuTypes.end(), outage.nanuType) == nanuTypes.end())
          return FALSE;
      }

      // NANU number
      if(iLine == 4)
        outage.nanuNumber = std::stoi(line.substr(20,7));

      // reference NANU number
      if(iLine == 6)
      {
        if(line.substr(23,3) == "N/A")
          outage.referenceNumber = NULLINDEX;
        else
          outage.referenceNumber = std::stoi(line.substr(23,7));
      }

      // read PRN
      if(iLine == 9)
        outage.prn = std::stoi(line.substr(12,2));

      // start day
      if(iLine == 10)
      {
        UInt idxStart = outage.nanuType == "DECOM" ? 28 : 19;
        outage.timeStart += seconds2time(86400)*(std::stoi(line.substr(idxStart,3))-1);;
      }

      // start time
      if(iLine == 11)
      {
        UInt idxStart = outage.nanuType == "DECOM" ? 33 : 24;
        outage.timeStart += seconds2time(std::stoi(line.substr(idxStart,2))*3600 + std::stoi(line.substr(idxStart+2,2))*60);
      }

      // start year
      if(iLine == 12)
      {
        UInt idxStart = outage.nanuType == "DECOM" ? 44 : 35;
        outage.timeStart += date2time(std::stoi(line.substr(idxStart,4)), 1, 1);
      }

      // end day
      if(iLine == 13)
      {
        if(outage.nanuType == "DECOM" || outage.nanuType == "USABINIT" || outage.nanuType == "FCSTCANC" || line.substr(18,3) == "UFN" || line.substr(18,3) == "N/A")
          outage.timeEnd = date2time(2500, 1, 1);
        else
          outage.timeEnd = seconds2time(86400)*(std::stoi(line.substr(18,3))-1);
      }

      // end time
      if(iLine == 14)
        if(outage.timeEnd < date2time(2500, 1, 1))
          outage.timeEnd += seconds2time(std::stoi(line.substr(23,2))*3600 + std::stoi(line.substr(25,2))*60);

      // end year
      if(iLine == 15)
        if(outage.timeEnd < date2time(2500, 1, 1))
          outage.timeEnd += date2time(std::stoi(line.substr(34,4)), 1, 1);
    }

    if(outage.timeStart > outage.timeEnd)
    {
      logWarning << "ignored NANU " << outage.nanuNumber << " because start time > end time: " << outage.timeStart.dateTimeStr() << " > " << outage.timeEnd.dateTimeStr() << Log::endl;
      return FALSE;
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
