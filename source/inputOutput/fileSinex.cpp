/***********************************************/
/**
* @file fileSinex.cpp
*
* @brief SINEX file representation.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-05-15
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "inputOutput/logging.h"
#include "inputOutput/file.h"
#include "inputOutput/system.h"
#include "config/config.h"
#include "fileSinex.h"

/***********************************************/

SinexBlockPtr Sinex::addBlock(const std::string &label)
{
  try
  {
    SinexBlockPtr block = std::make_shared<SinexBlock>();
    block->label = label;
    blocks.push_back(block);
    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SinexBlockPtr Sinex::findBlock(const std::string &label)
{
  try
  {
    auto iter = std::find_if(blocks.begin(), blocks.end(), [&](const auto &b) {return b->label == label;});
    if(iter == blocks.end())
      throw(Exception("SINEX block not found: "+label));
    return *iter;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::format(Double value, UInt length, UInt precision)
{
  try
  {
    std::string s = value%("%"+length%"%i"s+"."+precision%"%i"s+"f");
    if(s.size() > length && s.substr(0,3) == "-0.")
      return "-." + s.substr(3, s.size());
    return s;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::time2str(Time time, Bool fourDigitYear)
{
  try
  {
    if(time == Time() || time >= date2time(2500, 1, 1))
      return (fourDigitYear ? "0000" : "00") + ":000:00000"s;

    // round to full second including rollover
    UInt   year, month, day, hour, minute;
    Double second;
    time.date(year, month, day, hour, minute, second);
    time = date2time(year, month, day, hour, minute, std::round(second)+0.1);
    time.date(year, month, day, hour, minute, second);

    std::stringstream ss;
    if(fourDigitYear)
      ss<<year%"%04i"s<<":";
    else
      ss<<(year%100)%"%02i"s<<":";
    ss<<time.dayOfYear()%"%03i"s<<":";
    ss<<std::round(time.mjdMod()*86400)%"%05i"s;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Time Sinex::str2time(const std::string &line, std::size_t pos, Bool zeroIsMaxTime, Bool fourDigitYear)
{
  try
  {
    UInt posOffset = fourDigitYear ? 2 : 0;
    UInt year = static_cast<UInt>(String::toInt(line.substr(pos+0, 2+posOffset)));
    UInt day  = static_cast<UInt>(String::toInt(line.substr(pos+3+posOffset, 3)));
    UInt sec  = static_cast<UInt>(String::toInt(line.substr(pos+7+posOffset, 5)));
    if((year == 0) && (day == 0) && (sec == 0))
      return zeroIsMaxTime ? date2time(2500, 1, 1) : Time();
    if(!fourDigitYear)
      year += (year <= 50) ? 2000 : 1900;
    return date2time(year,1,1) + mjd2time(day-1.) + seconds2time(static_cast<Double>(sec));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileSinex(const FileName &fileName, const Sinex &sinex)
{
  try
  {
    OutFile file(fileName);
    file<<sinex.header<<std::endl;
    file<<"*"<<std::string(79, '-')<<std::endl;
    for(const auto &block : sinex.blocks)
    {
      file<<"+"<<block->label<<std::endl;
      file<<block->ss.str();
      file<<"-"<<block->label<<std::endl;
      file<<"*"<<std::string(79, '-')<<std::endl;
    }
    file<<"%ENDSNX";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileSinex(const FileName &fileName, Sinex &sinex)
{
  try
  {
    InFile file(fileName);

    std::string line;
    if(file.peek() == '%')
      std::getline(file, sinex.header);
    else
      logWarning<<"mandatory header in first line (starting with %) is missing"<<Log::endl;

    SinexBlockPtr block;
    while(std::getline(file, line))
    {
      line = String::trimRight(line); // trim from end
      if(line.empty() || (line.at(0) == '*')) // skip comments
        continue;
      else if(line.at(0) == '%')              // %ENDSNX
        break;
      else if(line.at(0) == '+')              // start data block
      {
        if(block && (block->label == "FILE/COMMENT"))
          continue;
        if(block)
          throw(Exception("New SINEX block starts unexpectedly: '"+line+"' in block '"+block->label+"'"));
        block = sinex.addBlock(String::trim(line.substr(1)));
      }
      else if(line.at(0) == '-') // end data block
      {
        // Do not close the FILE/COMMENT block in case of comment lines starting incorrectly with "-"
        if(block && (block->label == "FILE/COMMENT") && (block->label != String::trim(line.substr(1))))
          continue;
        if(!block || (block->label != String::trim(line.substr(1))))
          throw(Exception("SINEX block ends unexpectedly: '"+line+"'"));
        block = nullptr;
      }
      else if(!block || (line.at(0) != ' ')) // unknown line
      {
        if(!block || (block->label != "FILE/COMMENT"))
          logWarning<<"Unknown line identifier: '"<<line<<"'"<<Log::endl;
      }
      else
        block->lines.push_back(String::trimRight(line));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, Sinex &sinex, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    Time                     timeStart, timeEnd;
    std::string              agencyCode, observationCode, constraintCode, solutionContent;
    std::string              description, output, contact, software, hardware, input;
    std::vector<std::string> comments;
    FileName                 fileNameComment;

    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;
    readConfig(config, "agencyCode",       agencyCode,       Config::OPTIONAL, "TUG",    "identify the agency providing the data");
    readConfig(config, "timeStart",        timeStart,        Config::OPTIONAL, "",       "start time of the data");
    readConfig(config, "timeEnd",          timeEnd,          Config::OPTIONAL, "",       "end time of the data ");
    readConfig(config, "observationCode",  observationCode,  Config::OPTIONAL, "C",      "technique used to generate the SINEX solution");
    readConfig(config, "constraintCode",   constraintCode,   Config::OPTIONAL, "2",      "0: tight constraint, 1: siginficant constraint, 2: unconstrained");
    readConfig(config, "solutionContent",  solutionContent,  Config::OPTIONAL, "",       "solution types contained in the SINEX solution (S O E T C A)");
    readConfig(config, "description",      description,      Config::OPTIONAL, "",       "organizitions gathering/alerting the file contents");
    readConfig(config, "contact",          contact,          Config::OPTIONAL, "",       "Address of the relevant contact. e-mail");
    readConfig(config, "output",           output,           Config::OPTIONAL, "",       "Description of the file contents");
    readConfig(config, "input",            input,            Config::OPTIONAL, "",       "Brief description of the input used to generate this solution");
    readConfig(config, "software",         software,         Config::OPTIONAL, "GROOPS", "Software used to generate the file");
    readConfig(config, "hardware",         hardware,         Config::OPTIONAL, "",       "Computer hardware on which above software was run");
    readConfig(config, "inputfileComment", fileNameComment,  Config::OPTIONAL, "",       "comments in the comment block from a file (truncated at 80 characters)");
    readConfig(config, "comment",          comments,         Config::OPTIONAL, "",       "comments in the comment block");
    endSequence(config);
    if(isCreateSchema(config))
      return TRUE;

    // header line
    std::stringstream ss;
    ss<<"%=SNX 2.02 "<<Sinex::resize(agencyCode.substr(0,3),3)<<" "<<Sinex::time2str(System::now())<<" "<<Sinex::resize(agencyCode.substr(0,3),3);
    ss<<" "<<Sinex::time2str(timeStart)<<" "<<Sinex::time2str(timeEnd)<<" "<<Sinex::resize(observationCode.substr(0,1),1)<<" 00000";
    ss<<" "<<Sinex::resize(constraintCode.substr(0,1),1)<<" "<<solutionContent;
    sinex.header = ss.str();

    {
      SinexBlockPtr block = sinex.addBlock("FILE/REFERENCE");
      *block<<"*INFO_TYPE_________ INFO________________________________________________________"<<std::endl;
      if(!description.empty()) *block<<" DESCRIPTION        "<<description<<std::endl;
      if(!output.empty())      *block<<" OUTPUT             "<<output<<std::endl;
      if(!contact.empty())     *block<<" CONTACT            "<<contact<<std::endl;
      if(!software.empty())    *block<<" SOFTWARE           "<<software<<std::endl;
      if(!hardware.empty())    *block<<" HARDWARE           "<<hardware<<std::endl;
      if(!input.empty())       *block<<" INPUT              "<<input<<std::endl;
    }

    // comment block
    if(!fileNameComment.empty() || comments.size())
    {
      SinexBlockPtr block = sinex.addBlock("FILE/COMMENT");
      if(!fileNameComment.empty())
      {
        InFile commentFile(fileNameComment);
        std::string line;
        while(std::getline(commentFile, line))
          *block<<" "<<line<<std::endl;
      }
      for(const auto &line : comments)
        *block<<" "<<line<<std::endl;
    }

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
