/***********************************************/
/**
* @file system.cpp
*
* @brief Operating system related functions.
*
* @author Andreas Kvas
* @date 2017-03-21
*
*/
/***********************************************/

#ifdef _WIN32
  #define popen  _popen
  #define pclose _pclose
#endif

#include <ctime>
#include <filesystem>
#include "base/importStd.h"
#include "base/string.h"
#include "base/time.h"
#include "system.h"

/***********************************************/

Bool System::exec(const std::string &command, std::vector<std::string> &output)
{
  try
  {
    std::FILE *pipe = popen(command.c_str(), "r");
    if(!pipe)
      throw(Exception("Cannot open pipe"));

    char buffer[1024];
    while(std::fgets(buffer, sizeof(buffer), pipe) != nullptr)
    {
      std::string line(buffer);
      output.push_back(line.substr(0, line.find_last_of('\n')));
    }

    return (pclose(pipe) == 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool System::exec(const std::string &command)
{
  std::vector<std::string> output;
  return exec(command, output);
}

/***********************************************/
/***********************************************/

Bool System::createDirectories(const FileName &fileName)
{
  if(isDirectory(fileName))
    return TRUE;
  return std::filesystem::create_directories(fileName.str());
}

/***********************************************/

Bool System::remove(const FileName &fileName)
{
  return (std::filesystem::remove_all(fileName.str()) != 0);
}

/***********************************************/

Bool System::move(const FileName &fileNameOld, const FileName &fileNameNew)
{
  try
  {
    if(!exists(fileNameOld))
      return FALSE;
    if(isDirectory(fileNameNew))
      std::filesystem::rename(fileNameOld.str(), fileNameNew.append(fileNameOld.stripDirectory()).str());
    else
      std::filesystem::rename(fileNameOld.str(), fileNameNew.str());
    return TRUE;
  }
  catch(std::exception &/*e*/)
  {
    return FALSE;
  }
}

/***********************************************/

Bool System::exists(const FileName &fileName)
{
  try
  {
    return fileList(fileName).size();
  }
  catch(std::exception &/*e*/)
  {
    return FALSE;
  }
}

/***********************************************/

std::vector<FileName> System::fileList(const FileName &fileName)
{
  try
  {
    // contains no wildcards?
    if(fileName.str().find_first_of("*?") == std::string::npos)
      return (std::filesystem::exists(fileName.str()) ? std::vector<FileName>{fileName} : std::vector<FileName>{});

    // split path into parts
    std::vector<std::string> parts = String::split(fileName.str(), "/\\");

    // deepest directory without wildcards
    std::string path;
    UInt level = 0;
    while((level<parts.size()-1) && (parts.at(level).find_first_of("*?") == std::string::npos))
      path += parts.at(level++) + "/";
    if(path.empty())
      path = currentWorkingDirectory().str() + "/";
    if(!std::filesystem::exists(path))
      return std::vector<FileName>{};

    // search directory for matching entries
    std::vector<FileName> list;
    const std::regex pattern = String::wildcard2regex(parts.at(level));
    for(auto const &entry : std::filesystem::directory_iterator{path})
      if(std::regex_match(entry.path().filename().string(), pattern))
      {
        if(level == parts.size()-1)
          list.push_back(path+entry.path().filename().string());
        else if(std::filesystem::is_directory(entry))
        {
          // replace pattern by current entry
          std::string path2 = path + entry.path().filename().string() + "/";
          for(UInt i=level+1; i<parts.size()-1; i++)
            path2 += parts.at(i) + "/";
          std::vector<FileName> list2 = fileList(path2+parts.back());
          list.insert(list.end(), list2.begin(), list2.end());
        }
      }

    return list;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt System::fileSize(const FileName &fileName)
{
  try
  {
    return std::filesystem::file_size(fileName.str());
  }
  catch(std::exception &/*e*/)
  {
    return 0;
  }
}

/***********************************************/

Bool System::isDirectory(const FileName &fileName)
{
  return std::filesystem::is_directory(fileName.str());
}

/***********************************************/

FileName System::currentWorkingDirectory()
{
  return FileName(std::filesystem::current_path().string());
}

/***********************************************/

std::vector<FileName> System::directoryListing(const FileName &path, const std::regex &pattern)
{
  try
  {
    std::vector<FileName> files;
    for(auto const &entry : std::filesystem::directory_iterator{path.str()})
      if(std::filesystem::is_regular_file(entry) && std::regex_match(entry.path().filename().string(), pattern))
        files.push_back(FileName(entry.path().filename().string()));
    return files;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Time System::now()
{
  std::time_t tt = std::time(nullptr);
  std::tm     t  = *std::localtime(&tt);
  return date2time(t.tm_year+1900, t.tm_mon+1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec);
}

/***********************************************/
