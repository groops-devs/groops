/***********************************************/
/**
* @file program.h
*
* @brief Interface for applications in groops.
*
* @author Torsten Mayer-Guerr
* @date 2008-07-27
*
*/
/***********************************************/

#ifndef __GROOPS_PROGRAM__
#define __GROOPS_PROGRAM__

#include "base/import.h"
#include "inputOutput/logging.h"
#include "config/config.h"
#include "parallel/parallel.h"

/***** DEFINE **********************************/

#ifndef GROOPS_TAGS
#define GROOPS_TAGS(tag_, ...) enum Tags {tag_, __VA_ARGS__};
#endif

namespace Program
{
  constexpr Bool SINGLEPROCESS = TRUE;
  constexpr Bool PARALLEL      = FALSE;

  /** @brief Tags to desribe programs. Used in @a GROOPS_REGISTER_PROGRAM. */
  GROOPS_TAGS(NONE,
              Covariance,
              DoodsonHarmonics,
              Gnss,
              Grace,
              Gravityfield,
              Grid,
              Instrument,
              KalmanFilter,
              Matrix,
              Misc,
              Noise,
              NormalEquation,
              Orbit,
              Plot,
              PotentialCoefficients,
              Preprocessing,
              Residuals,
              Simulation,
              Slr,
              SpatialTimeSeries,
              Statistics,
              System,
              TimeSeries,
              TimeSplines,
              VariationalEquation,
              /* others */
              Conversion, Provisional, Deprecated)

  extern const char *tagStrings[];
} // namespace Program

/***** DEFINE **********************************/

#ifndef DOCSTRING
#define DOCSTRING "" //"documentation is missing.\n"
// #warning ******** documentation is missing  **********
#endif

/***** DEFINE **********************************/

/** @brief Register a program in schema and documentation.
* A list of tags (from @a GROOPS_TAGS) followed the @a _name of the program and the short @a _description.
* The first tag should agree with the first part of the name/directory.
* The list should also include tags for the output/input file formats. */
#define GROOPS_REGISTER_PROGRAM(_name, _parallel, _description, ...) \
namespace Program\
{\
  class _Program##_name : public Program\
  {\
  public:\
    std::string       name()              const {return #_name;}\
    std::string       description()       const {return _description;}\
    std::string       documentation()     const {return DOCSTRING;}\
    std::vector<Tags> tags()              const {return std::vector<Tags>{__VA_ARGS__};}\
    Bool              isSingleProcess()   const {return _parallel;}\
    void              run(Config &config, Parallel::CommunicatorPtr comm) const\
    {\
      ::_name tmp;\
      if(!isSingleProcess())\
        tmp.run(config, comm);\
      else if(::Parallel::isMaster(comm))\
        tmp.run(config, ::Parallel::selfCommunicator());\
    }\
  };\
  static _Program##_name program##_name;\
}\
// end macro

/***** DEFINE **********************************/

/** @brief Register an old name of a renamed program.
* This macro can be used to keep old config files working. */
#define GROOPS_RENAMED_PROGRAM(_oldname, _newname, _time) \
namespace Program\
{\
  class _RenamedProgram##_oldname : public RenamedProgram\
  {\
  public:\
    _RenamedProgram##_oldname(const Renamed &renamed) : RenamedProgram(renamed) {}\
  };\
  static _RenamedProgram##_oldname renamedprogram##_oldname(RenamedProgram::Renamed(#_oldname, #_newname, _time));\
}\
// end macro

/***** CLASS ***********************************/

namespace Program
{

/** @brief Interface for applications in groops.
* @ingroup programsGroup
* An application must implement a @a run(Config &config) function. */
class Program
{
public:
  Program() {programList(this);}
  virtual ~Program() {}

  virtual std::string       name()              const = 0;
  virtual std::string       description()       const = 0;
  virtual std::string       documentation()     const = 0;
  virtual std::vector<Tags> tags()              const = 0;
  virtual Bool              isSingleProcess()   const = 0;
  virtual void              run(Config &config, Parallel::CommunicatorPtr comm) const = 0;

  static std::vector<Program*> programList(Program *program=nullptr);
  static void sortList(std::vector<Program*> &list);
};

/***** CLASS ***********************************/

/** @brief Interface for renamed applications.
* @ingroup programsGroup */
class RenamedProgram
{
public:
  class Renamed
  {
  public:
    std::string oldName;
    std::string newName;
    Time        time;
    Renamed(const std::string &oldName_, const std::string &newName_, const Time &time_) : oldName(oldName_), newName(newName_), time(time_) {}
  };

  RenamedProgram(const Renamed &renamed) {renamedList(renamed);}
  virtual ~RenamedProgram() {}

  static std::vector<Renamed> renamedList(const Renamed &renamed = Renamed("", "", Time()));
};

} // namespace Program

/***********************************************/

#endif /* __GROOPS__ */

