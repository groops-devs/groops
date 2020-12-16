/***********************************************/
/**
* @file orbit2Sp3Format.cpp
*
* @brief Write orbits (position, velocity, covariance) to SP3 format.
*
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2012-12-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Writes orbits to \href{https://files.igs.org/pub/data/format/sp3d.pdf}{SP3 format}.

SP3 orbits are usually given in the terrestrial reference frame (TRF), so providing \configClass{earthRotation}{earthRotationType}
automatically rotates the orbits from the celestial reference frame (CRF) to the TRF.
Since SP3 orbits often use the center of Earth as a reference, a correction from center of mass to center
of Earth can be applied to the orbits by providing \configClass{gravityfield}{gravityfieldType} (e.g. ocean tides).

See also \program{Sp3Format2Orbit}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Write orbits (position, velocity, covariance) to SP3 format.
* @ingroup programsConversionGroup */
class Orbit2Sp3Format
{
public:
  void run(Config &config);

  class Satellite
  {
  public:
    FileName inNameOrbit, inNameClock, inNameCov;
    std::string identifier;
    OrbitArc        orbit;
    MiscValueArc    clock;
    Covariance3dArc cov;
    std::vector<Time> times;
    UInt idEpoch = 0;
    Double orbitAccuracy;
  };
};

GROOPS_REGISTER_PROGRAM(Orbit2Sp3Format, SINGLEPROCESS, "Write orbits (position, velocity, covariance) to SP3 format.", Conversion, Orbit, Covariance, Instrument)
GROOPS_RENAMED_PROGRAM(Orbit2Sp3, Orbit2Sp3Format, date2time(2020, 8, 4))

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, Orbit2Sp3Format::Satellite &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileOrbit",      var.inNameOrbit,    Config::MUSTSET,  "", "");
  readConfig(config, "inputfileClock",      var.inNameClock,    Config::OPTIONAL, "", "");
  readConfig(config, "inputfileCovariance", var.inNameCov,      Config::OPTIONAL, "", "");
  readConfig(config, "identifier",          var.identifier,     Config::MUSTSET,  "", "3 characters (e.g. GNSS PRN: G01)");
  readConfig(config, "orbitAccuracy",       var.orbitAccuracy,  Config::DEFAULT,  "0", "[m] used for accuracy codes in header (0 = unknown)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void Orbit2Sp3Format::run(Config &config)
{
  try
  {
    FileName         outName;
    Bool             writeVel, useSp3k;
    EarthRotationPtr earthRotation;
    GravityfieldPtr  gravityfield;
    std::string      textLine1;
    std::vector<std::string> commentLines;
    std::vector<Satellite>   satellites;

    readConfig(config, "outputfile",          outName,        Config::MUSTSET,  "", "");
    readConfig(config, "satellite",           satellites,     Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",       earthRotation,  Config::OPTIONAL, "", "rotate data into Earth-fixed frame");
    readConfig(config, "gravityfield",        gravityfield,   Config::DEFAULT,  "", "degree 1 fluid mantle for CM2CE correction (SP3 orbits should be in center of Earth)");
    readConfig(config, "comment",             commentLines,   Config::OPTIONAL, "", "comment lines (77 char max)");
    readConfig(config, "firstLine",           textLine1,      Config::DEFAULT,  "RAW   IGS14 FIT  TUG", "Text for first line e.g:  u+U  IGS14 KIN ITSG");
    readConfig(config, "writeVelocity",       writeVel,       Config::DEFAULT,  "0", "write velocity in addition to position");
    readConfig(config, "useSp3kFormat",       useSp3k,        Config::DEFAULT,  "0", "use the extended sp3k format");
    if(isCreateSchema(config)) return;

    std::set<Time> timesSet; // unique set of times covering all satellites

    auto iter = satellites.begin();
    while(iter != satellites.end())
    {
      try
      {
        logStatus<<"read orbit file <"<<iter->inNameOrbit<<">"<<Log::endl;
        iter->orbit = InstrumentFile::read(iter->inNameOrbit);
        iter->clock = InstrumentFile::read(iter->inNameClock);
        iter->cov   = InstrumentFile::read(iter->inNameCov);
      }
      catch(std::exception &e)
      {
        logWarning << e.what() << " continue..." << Log::endl;
        iter = satellites.erase(iter);
        continue;
      }
      iter->times = iter->orbit.times();
      Arc::checkSynchronized({iter->orbit, iter->clock, iter->cov});
      std::copy(iter->times.begin(), iter->times.end(), std::inserter(timesSet, timesSet.end()));
      iter++;
    }
    const std::vector<Time> times(timesSet.begin(), timesSet.end());
    const Double sampling = medianSampling(times).seconds();

    logStatus<<"write file <"<<outName<<">"<<Log::endl;
    OutFile file(outName);

    const Time timeStart = times.at(0);
    UInt year, month, day, hour, minute;
    Double second;
    timeStart.date(year, month, day, hour, minute, second);

    // first line
    file<<"#"<<(useSp3k?"k":"d")<<(writeVel?"V":"P")<<timeStart%"%y %m %d %H %M %011.8S "s<<times.size()%"%7i "s<<textLine1<<std::endl;

    // second line
    file<<"## "<<timeStart%"%4W"<<" "<<((timeStart.mjdInt()-44244)%7*86400+timeStart.mjdMod()*86400)%"%15.8f "s<<sampling%"%14.8f "s<<timeStart.mjdInt()%"%5i "s<<timeStart.mjdMod()%"%15.13f"s<<std::endl;

    // satellite identifier lines
    const UInt satPerLine = 17;
    const UInt lineCount = std::max<UInt>(5, (satellites.size()+satPerLine+1)/satPerLine);
    for(UInt i = 0; i < lineCount; i++)
    {
      if(i == 0)
        file<<"+  "<<satellites.size()%"% 3i"s<<"   ";
      else
        file<<"+        ";
      for(UInt j = 0; j < satPerLine; j++)
      {
        const UInt idSat = i*satPerLine+j;
        if(idSat < satellites.size())
          file<<satellites.at(i*satPerLine+j).identifier;
        else
          file<<"  0";
      }
      file<<std::endl;
    }

    // satellite accuracy lines
    for(UInt i = 0; i < lineCount; i++)
    {
      file<<"++       ";
      for(UInt j = 0; j < satPerLine; j++)
      {
        const UInt idSat = i*satPerLine+j;
        if(idSat < satellites.size() && satellites.at(i*satPerLine+j).orbitAccuracy > 0)
          file<<std::round(std::log2(satellites.at(i*satPerLine+j).orbitAccuracy*1e3))%"% 3i"s;
        else
          file<<"  0";
      }
      file<<std::endl;
    }

    file<<"%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc"<<std::endl;
    file<<"%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc"<<std::endl;
    file<<"%f  1.2500000  1.025000000  0.00000000000  0.000000000000000"<<std::endl;
    file<<"%f  0.0000000  0.000000000  0.00000000000  0.000000000000000"<<std::endl;
    file<<"%i    0    0    0    0      0      0      0      0         0"<<std::endl;
    file<<"%i    0    0    0    0      0      0      0      0         0"<<std::endl;

    // comment lines
    for(UInt i = 0; i < std::max<UInt>(commentLines.size(), 4); i++)
    {
      if(i < commentLines.size())
        file<<"/* "<<commentLines.at(i)<<std::endl;
      else
        file<<"/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"<<std::endl;
    }

    logTimerStart;
    for(UInt i=0; i<times.size(); i++)
    {
      logTimerLoop(i, times.size());

      file<<"*  "<<times.at(i)%"%y %m %d %H %M %011.8S"s<<std::endl;

      Rotary3d rot;
      Vector3d omega;
      if(earthRotation)
      {
        rot   = earthRotation->rotaryMatrix(times.at(i));
        if(writeVel)
          omega = earthRotation->rotaryAxis(times.at(i));
      }

      const SphericalHarmonics harmonics = gravityfield->sphericalHarmonics(times.at(i), 1, 1);
      const Vector coeff = harmonics.x(); // [c00, c10, c11, s11]
      const Vector3d cm2ceCorrection = std::sqrt(3.) * harmonics.R() * Vector3d(coeff(2), coeff(3), coeff(1));

      for(auto &&satellite : satellites)
      {
        if(satellite.idEpoch >= satellite.times.size() || satellite.times.at(satellite.idEpoch) != times.at(i))
          continue; // no data for this epoch

        const Vector3d position = rot.rotate(satellite.orbit.at(satellite.idEpoch).position) + cm2ceCorrection;
        const std::string clock = satellite.clock.size() ? (satellite.clock.at(satellite.idEpoch).value*1e6)%"%14.6f"s : " 999999.999999";
        if(useSp3k)
          file<<"P"<<satellite.identifier<<(position.x()/1000)%"%14.7f"s<<(position.y()/1000)%"%14.7f"s<<(position.z()/1000)%"%14.7f"s<<clock<<"             "<<std::endl;
        else
          file<<"P"<<satellite.identifier<<(position.x()/1000)%"%14.6f"s<<(position.y()/1000)%"%14.6f"s<<(position.z()/1000)%"%14.6f"s<<clock<<"             "<<std::endl;

        if(writeVel)
        {
          const Vector3d velocity = rot.rotate(satellite.orbit.at(satellite.idEpoch).velocity-crossProduct(omega, satellite.orbit.at(satellite.idEpoch).position));
          if(useSp3k)
            file<<"V"<<satellite.identifier<<(velocity.x()*10)%"%14.7f"s<<(velocity.y()*10)%"%14.7f"s<<(velocity.z()*10)%"%14.7f 999999.999999             "s<<std::endl;
          else
            file<<"V"<<satellite.identifier<<(velocity.x()*10)%"%14.6f"s<<(velocity.y()*10)%"%14.6f"s<<(velocity.z()*10)%"%14.6f 999999.999999             "s<<std::endl;
        }

        if(satellite.cov.size())
        {
          Tensor3d cv = earthRotation ? rot.rotate(satellite.cov.at(satellite.idEpoch).covariance) : satellite.cov.at(satellite.idEpoch).covariance;
          if(useSp3k)
          {
            file<<"EPx "
                <<std::min(1000*std::sqrt(cv.xx()), 99.)%"%4.1f "s
                <<std::min(1000*std::sqrt(cv.yy()), 99.)%"%4.1f "s
                <<std::min(1000*std::sqrt(cv.zz()), 99.)%"%4.1f "s;
          }
          else
          {
            file<<"EP  "
                <<std::min(std::round(1000*std::sqrt(cv.xx())), 9999.)%"%4i "s
                <<std::min(std::round(1000*std::sqrt(cv.yy())), 9999.)%"%4i "s
                <<std::min(std::round(1000*std::sqrt(cv.zz())), 9999.)%"%4i "s;
          }
          file<<"        " // clk-sdev
              <<std::min(std::round(10000000*cv.xy()/(std::sqrt(cv.xx())*std::sqrt(cv.yy()))), 99999999.)%"%8.0f "s
              <<std::min(std::round(10000000*cv.xz()/(std::sqrt(cv.xx())*std::sqrt(cv.zz()))), 99999999.)%"%8.0f "s
              <<"         " // xc
              <<std::min(std::round(10000000*cv.yz()/(std::sqrt(cv.yy())*std::sqrt(cv.zz()))), 99999999.)%"%8.0f "s
              <<"         " // yc
              <<"         " // zc
              <<std::endl;
        }

        satellite.idEpoch++;
      }
    }
    logTimerLoopEnd(times.size());

    file<<"EOF"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
