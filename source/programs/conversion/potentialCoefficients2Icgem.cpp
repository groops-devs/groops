/***********************************************/
/**
* @file potentialCoefficients2Icgem.cpp
*
* @brief Write spherical harmonics in ICGEM format.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2005-01-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write spherical harmonics in ICGEM format.
GROOPS uses this format as default but this program enables
the possibility to include comments and set the modelname.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileSphericalHarmonics.h"

/***** CLASS ***********************************/

/** @brief Write spherical harmonics in ICGEM format.
* @ingroup programsConversionGroup */
class PotentialCoefficients2Icgem
{
public:
  void run(Config &config);

  class Oscillation
  {
  public:
    FileName    fileNameCos, fileNameSin;
    std::string period;
    SphericalHarmonics harmCos, harmSin;
  };
};

GROOPS_REGISTER_PROGRAM(PotentialCoefficients2Icgem, SINGLEPROCESS, "write spherical harmonics in ICGEM format", Conversion, PotentialCoefficients)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, PotentialCoefficients2Icgem::Oscillation &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;
    readConfig(config, "inputfileCosPotentialCoefficients", var.fileNameCos, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSinPotentialCoefficients", var.fileNameSin, Config::MUSTSET,  "", "");
    readConfig(config, "period",                            var.period,      Config::MUSTSET, "1.0", "period of oscillation [year]");
    endSequence(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PotentialCoefficients2Icgem::run(Config &config)
{
  try
  {
    FileName                 fileNameOut;
    FileName                 fileNameStatic, fileNameTrend;
    FileName                 fileNameComments;
    std::string              modelname;
    std::vector<std::string> comment;
    std::vector<Oscillation> oscillation;
    std::string              choice;
    std::string              tides;
    UInt                     minDegree, maxDegree = INFINITYDEGREE;
    Double                   GM, R;
    Time                     time;

    readConfig(config, "outputfile",                     fileNameOut,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfilePotentialCoefficients", fileNameStatic,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTrend",                 fileNameTrend,    Config::OPTIONAL, "", "");
    readConfig(config, "oscillation",                    oscillation,      Config::OPTIONAL, "", "");
    readConfig(config, "inputfileComment",               fileNameComments, Config::OPTIONAL, "", "file containing comments for header");
    readConfig(config, "comment",                        comment,          Config::OPTIONAL, "", "comment in header");
    readConfig(config, "modelname",                      modelname,        Config::MUSTSET,  "", "name of the model");
    if(readConfigChoice(config, "tideSystem", choice,  Config::OPTIONAL, "", "tide system of model"))
    {
      if(readConfigChoiceElement(config, "zero_tide", choice)) tides = "zero_tide";
      if(readConfigChoiceElement(config, "tide_free", choice)) tides = "tide_free";
      endChoice(config);
    }
    readConfig(config, "minDegree",       minDegree,      Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",       maxDegree,      Config::OPTIONAL, "",  "");
    readConfig(config, "GM",              GM,             Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",               R,              Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "time",            time,           Config::OPTIONAL, "", "reference time");
    if(isCreateSchema(config)) return;

    // read potential coefficients
    // ---------------------------
    logStatus<<"read potential coefficients from file <"<<fileNameStatic<<">"<<Log::endl;
    SphericalHarmonics harm;
    readFileSphericalHarmonics(fileNameStatic, harm);
    harm = harm.get(maxDegree, minDegree, GM, R);

    SphericalHarmonics harmTrend;
    if(!fileNameTrend.empty())
    {
      logStatus<<"read potential coefficients from file <"<<fileNameTrend<<">"<<Log::endl;
      readFileSphericalHarmonics(fileNameTrend, harmTrend);
    }
    harmTrend = harmTrend.get(harm.maxDegree(), minDegree, GM, R);

    for(UInt i=0; i<oscillation.size(); i++)
    {
      logStatus<<"read potential coefficients from file <"<<oscillation.at(i).fileNameCos<<">"<<Log::endl;
      readFileSphericalHarmonics(oscillation.at(i).fileNameCos, oscillation.at(i).harmCos);
      oscillation.at(i).harmCos = oscillation.at(i).harmCos.get(harm.maxDegree(), minDegree, GM, R);

      logStatus<<"read potential coefficients from file <"<<oscillation.at(i).fileNameSin<<">"<<Log::endl;
      readFileSphericalHarmonics(oscillation.at(i).fileNameSin, oscillation.at(i).harmSin);
      oscillation.at(i).harmSin = oscillation.at(i).harmSin.get(harm.maxDegree(), minDegree, GM, R);
    }

    const Bool hasSigmas = (harm.sigma2cnm().size()) && ((quadsum(harm.sigma2cnm())+quadsum(harm.sigma2snm())) != 0);

    // description of data section
    // ---------------------------
    std::stringstream ssHeader;
    ssHeader<<"key     L    M         C                   S          ";
    if(hasSigmas)
      ssHeader<<"      sigma C             sigma S       ";
    if((!fileNameTrend.empty()) || oscillation.size())
      ssHeader<<" t0 [yyyymmdd]";
    if(oscillation.size())
      ssHeader<<"/period [y]";

    // Potentialkoeffizienten speichern
    // --------------------------------
    logStatus<<"writing potential coefficients to file <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);
    if(!fileNameComments.empty())
    {
      InFile commentFile(fileNameComments);
      std::string line;
      while(std::getline(commentFile, line))
        file<<line<<std::endl;
      file<<std::endl;
    }
    for(UInt i=0; i<comment.size(); i++)
      file<<comment.at(i)<<std::endl;
    if(comment.size())
      file<<std::endl;
    file<<"begin_of_head "<<std::string(ssHeader.str().size()-14, '=')<<std::endl<<std::endl;
    file<<"modelname              "<<modelname<<std::endl;
    file<<"product_type           gravity_field"<<std::endl;
    file<<"earth_gravity_constant "<<harm.GM()%"%16.10e"s<<std::endl;
    file<<"radius                 "<<harm.R() %"%16.10e"s<<std::endl;
    file<<"max_degree             "<<harm.maxDegree()<<std::endl;
    file<<"norm                   fully_normalized"<<std::endl;
    if(!tides.empty())
      file<<"tide_system            "<<tides<<std::endl;
    file<<"errors                 "<<((hasSigmas) ? "formal" : "no")<<std::endl;
    file<<std::endl;
    // end_of_head
    file<<ssHeader.str()<<std::endl;
    file<<"end_of_head "<<std::string(ssHeader.str().size()-12, '=')<<std::endl;

    std::string fdouble = "%20.12e";
    std::string fint    = "%5i";

    // data section
    // ------------
    for(UInt n=0; n<=harm.maxDegree(); n++)
      for(UInt m=0; m<=n; m++)
      {
        Bool hasTemporal = ((harmTrend.cnm()(n,m)!= 0.) || (harmTrend.snm()(n,m)!= 0.));
        for(UInt i=0; i<oscillation.size(); i++)
        {
          hasTemporal = hasTemporal || (oscillation.at(i).harmCos.cnm()(n,m)!= 0.) || (oscillation.at(i).harmCos.snm()(n,m)!= 0.);
          hasTemporal = hasTemporal || (oscillation.at(i).harmSin.cnm()(n,m)!= 0.) || (oscillation.at(i).harmSin.snm()(n,m)!= 0.);
        }

        // static coefficients
        file<<((hasTemporal) ? "gfct" : "gfc ");
        file<<n%fint<<m%fint;
        file<<harm.cnm()(n,m)%fdouble<<harm.snm()(n,m)%fdouble;
        if(hasSigmas)
          file<<sqrt(harm.sigma2cnm()(n,m))%fdouble<<sqrt(harm.sigma2snm()(n,m))%fdouble;
        if(hasTemporal)
          file<<" "<<time%"%y%m%d"s;
        file<<std::endl;

        // trend
        if((harmTrend.cnm()(n,m)!= 0.) || (harmTrend.snm()(n,m)!= 0.))
        {
          file<<"trnd"<<n%fint<<m%fint;
          file<<harmTrend.cnm()(n,m)%fdouble<<harmTrend.snm()(n,m)%fdouble;
          if(hasSigmas)
            file<<sqrt(harmTrend.sigma2cnm()(n,m))%fdouble<<sqrt(harmTrend.sigma2snm()(n,m))%fdouble;
          file<<std::endl;
        }

        // oscillations
        for(UInt i=0; i<oscillation.size(); i++)
          if((oscillation.at(i).harmCos.cnm()(n,m)!= 0.) ||
             (oscillation.at(i).harmCos.snm()(n,m)!= 0.) ||
             (oscillation.at(i).harmSin.cnm()(n,m)!= 0.) ||
             (oscillation.at(i).harmSin.snm()(n,m)!= 0.))
          {
            file<<"acos"<<n%fint<<m%fint;
            file<<oscillation.at(i).harmCos.cnm()(n,m)%fdouble<<oscillation.at(i).harmCos.snm()(n,m)%fdouble;
            if(hasSigmas)
              file<<sqrt(oscillation.at(i).harmCos.sigma2cnm()(n,m))%fdouble<<sqrt(oscillation.at(i).harmCos.sigma2snm()(n,m))%fdouble;
            file<<" "<<oscillation.at(i).period<<std::endl;

            file<<"asin"<<n%fint<<m%fint;
            file<<oscillation.at(i).harmSin.cnm()(n,m)%fdouble<<oscillation.at(i).harmSin.snm()(n,m)%fdouble;
            if(hasSigmas)
              file<<sqrt(oscillation.at(i).harmSin.sigma2cnm()(n,m))%fdouble<<sqrt(oscillation.at(i).harmSin.sigma2snm()(n,m))%fdouble;
            file<<" "<<oscillation.at(i).period<<std::endl;
          }
      } // for(n,m)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
