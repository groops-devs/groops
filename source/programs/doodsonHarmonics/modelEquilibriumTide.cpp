/***********************************************/
/**
* @file modelEquilibriumTide.cpp
*
* @brief Equilibrium tide, taking the effect of loading and self attraction into account.
*
* @author Torsten Mayer-Guerr
* @date 2010-05-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes the equilibrium ocean tide of the long periodic \config{tideGeneratingPotential}.

The ocean surface is represented by \configClass{gridOcean}{gridType} and the gravitational
effect is numerical integrated to spherical harmonics using \config{maxDegree}, \config{GM},
and \config{R}.

It takes self attraction and loading into account using the Love numbers
\configFile{inputfilePotentialLoadLoveNumber}{matrix} and
\configFile{inputfileDeformationLoadLoveNumber}{matrix}.

Additionally the effects of the solid Earth tide are considered,
both the gravitational (Love numbers \config{k20}, \config{k20plus})
and the geometrical (Love numbers \config{h20,0}, \config{h20,2}) effect.

See also \program{PotentialCoefficients2DoodsonHarmonics}.

\fig{!hb}{0.8}{modelEquilibriumTide}{fig:modelEquilibriumTide}{Equilibrium tide of SA constituent}
)";

/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "base/planets.h"
#include "files/fileMatrix.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Equilibrium tide, taking the effect of loading and self attraction into account.
* @ingroup programsGroup */
class ModelEquilibriumTide
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(ModelEquilibriumTide, PARALLEL, "equilibrium tide, taking the effect of loading and self attraction into account.", DoodsonHarmonics, Simulation)

/***********************************************/

void ModelEquilibriumTide::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName  potentialCoefficientsName;
    FileName  loadName, deformationName;
    GridPtr   gridOcean;
    UInt      maxDegree;
    Double    TGP, k20, k20plus, h20_0, h20_2;
    Double    density;
    Double    GM, R;

    readConfig(config, "outputfilePotentialCoefficients",    potentialCoefficientsName, Config::MUSTSET, "",       "includes the loading");
    readConfig(config, "gridOcean",                          gridOcean,       Config::MUSTSET,  "",                "");
    readConfig(config, "maxDegree",                          maxDegree,       Config::DEFAULT,  "120",             "");
    readConfig(config, "GM",                                 GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                  R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "density",                            density,         Config::DEFAULT,  "1025",            "[kg/m**3] density of sea water");
    readConfig(config, "tideGeneratingPotential",            TGP,             Config::MUSTSET,  "-0.0856534056",   "[m**2/s**2]");
    readConfig(config, "k20",                                k20,             Config::DEFAULT,  "0.30190",         "earth tide love number");
    readConfig(config, "k20plus",                            k20plus,         Config::DEFAULT,  "-0.00089",        "earth tide love number");
    readConfig(config, "h20_0",                              h20_0,           Config::DEFAULT,  "0.6078",          "earth tide love number");
    readConfig(config, "h20_2",                              h20_2,           Config::DEFAULT,  "-0.0006",         "earth tide love number");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   loadName,        Config::MUSTSET,  "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName, Config::MUSTSET,  "{groopsDataDir}/loading/deformationLoveNumbers_CE_Gegout97.txt", "");
    if(isCreateSchema(config)) return;

    // =======================================

    // init Ocean
    // ----------
    logStatus<<"init grid"<<Log::endl;
    std::vector<Vector3d> points = gridOcean->points();
    std::vector<Double>   areas  = gridOcean->areas();
    Double totalArea = 0;
    for(UInt i=0; i<points.size(); i++)
        totalArea += areas.at(i);
    logInfo<<"  number of points: "<<points.size()<<Log::endl;
    logInfo<<"  total area:       "<<totalArea/(4*PI)*100<<"% of earth surface"<<Log::endl;

    // deformation load love numbers
    // -----------------------------
    Matrix load_kn, load_hnln;
    readFileMatrix(deformationName, load_hnln);
    readFileMatrix(loadName,        load_kn);

    // =======================================

    // init tgp
    // --------
    logStatus<<"direct tide and earth tide"<<Log::endl;
    std::vector<Double> height0(points.size(), 0.0);
    Parallel::forEach(height0, [&](UInt i)
    {
      const Double r       = points.at(i).r();
      const Double gravity = Planets::normalGravity(points.at(i));
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm((1./r)*points.at(i), 4, Cnm, Snm);

      // direct tide
      Double height0 = std::pow(r/DEFAULT_R, 2) * TGP * Cnm(2,0) / gravity;
      // earth tide potential
      height0 += k20     * std::pow(DEFAULT_R/r, 3) * TGP * Cnm(2,0) / gravity;
      height0 += k20plus * std::pow(DEFAULT_R/r, 5) * TGP * Cnm(4,0) / gravity;
      // displacement
      const Double h20 = h20_0 + h20_2*Cnm(2,0)/std::sqrt(5);
      height0 -= h20 * std::pow(DEFAULT_R, 2)/GM * TGP * Cnm(2,0);
      return height0;
    }, comm);
    Parallel::broadCast(height0, 0, comm);

    // =======================================

    logStatus<<"loading and self attraction"<<Log::endl;
    std::vector<Double> tide = height0;
    Double massDiff = 0;
    for(UInt i=0; i<points.size(); i++)
      massDiff += tide.at(i) * areas.at(i);
    for(UInt i=0; i<points.size(); i++)
      tide.at(i) -= massDiff/totalArea;

    SphericalHarmonics harm;
    for(UInt iter=0; iter<3; iter++)
    {
      logStatus<<"compute spherical harmonics"<<Log::endl;
      Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Parallel::forEach(points.size(), [&](UInt i)
      {
        Matrix Cnm, Snm;
        SphericalHarmonics::CnmSnm((1/R)*points.at(i), maxDegree, Cnm, Snm);
        const Double factor = GRAVITATIONALCONSTANT*density*points.at(i).r()*R/GM;
        for(UInt n=0; n<=maxDegree; n++)
        {
          axpy(factor/(2*n+1) * tide.at(i) * areas.at(i), Cnm.row(n), cnm.row(n));
          axpy(factor/(2*n+1) * tide.at(i) * areas.at(i), Snm.row(n), snm.row(n));
        }
      }, comm);
      Parallel::reduceSum(cnm, 0, comm);
      Parallel::reduceSum(snm, 0, comm);
      Parallel::broadCast(cnm, 0, comm);
      Parallel::broadCast(snm, 0, comm);
      logInfo<<"  total Mass change =  "<<R*cnm(0,0)<<" m"<<Log::endl;
      cnm(0,0) = 0;
      harm = SphericalHarmonics(GM, R, cnm, snm);

      logStatus<<"compute loading and self attraction"<<Log::endl;
      std::vector<Double> lsa(points.size());
      Parallel::forEach(lsa, [&](UInt i)
      {
        const Vector Yn = harm.Yn(points.at(i));
        const Double gravity = Planets::normalGravity(points.at(i));
        Double sum = 0;
        for(UInt n=0; n<Yn.rows(); n++)
          sum += (1.+load_kn(n,0)-load_hnln(n,0))/gravity * Yn(n);
        return sum;
      }, comm);
      Parallel::broadCast(lsa, 0, comm);

      // compute new tide
      std::vector<Double> tideOld = tide;
      for(UInt i=0; i<points.size(); i++)
        tide.at(i) = height0.at(i) + lsa.at(i);

      // mass conservation
      Double massDiff = 0;
      for(UInt i=0; i<points.size(); i++)
        massDiff += tide.at(i) * areas.at(i);
      for(UInt i=0; i<points.size(); i++)
        tide.at(i) -= massDiff/totalArea;

      Double maxChange = 0;
      for(UInt i=0; i<points.size(); i++)
        maxChange = std::max(maxChange, std::fabs(tide.at(i)-tideOld.at(i)));
      logInfo<<"  max. change in iteration =  "<<maxChange<<" m"<<Log::endl;
    }

    // add loading potential
    // ---------------------
    for(UInt n=0; n<=maxDegree; n++)
    {
      harm.cnm().row(n) *= 1.+load_kn(n,0);
      harm.snm().row(n) *= 1.+load_kn(n,0);
    }

    // write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"writing potential coefficients to file <"<<potentialCoefficientsName<<">"<<Log::endl;
      writeFileSphericalHarmonics(potentialCoefficientsName, harm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
