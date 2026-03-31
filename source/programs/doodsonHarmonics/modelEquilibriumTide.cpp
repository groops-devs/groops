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
The spherical harmonics expansion up to \config{maxDegree} with \config{GM} and \config{R}
is estimated using a least squares adjustment.

The \configFile{inputfileDensityGrid}{griddedData} must be a global regular grid with the
vertically averaged seawater density over the ocean and zero over land.

It takes iteratively self attraction and loading into account using the Love numbers
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
#include "files/fileGriddedData.h"
#include "files/fileSphericalHarmonics.h"
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
    FileName  fileNameCoeff;
    FileName  fileNameDensityGrid;
    FileName  loadName, deformationName;
    UInt      maxDegree;
    Double    GM, R;
    Double    TGP, k20, k20plus, h20_0, h20_2;
    UInt      iterationCount;

    readConfig(config, "outputfilePotentialCoefficients",    fileNameCoeff,       Config::MUSTSET, "",                "includes the loading");
    readConfig(config, "maxDegree",                          maxDegree,           Config::DEFAULT, "120",             "");
    readConfig(config, "GM",                                 GM,                  Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                  R,                   Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "inputfileDensityGrid",               fileNameDensityGrid, Config::MUSTSET, "",                "[kg/m^3] density of sea water, zero over land");
    readConfig(config, "tideGeneratingPotential",            TGP,                 Config::MUSTSET, "-0.0856534056",   "[m^2/s^2]");
    readConfig(config, "k20",                                k20,                 Config::DEFAULT, "0.30190",         "earth tide love number");
    readConfig(config, "k20plus",                            k20plus,             Config::DEFAULT, "-0.00089",        "earth tide love number");
    readConfig(config, "h20_0",                              h20_0,               Config::DEFAULT, "0.6078",          "earth tide love number");
    readConfig(config, "h20_2",                              h20_2,               Config::DEFAULT, "-0.0006",         "earth tide love number");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   loadName,            Config::MUSTSET, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName,     Config::MUSTSET, "{groopsDataDir}/loading/deformationLoveNumbers_CE_Gegout97.txt", "");
    readConfig(config, "iterationCount",                     iterationCount,      Config::DEFAULT, "5",               "");
    if(isCreateSchema(config)) return;

    // =======================================

    // init grid
    // ---------
    logStatus<<"read density grid from file <"<<fileNameDensityGrid<<">"<<Log::endl;
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameDensityGrid, grid);
    MiscGriddedData::printStatistics(grid);
    Matrix density;
    std::swap(density, grid.values.at(0));

    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius, dLambda, dPhi;
    grid.geocentric(lambda, phi, radius);
    grid.areaElements(dLambda, dPhi);

    // check assumption: same area for all longitudes
    Bool failed = FALSE;
    for(UInt k=0; k<lambda.size()-1; k++)
      if(!failed && (std::fabs(dLambda.at(k+1)-dLambda.at(k)) > 1e-3*std::fabs(dLambda.at(k))))
      {
        logWarningOnce<<"assumption of the same area/weight for all longitudes not fulfilled"<<Log::endl;
        failed = TRUE;
      }

    Double totalArea = 0;
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        totalArea += dPhi.at(i) * dLambda.at(k) * std::pow(radius.at(i), 2);

    Double oceanArea = 0;
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        if(density(i, k) > 0)
          oceanArea += dPhi.at(i) * dLambda.at(k) * std::pow(radius.at(i), 2);
    logInfo<<"  ocean area:   "<<oceanArea/totalArea*100<<"% of earth surface"<<Log::endl;

    Double meanDensity = 0;
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        if(density(i, k) > 0)
          meanDensity += density(i, k) * dPhi.at(i) * dLambda.at(k) * std::pow(radius.at(i), 2);
    meanDensity /= oceanArea;
    logInfo<<"  mean density: "<<meanDensity<<" kg/m^2"<<Log::endl;

    // =======================================

    // deformation load love numbers
    // -----------------------------
    Matrix load_kn, load_hnln;
    readFileMatrix(deformationName, load_hnln);
    readFileMatrix(loadName,        load_kn);

    // =======================================

    // init tgp
    // --------
    logStatus<<"Tide raisung potential: TGP + Solid Earth Tide"<<Log::endl;
    Matrix height0(phi.size(), lambda.size());
    Parallel::forEach(phi.size(), [&](UInt i)
    {
      const Vector3d p       = polar(Angle(0), phi.at(i), radius.at(i));
      const Double   gravity = Planets::normalGravity(p);
      Matrix         Cnm, Snm;
      SphericalHarmonics::CnmSnm((1./radius.at(i))*p, 4, Cnm, Snm);

      Double h0 =     std::pow(radius.at(i)/DEFAULT_R, 2) * TGP * Cnm(2,0) / gravity;           // direct tide
      h0 += k20     * std::pow(DEFAULT_R/radius.at(i), 3) * TGP * Cnm(2,0) / gravity;           // earth tide potential
      h0 += k20plus * std::pow(DEFAULT_R/radius.at(i), 5) * TGP * Cnm(4,0) / gravity;           // earth tide potential
      h0 -= (h20_0 + h20_2*Cnm(2,0)/std::sqrt(5)) * std::pow(DEFAULT_R, 2)/GM * TGP * Cnm(2,0); // displacement
      for(UInt k=0; k<lambda.size(); k++)
        if(density(i, k) > 0)
          height0(i, k) = h0;
    }, comm);
    Parallel::reduceSum(height0, 0, comm);
    Parallel::broadCast(height0, 0, comm);

    // =======================================

    Matrix tide = height0;
    Double massDiff = 0;
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        if(density(i, k) > 0)
          massDiff += tide(i, k) * density(i, k) * dPhi.at(i) * dLambda.at(k) * std::pow(radius.at(i), 2);
    logInfo<<"  total Mass change =  "<<1000*massDiff/meanDensity/oceanArea<<" mm"<<Log::endl;
    for(UInt i=0; i<phi.size(); i++)
      for(UInt k=0; k<lambda.size(); k++)
        if(density(i, k) > 0)
          tide(i,k) -= massDiff/meanDensity/oceanArea;

    // ========================================================

    // precompute cos(m*lambda), sin(m*lambda)
    Matrix cosml(lambda.size(), maxDegree+1);
    Matrix sinml(lambda.size(), maxDegree+1);
    for(UInt m=0; m<=maxDegree; m++)
      for(UInt k=0; k<lambda.size(); k++)
      {
        cosml(k, m) = std::cos(m*static_cast<Double>(lambda.at(k)));
        sinml(k, m) = std::sin(m*static_cast<Double>(lambda.at(k)));
      }

    // system of normal equations (order by order)
    std::vector<Matrix> N, n;
    N.push_back(Matrix(maxDegree+1, Matrix::SYMMETRIC));
    n.push_back(Matrix(maxDegree+1, 1));
    for(UInt m=1; m<=maxDegree; m++)
    {
      N.push_back(Matrix(2*(maxDegree+1-m), Matrix::SYMMETRIC));
      n.push_back(Matrix(2*(maxDegree+1-m), 1));
    }

    SphericalHarmonics harm;
    for(UInt iter=0; iter<iterationCount; iter++)
    {
      logStatus<<iter+1<<". iteration"<<Log::endl;
      logStatus<<"  compute spherical harmonics"<<Log::endl;
      // least squares adjustment order by order
      for(UInt m=0; m<n.size(); m++)
        n.at(m).setNull();

      Parallel::forEach(phi.size(), [&](UInt i)
      {
        Vector l = tide.row(i).trans();
        for(UInt k=0; k<lambda.size(); k++)
          l(k) *= density(i, k);

        // legendre functions with kernel coefficients
        Matrix Pnm = SphericalHarmonics::Pnm(Angle(0.5*PI-phi.at(i)), radius.at(i)/R, maxDegree);
        for(UInt n=0; n<=maxDegree; n++)
          Pnm.slice(n, 0, 1, n+1) *= GM/R/(4.*PI*GRAVITATIONALCONSTANT*radius.at(i)) * (2*n+1);

        const Double weight = dPhi.at(i) * dLambda.at(0)/(4*PI) * std::pow(radius.at(i)/R, 2); // assume same area for all longitudes
        Matrix A = cosml.column(0) * Pnm.column(0).trans();
        if(iter == 0)
          rankKUpdate(weight, A, N.at(0));
        matMult(weight, A.trans(), l, n.at(0));
        for(UInt m=1; m<=maxDegree; m++)
        {
          Matrix A(lambda.size(), 2*(maxDegree+1-m));
          matMult(1.0, cosml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(0, maxDegree+1-m));
          matMult(1.0, sinml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(maxDegree+1-m, maxDegree+1-m));
          if(iter == 0)
            rankKUpdate(weight, A, N.at(m));
          matMult(weight, A.trans(), l, n.at(m));
        }
      }, comm);

      if(iter == 0)
        for(UInt m=0; m<=maxDegree; m++)
          Parallel::reduceSum(N.at(m), 0, comm);
      for(UInt m=0; m<=maxDegree; m++)
        Parallel::reduceSum(n.at(m), 0, comm);

      // solve normals
      // -------------
      Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      if(Parallel::isMaster(comm))
      {
        copy(solve(Matrix(N.at(0)), n.at(0)), cnm.column(0));
        for(UInt m=1; m<=maxDegree; m++)
        {
          const Vector x = solve(Matrix(N.at(m)), n.at(m));
          copy(x.row(0,             maxDegree+1-m), cnm.slice(m, m, maxDegree+1-m, 1));
          copy(x.row(maxDegree+1-m, maxDegree+1-m), snm.slice(m, m, maxDegree+1-m, 1));
        }
      } // if(Parallel::isMaster(comm))
      Parallel::broadCast(cnm, 0, comm);
      Parallel::broadCast(snm, 0, comm);
      logInfo<<"  c00 change =  "<<1000*R*cnm(0,0)<<" mm"<<Log::endl;
      cnm(0,0) = 0;
      harm = SphericalHarmonics(GM, R, cnm, snm);

      // synthesis spherical harmonics
      // -----------------------------
      logStatus<<"  compute loading and self attraction"<<Log::endl;
      Matrix lsa(phi.size(), lambda.size());
      Parallel::forEach(phi.size(), [&](UInt i) // Compute Legendre functions for each phi (row)
      {
        const Double gravity = Planets::normalGravity(polar(Angle(0), phi.at(i), radius.at(i)));
        Matrix Pnm = SphericalHarmonics::Pnm(Angle(PI/2-phi.at(i)), radius.at(i)/R, maxDegree);
        for(UInt n=0; n<=maxDegree; n++)
          Pnm.slice(n,0,1,n+1) *= GM/R*(1.+load_kn(n,0)-load_hnln(n,0))/gravity;

        Vector cm(maxDegree+1), sm(maxDegree+1);
        cm(0) = inner(cnm.column(0), Pnm.column(0));
        for(UInt m=1; m<=maxDegree; m++)
        {
          cm(m) = inner(cnm.slice(m, m, maxDegree-m+1, 1), Pnm.slice(m, m, maxDegree-m+1, 1));
          sm(m) = inner(snm.slice(m, m, maxDegree-m+1, 1), Pnm.slice(m, m, maxDegree-m+1, 1));
        }
        Vector row = cosml * cm + sinml * sm;

        for(UInt k=0; k<lambda.size(); k++)
          if(density(i, k) > 0)
            lsa(i, k) += row(k);
      }, comm);
      Parallel::reduceSum(lsa, 0, comm);

      if(Parallel::isMaster(comm))
      {
        // compute new tide
        Matrix tideOld = tide;
        tide = height0 + lsa;

        // mass conservation
        Double massDiff = 0;
        for(UInt i=0; i<phi.size(); i++)
          for(UInt k=0; k<lambda.size(); k++)
            if(density(i, k) > 0)
              massDiff += tide(i, k) * density(i, k) * dPhi.at(i) * dLambda.at(k) * std::pow(radius.at(i), 2);
        for(UInt i=0; i<phi.size(); i++)
          for(UInt k=0; k<lambda.size(); k++)
            if(density(i, k) > 0)
              tide(i,k) -= massDiff/meanDensity/oceanArea;

        logInfo<<"  max. change in iteration =  "<<1000*maxabs(tide-tideOld)<<" mm"<<Log::endl;
      }

      Parallel::broadCast(tide, 0, comm);
    } // for(iter)

    // =======================================

    // add loading potential
    // ---------------------
    if(Parallel::isMaster(comm))
      for(UInt n=0; n<=maxDegree; n++)
      {
        harm.cnm().row(n) *= 1.+load_kn(n,0);
        harm.snm().row(n) *= 1.+load_kn(n,0);
      }

    // write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"writing potential coefficients to file <"<<fileNameCoeff<<">"<<Log::endl;
      writeFileSphericalHarmonics(fileNameCoeff, harm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
