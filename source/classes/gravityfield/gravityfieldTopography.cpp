/***********************************************/
/**
* @file gravityfieldTopography.cpp
*
* @brief Gravity field from topographic masses.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2013-02-10
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "parser/expressionParser.h"
#include "parser/dataVariables.h"
#include "parallel/parallel.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldTopography.h"

/***********************************************/

GravityfieldTopography::GravityfieldTopography(Config &config)
{
  try
  {
    FileName gridName;
    ExpressionVariablePtr expressionUpper, expressionLower, expressionRho;
    cosPsiMax = 1e99;

    readConfig(config, "inputfileGridRectangular", gridName,        Config::MUSTSET,  "",      "Digital Terrain Model");
    readConfig(config, "density",                  expressionRho,   Config::DEFAULT,  "2670",  "expression [kg/m**3]");
    readConfig(config, "radialUpperBound",         expressionUpper, Config::DEFAULT,  "data0", "expression (variables 'height', 'data', 'L', 'B' and, 'area' are taken from the gridded data");
    readConfig(config, "radialLowerBound",         expressionLower, Config::DEFAULT,  "0",     "expression (variables 'height', 'data', 'L', 'B' and, 'area' are taken from the gridded data");
    readConfig(config, "distanceMin",              cosPsiMin,       Config::DEFAULT,  "0",     "[km] min. influence distance (ignore near zone)");
    readConfig(config, "distancePrism",            cosPsiPrism,     Config::DEFAULT,  "15",    "[km] max. distance for prism formular");
    readConfig(config, "distanceLine",             cosPsiLine,      Config::DEFAULT,  "100",   "[km] max. distance for radial integration");
    readConfig(config, "distanceMax",              cosPsiMax,       Config::OPTIONAL, "",      "[km] max. influence distance (ignore far zone)");
    readConfig(config, "factor",                   factor,          Config::DEFAULT,  "1.0",   "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;

    cosPsiMin   = std::cos(std::min(PI, 1e3/DEFAULT_R * cosPsiMin));
    cosPsiPrism = std::cos(std::min(PI, 1e3/DEFAULT_R * cosPsiPrism));
    cosPsiLine  = std::cos(std::min(PI, 1e3/DEFAULT_R * cosPsiLine));
    cosPsiMax   = std::cos(std::min(PI, 1e3/DEFAULT_R * cosPsiMax));

    // read rectangular grid
    // ---------------------
    GriddedDataRectangular grid;
    readFileGriddedData(gridName, grid);
    std::vector<Double> radius;
    grid.geocentric(lambda0, phi0, radius);
    grid.areaElements(dLambda, dPhi);
    grid.cellBorders(lambda, phi);
    for(UInt i=0; i<phi.size(); i++)
      phi.at(i) = grid.ellipsoid(Angle(0.), Angle(phi.at(i)), 0.).phi();  // geocentric

    // evaluate upper and lower height
    // -------------------------------
    VariableList varList;
    addDataVariables(grid, varList);
    expressionUpper->simplify(varList);
    expressionLower->simplify(varList);
    expressionRho  ->simplify(varList);

    rLower = Matrix(phi0.size(), lambda0.size());
    rUpper = Matrix(phi0.size(), lambda0.size());
    rho    = Matrix(phi0.size(), lambda0.size());
    for(UInt i=0; i<phi0.size(); i++)
      for(UInt k=0; k<lambda0.size(); k++)
      {
        evaluateDataVariables(grid, i, k, varList);
        rUpper(i,k) = radius.at(i) + expressionUpper->evaluate(varList);
        rLower(i,k) = radius.at(i) + expressionLower->evaluate(varList);
        rho(i,k)    = expressionRho->evaluate(varList);
      }

    // precompute sin & cos terms
    // --------------------------
    sinL.resize(lambda0.size());
    cosL.resize(lambda0.size());
    for(UInt k=0; k<lambda0.size(); k++)
    {
      sinL.at(k) = std::sin(lambda0.at(k));
      cosL.at(k) = std::cos(lambda0.at(k));
    }
    sinPhi.resize(phi0.size());
    cosPhi.resize(phi0.size());
    for(UInt i=0; i<phi0.size(); i++)
    {
      sinPhi.at(i) = std::sin(phi0.at(i));
      cosPhi.at(i) = std::cos(phi0.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GravityfieldTopography::findRectangle(const Vector3d &point, UInt &colsMin, UInt &rowsMin, UInt &colsMax, UInt &rowsMax) const
{
  try
  {
    colsMin = 0;
    rowsMin = 0;
    colsMax = lambda0.size();
    rowsMax = phi0.size();
    if(cosPsiMax < -0.999)
      return;

    const Double phi    = point.phi();
    const Double lambda = point.lambda();
    const Double deltaP = std::acos(cosPsiMax);
    const Double deltaL = deltaP/std::cos(phi);

    rowsMin = std::distance(phi0.begin(),    std::find_if(phi0.begin(),    phi0.end(),    [&](auto x){return std::fabs(x-phi)    < deltaP;}));
    colsMin = std::distance(lambda0.begin(), std::find_if(lambda0.begin(), lambda0.end(), [&](auto x){return std::fabs(x-lambda) < deltaL;}));
    rowsMax = std::distance(std::find_if(phi0.rbegin(),    phi0.rend(),    [&](auto x){return std::fabs(x-phi)    < deltaP;}), phi0.rend());
    colsMax = std::distance(std::find_if(lambda0.rbegin(), lambda0.rend(), [&](auto x){return std::fabs(x-lambda) < deltaL;}), lambda0.rend());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

Double GravityfieldTopography::potential(const Time &/*time*/, const Vector3d &point) const
{
  try
  {
    const Double r = point.r();

    // find rectangle border
    UInt colsMin, rowsMin, colsMax, rowsMax;
    findRectangle(point, colsMin, rowsMin, colsMax, rowsMax);

    Double sum = 0.0;
    for(UInt k=colsMin; k<colsMax; k++)
    {
      const Vector3d p = rotaryZ(Angle(lambda0.at(k))).rotate(point);

      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr) < 0.001)
          continue;

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinPhi[i]*p.x() + cosPhi[i]*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosPhi[i]*p.x() + sinPhi[i]*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;


        // point mass
        // ----------
        if(cosPsi < cosPsiLine)
        {
          const Double r0 = 0.5*(r2+r1);
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          sum += rho(i,k)*r0*r0*dr*dLambda[k]*dPhi[i]/l0;
          continue;
        }

        // integration of a radial line
        // ----------------------------
        if(cosPsi < cosPsiPrism)
        {
          z[0] -= r1;
          z[1] -= r2;
          const Double hd = x[0]*x[0]+y[0]*y[0];
          const Double l1 = std::sqrt(hd + z[0]*z[0]);
          const Double l2 = std::sqrt(hd + z[1]*z[1]);
          sum += 0.5*(l2*r2 - l1*r1 + 3*rcosPsi * (l2-l1) + (3*rcosPsi*rcosPsi-r*r) * std::log((l2+r2-rcosPsi)/(l1+r1-rcosPsi)))*rho(i,k)*dLambda[k]*dPhi[i];
          continue;
        }

        // prism - Franziska Wild-Pfeiffer Page 26
        // ---------------------------------------
        // coordinates relative to prism corners
        const Double r0 = 0.5*(r2+r1);
        const Double dy = r0*dLambda[k]*dPhi[i] / std::fabs(phi[i+1]-phi[i]); // Volume of tesseroid: r0*r0*dr*dLambda[k]*dPhi[i] = dx*dy*dz
        x[0] += r0*std::fabs(phi[i]-phi0[i]);
        x[1]  = x[0]-r0*std::fabs(phi[i+1]-phi[i]);
        y[0] += 0.5*dy;
        y[1]  = y[0]-dy;
        z[0] -= r1;
        z[1] -= r2;

        Double sum2 = 0.;
        for(UInt ix=0; ix<2; ix++)
          for(UInt iy=0; iy<2; iy++)
            for(UInt iz=0; iz<2; iz++)
            {
              const Double l = std::sqrt(x[ix]*x[ix]+y[iy]*y[iy]+z[iz]*z[iz]);
              sum2 += std::pow(-1, ix+iy+iz) * (x[ix]*y[iy]*std::log(z[iz]+l) + y[iy]*z[iz]*std::log(x[ix]+l) + z[iz]*x[ix]*std::log(y[iy]+l)
                     -0.5*(x[ix]*x[ix]*std::atan(y[iy]*z[iz]/(x[ix]*l)) + y[iy]*y[iy]*std::atan(z[iz]*x[ix]/(y[iy]*l)) + z[iz]*z[iz]*std::atan(x[ix]*y[iy]/(z[iz]*l))));
            }
        sum += sum2 * rho(i,k);
      } // for(i=0..rows)
    } // for(k=0..cols)

    return factor * GRAVITATIONALCONSTANT * sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldTopography::radialGradient(const Time &/*time*/, const Vector3d &point) const
{
  try
  {
    const Double r = point.r();

    // find rectangle border
    UInt colsMin, rowsMin, colsMax, rowsMax;
    findRectangle(point, colsMin, rowsMin, colsMax, rowsMax);

    Double sum = 0.0;
    for(UInt k=colsMin; k<colsMax; k++)
    {
      const Vector3d p = rotaryZ(Angle(lambda0.at(k))).rotate(point);

      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr)<0.001)
          continue;

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinPhi[i]*p.x() + cosPhi[i]*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosPhi[i]*p.x() + sinPhi[i]*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;

        // point mass
        // ----------
        if(cosPsi < cosPsiLine)
        {
          const Double r0 = 0.5*(r2+r1);
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          sum -= (r*r-r0*rcosPsi)/(l0*l0*l0)*rho(i,k)*r0*r0*dr*dLambda[k]*dPhi[i];
          continue;
        }

        // integration of a radial line
        // ----------------------------
        if(cosPsi < cosPsiPrism)
        {
          z[0] -= r1;
          z[1] -= r2;
          const Double hd = x[0]*x[0]+y[0]*y[0];
          const Double l1 = std::sqrt(hd + z[0]*z[0]);
          const Double l2 = std::sqrt(hd + z[1]*z[1]);
          sum += (r1*r1*r1/l1 - r2*r2*r2/l2 + l2*r2 - l1*r1 + 3*rcosPsi * (l2-l1)
                  +(3*rcosPsi*rcosPsi-r*r) * std::log((l2+r2-rcosPsi)/(l1+r1-rcosPsi)))
                 * rho(i,k)*dLambda[k]*dPhi[i];
          continue;
        }

        // computation point in local system
        const Vector3d p1(x[0], y[0], z[0]);

        // prism - Franziska Wild-Pfeiffer Page 26
        // ---------------------------------------
        // coordinates relative to prism corners
        const Double r0 = 0.5*(r2+r1);
        const Double dy = r0*dLambda[k]*dPhi[i] / std::fabs(phi[i+1]-phi[i]); // Volume of tesseroid: r0^2*dr*dLambda[k]*dPhi[i] = dx*dy*dz
        x[0] += r0*std::fabs(phi[i]-phi0[i]);
        x[1]  = x[0]-r0*std::fabs(phi[i+1]-phi[i]);
        y[0] += 0.5*dy;
        y[1]  = y[0]-dy;
        z[0] -= r1;
        z[1] -= r2;

        Vector3d dgxyz;
        for(UInt ix=0; ix<2; ix++)
          for(UInt iy=0; iy<2; iy++)
            for(UInt iz=0; iz<2; iz++)
            {
              const Double sign = std::pow(-1, ix+iy+iz);
              const Double l    = std::sqrt(x[ix]*x[ix]+y[iy]*y[iy]+z[iz]*z[iz]);
              const Double logx = std::log(x[ix]+l);
              const Double logy = std::log(y[iy]+l);
              const Double logz = std::log(z[iz]+l);

              //gradients prism
              dgxyz.x() += sign * (y[iy]*logz + z[iz]*logy - x[ix]*std::atan(y[iy]*z[iz]/(x[ix]*l)));
              dgxyz.y() += sign * (z[iz]*logx + x[ix]*logz - y[iy]*std::atan(z[iz]*x[ix]/(y[iy]*l)));
              dgxyz.z() += sign * (x[ix]*logy + y[iy]*logx - z[iz]*std::atan(x[ix]*y[iy]/(z[iz]*l)));
            }
        sum += rho(i,k) * inner(dgxyz, p1); // projection into radial direction

      } // for(i=0..rows)
    } //for(k=0..cols)

    return factor * GRAVITATIONALCONSTANT * sum/r;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GravityfieldTopography::gravity(const Time &/*time*/, const Vector3d &point) const
{
  try
  {
    const Double r = point.r();

    // find rectangle border
    UInt colsMin, rowsMin, colsMax, rowsMax;
    findRectangle(point, colsMin, rowsMin, colsMax, rowsMax);

    Vector3d sum;
    for(UInt k=colsMin; k<colsMax; k++)
    {
      const Vector3d p = rotaryZ(Angle(lambda0.at(k))).rotate(point);

      Vector3d sum_local;
      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr) < 0.001)
          continue;

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinPhi[i]*p.x() + cosPhi[i]*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosPhi[i]*p.x() + sinPhi[i]*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;

        Vector3d dg;
        if(cosPsi < cosPsiLine)
        {
          // point mass
          // ----------
          const Double r0 = 0.5*(r2+r1);
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          const Double factor = -r0*r0*dr*dLambda[k]*dPhi[i]/(l0*l0*l0);
          dg = Vector3d(factor*x[0], factor*y[0], factor*z[0]);
        }
        else if(cosPsi < cosPsiPrism)
        {
          // integration of a radial line
          // ----------------------------
          const Vector3d dr(x[0]/r, y[0]/r, z[0]/r);
          z[0] -= r1;
          z[1] -= r2;
          const Double hd = x[0]*x[0]+y[0]*y[0];
          const Double l1 = std::sqrt(hd + z[0]*z[0]);
          const Double l2 = std::sqrt(hd + z[1]*z[1]);
          const Vector3d dl1(x[0]/l1, y[0]/l1, z[0]/l1);
          const Vector3d dl2(x[0]/l2, y[0]/l2, z[1]/l2);
          const Double term1 = l1+r1-rcosPsi;
          const Double term2 = l2+r2-rcosPsi;
          const Double term3 = 3*rcosPsi*rcosPsi-r*r;
          const Double logTerm2Term1 = std::log(term2/term1);
          const Double dV_dr  = -r * logTerm2Term1*dLambda[k]*dPhi[i];
          const Double dV_dl1 = -(r1 + 3*rcosPsi + term3/term1)*0.5*dLambda[k]*dPhi[i];
          const Double dV_dl2 =  (r2 + 3*rcosPsi + term3/term2)*0.5*dLambda[k]*dPhi[i];
          dg      = dV_dr*dr + dV_dl2*dl2 + dV_dl1*dl1;
          dg.z() += (3*(l2-l1) + 6*rcosPsi * logTerm2Term1 - term3/term2 + term3/term1)*0.5*dLambda[k]*dPhi[i]; // dV_drcosPsi = dV_dz
        }
        else
        {
          // prism - Franziska Wild-Pfeiffer Page 26
          // ---------------------------------------
          // coordinates relative to prism corners
          const Double r0 = 0.5*(r2+r1);
          const Double dy = r0*dLambda[k]*dPhi[i] / std::fabs(phi[i+1]-phi[i]); // Volume of tesseroid: r0*r0*dr*dLambda[k]*dPhi[i] = dx*dy*dz
          x[0] += r0*std::fabs(phi[i]-phi0[i]);
          x[1]  = x[0]-r0*std::fabs(phi[i+1]-phi[i]);
          y[0] += 0.5*dy;
          y[1]  = y[0]-dy;
          z[0] -= r1;
          z[1] -= r2;

          for(UInt ix=0; ix<2; ix++)
            for(UInt iy=0; iy<2; iy++)
              for(UInt iz=0; iz<2; iz++)
              {
                const Double sign = std::pow(-1, ix+iy+iz);
                const Double l    = std::sqrt(x[ix]*x[ix]+y[iy]*y[iy]+z[iz]*z[iz]);
                const Double logx = std::log(x[ix]+l);
                const Double logy = std::log(y[iy]+l);
                const Double logz = std::log(z[iz]+l);

                //gradients prism
                dg.x() += sign * (y[iy]*logz + z[iz]*logy - x[ix]*std::atan(y[iy]*z[iz]/(x[ix]*l)));
                dg.y() += sign * (z[iz]*logx + x[ix]*logz - y[iy]*std::atan(z[iz]*x[ix]/(y[iy]*l)));
                dg.z() += sign * (x[ix]*logy + y[iy]*logx - z[iz]*std::atan(x[ix]*y[iy]/(z[iz]*l)));
              }
        } // prim

        // rotate back (step 1)
        sum_local.x() += rho(i,k) * (-sinPhi[i]*dg.x() + cosPhi[i]*dg.z());
        sum_local.y() += rho(i,k) * dg.y();
        sum_local.z() += rho(i,k) * (cosPhi[i]*dg.x() + sinPhi[i]*dg.z());
      } // for(i=0..rows)

      // rotate back (step 2)
      sum += rotaryZ(Angle(-lambda0.at(k))).rotate(sum_local);
    } //for(k=0..cols)

    return factor * GRAVITATIONALCONSTANT * sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Tensor3d GravityfieldTopography::gravityGradient(const Time &/*time*/, const Vector3d &/*point*/) const
{
  try
  {
    throw(Exception("not implemented yet"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GravityfieldTopography::deformation(const Time &/*time*/, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  SphericalHarmonics harm = sphericalHarmonics(Time(), hn.rows()-1, 0, DEFAULT_GM, DEFAULT_R);
  return harm.deformation(point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldTopography::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                         const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  SphericalHarmonics harm = sphericalHarmonics(Time(), hn.rows()-1, 0, DEFAULT_GM, DEFAULT_R);
  for(UInt k=0; k<point.size(); k++)
  {
    Vector3d d = harm.deformation(point.at(k), gravity.at(k), hn, ln);
    for(UInt i=0; i<time.size(); i++)
      disp.at(k).at(i) = d;
  }
}

/***********************************************/

void GravityfieldTopography::variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const
{
}

/***********************************************/

SphericalHarmonics GravityfieldTopography::sphericalHarmonics(const Time &/*time*/, UInt /*maxDegree*/, UInt /*minDegree*/, Double /*GM*/, Double /*R*/) const
{
  try
  {
    throw(Exception("Conversion to spherical harmonics not implemented\nPlease use program GriddedTopography2PotentialCoefficients before."));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
