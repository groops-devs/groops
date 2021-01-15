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
    readConfig(config, "factor",                   factor,          Config::DEFAULT,  "1.0",   "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;

    cosPsiMin   /= DEFAULT_R*1e-3; if(cosPsiMin  >PI) cosPsiMin   = PI; cosPsiMin   = std::cos(cosPsiMin);
    cosPsiPrism /= DEFAULT_R*1e-3; if(cosPsiPrism>PI) cosPsiPrism = PI; cosPsiPrism = std::cos(cosPsiPrism);
    cosPsiLine  /= DEFAULT_R*1e-3; if(cosPsiLine >PI) cosPsiLine  = PI; cosPsiLine  = std::cos(cosPsiLine);
    cosPsiMax   /= DEFAULT_R*1e-3; if(cosPsiMax  >PI) cosPsiMax   = PI; cosPsiMax   = std::cos(cosPsiMax);

    // read rectangular grid
    // ---------------------
    GriddedDataRectangular grid;
    readFileGriddedData(gridName, grid);
    std::vector<Double> radius;
    grid.geocentric(lambda, phi, radius, dLambda, dPhi);
    rows = phi.size();
    cols = lambda.size();

    // evaluate upper and lower height
    // -------------------------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    expressionUpper->usedVariables(varList, usedVariables);
    expressionLower->usedVariables(varList, usedVariables);
    expressionRho  ->usedVariables(varList, usedVariables);
    addDataVariables(grid, varList, usedVariables);
    addVariable("area", varList);
    expressionUpper->simplify(varList);
    expressionLower->simplify(varList);
    expressionRho  ->simplify(varList);

    rLower = Matrix(rows,cols);
    rUpper = Matrix(rows,cols);
    rho    = Matrix(rows,cols);
    for(UInt z=0; z<rows; z++)
      for(UInt s=0; s<cols; s++)
      {
        evaluateDataVariables(grid, z, s, varList);
        varList["area"]->setValue( dLambda.at(s)*dPhi.at(z)*std::cos(phi.at(z)) ); // area
        rUpper(z,s) = radius.at(z) + expressionUpper->evaluate(varList);
        rLower(z,s) = radius.at(z) + expressionLower->evaluate(varList);
        rho(z,s)    = expressionRho->evaluate(varList);
      }

    // precompute sin & cos terms
    // --------------------------
    sinL = cosL = Vector(cols);
    for(UInt s=0; s<cols; s++)
    {
      sinL(s) = std::sin(lambda.at(s));
      cosL(s) = std::cos(lambda.at(s));
    }
    sinB = cosB = Vector(rows);
    for(UInt z=0; z<rows; z++)
    {
      sinB(z) = std::sin(phi.at(z));
      cosB(z) = std::cos(phi.at(z));
    }

    // Restrict computation to rectangle possible?
    // -------------------------------------------
    testRectangle = FALSE;
    if(cosPsiMax>-0.999)
    {
      testRectangle = TRUE;

      Bool ascending  = TRUE;
      Bool descending = TRUE;
      for(UInt i=1; i<rows; i++)
      {
        if(phi.at(i) > phi.at(i-1))
          descending = FALSE;
        if(phi.at(i) < phi.at(i-1))
          ascending = FALSE;
      }
      if(!ascending && !descending)
        testRectangle = FALSE;
      isPhiAscending = ascending;

      ascending  = descending = TRUE;
      for(UInt k=1; k<cols; k++)
      {
        if(lambda.at(k) > lambda.at(k-1))
          descending = FALSE;
        if(lambda.at(k) < lambda.at(k-1))
          ascending = FALSE;
      }
      if(!ascending && !descending)
        testRectangle = FALSE;
      isLambdaAscending = ascending;

      if(!testRectangle)
        logWarningOnce<<"Grid is not strictly ordered => acceleration of computation not possible"<<Log::endl;
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
    colsMax = cols;
    rowsMax = rows;
    if(!testRectangle)
      return;

    const Double delta = std::acos(cosPsiMax);

    const Double phiMin = point.phi() - delta;
    const Double phiMax = phiMin + 2*delta;
    if(isPhiAscending)
    {
      for(rowsMin=0;    (rowsMin<rows) && (phi.at(rowsMin) < phiMin); rowsMin++);
      for(rowsMax=rows; (rowsMax-->0)  && (phi.at(rowsMax) > phiMax););
    }
    else
    {
      for(rowsMin=0;    (rowsMin<rows) && (phi.at(rowsMin) > phiMax); rowsMin++);
      for(rowsMax=rows; (rowsMax-->0)  && (phi.at(rowsMax) < phiMin););
    }
    rowsMax++;

    const Double lambdaMin = point.lambda() - delta;
    const Double lambdaMax = lambdaMin + 2*delta;
    if(isLambdaAscending)
    {
      for(colsMin=0;    (colsMin<cols) && (lambda.at(colsMin) < lambdaMin); colsMin++);
      for(colsMax=cols; (colsMax-->0)  && (lambda.at(colsMax) > lambdaMax););
    }
    else
    {
      for(colsMin=0;    (colsMin<cols) && (lambda.at(colsMin) > lambdaMax); colsMin++);
      for(colsMax=cols; (colsMax-->0)  && (lambda.at(colsMax) < lambdaMin););
    }
    colsMax++;
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
      const Double   dL = this->dLambda.at(k);
      const Vector3d p  = rotaryZ(Angle(lambda.at(k))).rotate(point);

      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr)<0.001)
          continue;
        const Double dB    = this->dPhi.at(i);
        const Double cosB0 = this->cosB(i);
        const Double sinB0 = this->sinB(i);

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinB0*p.x() + cosB0*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosB0*p.x() + sinB0*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;


        // point mass
        // ----------
        if(cosPsi < cosPsiLine)
        {
          const Double r0 = (r2+r1)/2;
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          sum += rho(i,k)*r0*r0*cosB0*dL*dB*dr/l0;
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
          sum += (l2*r2 - l1*r1 + 3*rcosPsi * (l2-l1) + (3*rcosPsi*rcosPsi-r*r) * std::log((l2+r2-rcosPsi)/(l1+r1-rcosPsi)))*0.5*cosB0*dL*dB*rho(i,k);
          continue;
        }

        // prism - Franziska Wild-Pfeiffer Page 26
        // ---------------------------------------
        // coordinates relative to prism corners
        const Double r0 = (r2+r1)/2;
        x[0] +=  0.5*r0*dB;
        x[1]  = x[0]-r0*dB;
        y[0] +=  0.5*r0*dL*cosB0;
        y[1]  = y[0]-r0*dL*cosB0;
        z[0] -= r1;
        z[1] -= r2;

        for(UInt i=0; i<2; i++)
          for(UInt j=0; j<2; j++)
            for(UInt h=0; h<2; h++)
            {
              const Double l = std::sqrt(x[i]*x[i]+y[j]*y[j]+z[h]*z[h]);
              sum += std::pow(-1, i+j+h) * (x[i]*y[j]*std::log(z[h]+l) + y[j]*z[h]*std::log(x[i]+l) + z[h]*x[i]*std::log(y[j]+l)
                       -0.5*(x[i]*x[i]*std::atan(y[j]*z[h]/(x[i]*l)) + y[j]*y[j]*std::atan(z[h]*x[i]/(y[j]*l)) + z[h]*z[h]*std::atan(x[i]*y[j]/(z[h]*l)))) * rho(i,k);
            }
      } // for(i=0..rows)
    } //for(k=0..cols)

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
      const Double   dL = this->dLambda.at(k);
      const Vector3d p  = rotaryZ(Angle(lambda.at(k))).rotate(point);

      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr)<0.001)
          continue;
        const Double dB    = this->dPhi.at(i);
        const Double cosB0 = this->cosB(i);
        const Double sinB0 = this->sinB(i);

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinB0*p.x() + cosB0*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosB0*p.x() + sinB0*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;

        // point mass
        // ----------
        if(cosPsi < cosPsiLine)
        {
          const Double r0 = (r2+r1)/2;
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          sum -= (r*r-r0*rcosPsi)/(l0*l0*l0)*r0*r0*cosB0*dL*dB*dr*rho(i,k);
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
                 * rho(i,k) * cosB0*dL*dB;
          continue;
        }

        // computation point in local system
        const Vector3d p1(x[0], y[0], z[0]);

        // prism - Franziska Wild-Pfeiffer Page 26
        // ---------------------------------------
        // coordinates relative to prism corners
        const Double r0 = (r2+r1)/2;
        x[0] +=  0.5*r0*dB;
        x[1]  = x[0]-r0*dB;
        y[0] +=  0.5*r0*dL*cosB0;
        y[1]  = y[0]-r0*dL*cosB0;
        z[0] -= r1;
        z[1] -= r2;

        Vector3d dgxyz;
        for(UInt i=0; i<2; i++)
          for(UInt j=0; j<2; j++)
            for(UInt h=0; h<2; h++)
            {
              const Double sign = std::pow(-1, i+j+h);
              const Double l    = std::sqrt(x[i]*x[i]+y[j]*y[j]+z[h]*z[h]);
              const Double logx = std::log(x[i]+l);
              const Double logy = std::log(y[j]+l);
              const Double logz = std::log(z[h]+l);

              //gradients prism
              dgxyz.x() += sign * (y[j]*logz + z[h]*logy - x[i]*std::atan(y[j]*z[h]/(x[i]*l)));
              dgxyz.y() += sign * (z[h]*logx + x[i]*logz - y[j]*std::atan(z[h]*x[i]/(y[j]*l)));
              dgxyz.z() += sign * (x[i]*logy + y[j]*logx - z[h]*std::atan(x[i]*y[j]/(z[h]*l)));
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
      const Double   dL = this->dLambda.at(k);
      const Vector3d p  = rotaryZ(Angle(lambda.at(k))).rotate(point);

      Vector3d sum_local;
      for(UInt i=rowsMin; i<rowsMax; i++)
      {
        const Double r1 = rLower(i,k);
        const Double r2 = rUpper(i,k);
        const Double dr = r2-r1;
        if(std::fabs(dr)<0.001)
          continue;
        const Double dB    = this->dPhi.at(i);
        const Double cosB0 = this->cosB(i);
        const Double sinB0 = this->sinB(i);

        // transformation into local system
        Double x[2], y[2], z[2];
        x[0] = -sinB0*p.x() + cosB0*p.z();
        y[0] =  p.y();
        z[0] = z[1] = cosB0*p.x() + sinB0*p.z();
        const Double rcosPsi = z[0];
        const Double cosPsi  = rcosPsi/r;
        if((cosPsi < cosPsiMax) || (cosPsi > cosPsiMin))
          continue;

        Vector3d dg;
        if(cosPsi < cosPsiLine)
        {
          // point mass
          // ----------
          const Double r0 = (r2+r1)/2;
          z[0] -= r0;
          const Double l0 = std::sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
          const Double factor = -r0*r0*cosB0*dL*dB*dr/(l0*l0*l0);
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
          const Double dV_dr  = -r * logTerm2Term1*cosB0*dL*dB;
          const Double dV_dl1 = -(r1 + 3*rcosPsi + term3/term1)*0.5*cosB0*dL*dB;
          const Double dV_dl2 =  (r2 + 3*rcosPsi + term3/term2)*0.5*cosB0*dL*dB;
          dg      = dV_dr*dr + dV_dl2*dl2 + dV_dl1*dl1;
          dg.z() += (3*(l2-l1) + 6*rcosPsi * logTerm2Term1 - term3/term2 + term3/term1)*0.5*cosB0*dL*dB; // dV_drcosPsi = dV_dz
        }
        else
        {
          // prism - Franziska Wild-Pfeiffer Page 26
          // ---------------------------------------
          // coordinates relative to prism corners
          const Double r0 = (r2+r1)/2;
          x[0] +=  0.5*r0*dB;
          x[1]  = x[0]-r0*dB;
          y[0] +=  0.5*r0*dL*cosB0;
          y[1]  = y[0]-r0*dL*cosB0;
          z[0] -= r1;
          z[1] -= r2;

          for(UInt i=0; i<2; i++)
            for(UInt j=0; j<2; j++)
              for(UInt h=0; h<2; h++)
              {
                const Double sign = std::pow(-1, i+j+h);
                const Double l    = std::sqrt(x[i]*x[i]+y[j]*y[j]+z[h]*z[h]);
                const Double logx = std::log(x[i]+l);
                const Double logy = std::log(y[j]+l);
                const Double logz = std::log(z[h]+l);

                //gradients prism
                dg.x() += sign * (y[j]*logz + z[h]*logy - x[i]*std::atan(y[j]*z[h]/(x[i]*l)));
                dg.y() += sign * (z[h]*logx + x[i]*logz - y[j]*std::atan(z[h]*x[i]/(y[j]*l)));
                dg.z() += sign * (x[i]*logy + y[j]*logx - z[h]*std::atan(x[i]*y[j]/(z[h]*l)));
              }
        } // prim

        // rotate back (step 1)
        sum_local.x() += rho(i,k) * (-sinB0*dg.x() + cosB0*dg.z());
        sum_local.y() += rho(i,k) * dg.y();
        sum_local.z() += rho(i,k) * (cosB0*dg.x() + sinB0*dg.z());
      } // for(i=0..rows)

      // rotate back (step 2)
      sum += rotaryZ(Angle(-lambda.at(k))).rotate(sum_local);
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
                                         const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const
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

SphericalHarmonics GravityfieldTopography::sphericalHarmonics(const Time &/*time*/, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if(maxDegree==INFINITYDEGREE)
      throw(Exception("INFINITYDEGREE requested"));

    Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    Matrix cosm(cols, maxDegree+1);
    Matrix sinm(cols, maxDegree+1);
    for(UInt i=0; i<cols; i++)
      for(UInt m=0; m<=maxDegree; m++)
      {
        cosm(i,m) = std::cos(m*static_cast<Double>(lambda.at(i)));
        sinm(i,m) = std::sin(m*static_cast<Double>(lambda.at(i)));
      }

    for(UInt i=0; i<rows; i++)
    {
      Matrix f(cols, maxDegree+1);
      for(UInt k=0; k<cols; k++)
      {
        Double area = dLambda.at(k)*dPhi.at(i)*std::cos(phi.at(i));
        Double r1R  = rLower(i,k)/R;
        Double r2R  = rUpper(i,k)/R;
        Double r1Rn = factor * rho(i,k) * GRAVITATIONALCONSTANT/GM * area * R*R*R * std::pow(r1R, minDegree+3);
        Double r2Rn = factor * rho(i,k) * GRAVITATIONALCONSTANT/GM * area * R*R*R * std::pow(r2R, minDegree+3);

        for(UInt n=minDegree; n<=maxDegree; n++)
        {
          f(k,n) = (r2Rn-r1Rn)/((2*n+1)*(n+3));
          r1Rn *= r1R;
          r2Rn *= r2R;
        }
      }

      Matrix Pnm = SphericalHarmonics::Pnm(Angle(PI/2-phi.at(i)), 1.0, maxDegree);
      for(UInt n=minDegree; n<=maxDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          cnm(n,m) += Pnm(n,m) * inner(cosm.column(m), f.column(n));
          snm(n,m) += Pnm(n,m) * inner(sinm.column(m), f.column(n));
        }
    }

    return SphericalHarmonics(GM, R, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
