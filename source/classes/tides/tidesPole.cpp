/***********************************************/
/**
* @file tidesPole.cpp
*
* @brief Centrifugal effect of polar motion.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2014-05-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileMeanPolarMotion.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/tides/tidesPole.h"

/***********************************************/

TidesPole::TidesPole(Config &config)
{
  try
  {
    FileName fileNameMeanPole;

    readConfig(config, "scale",                  scale,                  Config::DEFAULT, "1.333e-9", "");
    readConfig(config, "outPhase",               outPhase,               Config::DEFAULT, "0.0115",   "");
    readConfig(config, "inputfileMeanPole",      fileNameMeanPole,       Config::MUSTSET, "{groopsDataDir}/tides/secularPole2018.xml", "");
    readConfig(config, "horizontalDisplacement", horizontalDisplacement, Config::DEFAULT, "0.009",    "[m]");
    readConfig(config, "verticalDisplacement",   verticalDisplacement,   Config::DEFAULT, "0.033",    "[m]");
    readConfig(config, "factor",                 factor,                 Config::DEFAULT, "1.0",      "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;

    readFileMeanPolarMotion(fileNameMeanPole, meanPole);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TidesPole::pole(const Time &time, EarthRotationPtr earthRotation, Double &m1, Double &m2) const
{
  try
  {
    Double xBar, yBar;
    meanPole.compute(time, xBar, yBar);

    // pole (IERS2010, eq. (7.24))
    Double xp, yp, sp, deltaUT, LOD, X, Y, S;
    earthRotation->earthOrientationParameter(time, xp, yp, sp, deltaUT, LOD, X, Y, S);
    m1 =  (xp*RAD2DEG*3600 - xBar);
    m2 = -(yp*RAD2DEG*3600 - yBar);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

SphericalHarmonics TidesPole::sphericalHarmonics(const Time &time, const Rotary3d &/*rotEarth*/, EarthRotationPtr earthRotation, EphemeridesPtr /*ephemerides*/,
                                                 UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    Double m1, m2;
    pole(time, earthRotation, m1, m2);

    Matrix cnm(3, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(3, Matrix::TRIANGULAR, Matrix::LOWER);
    cnm(2,1) = -factor * scale * (m1 + outPhase*m2);
    snm(2,1) = -factor * scale * (m2 - outPhase*m1);

    return SphericalHarmonics(GM_Earth, R_Earth, cnm, snm).get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d TidesPole::deformation(const Time &time, const Vector3d &point, const Rotary3d &/*rotEarth*/, EarthRotationPtr earthRotation, EphemeridesPtr /*ephemerides*/,
                                Double /*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/) const
{
  try
  {
    Double m1, m2;
    pole(time, earthRotation, m1, m2);

    Angle lambda = point.lambda();
    Angle theta  = point.theta();

    Vector3d disp(horizontalDisplacement*cos(2.*theta)*(m1*cos(lambda)+m2*sin(lambda)),  // north
                  horizontalDisplacement*cos(theta)   *(m1*sin(lambda)-m2*cos(lambda)),  // east
                 -verticalDisplacement  *sin(2.*theta)*(m1*cos(lambda)+m2*sin(lambda))); // up

    return factor * localNorthEastUp(point).transform(disp);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TidesPole::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &/*rotEarth*/,
                            EarthRotationPtr rotation, EphemeridesPtr /*ephemerides*/, const std::vector<Double> &/*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/,
                            std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    for(UInt i=0; i<time.size(); i++)
    {
      Double m1, m2;
      pole(time.at(i), rotation, m1, m2);

      for(UInt k=0; k<point.size(); k++)
      {
        Angle lambda = point.at(k).lambda();
        Angle theta  = point.at(k).theta();

        Vector3d displ(horizontalDisplacement*cos(2.*theta)*(m1*cos(lambda)+m2*sin(lambda)),  // north
                       horizontalDisplacement*cos(theta)   *(m1*sin(lambda)-m2*cos(lambda)),  // east
                      -verticalDisplacement  *sin(2.*theta)*(m1*cos(lambda)+m2*sin(lambda))); // down

        disp.at(k).at(i) += factor * localNorthEastUp(point.at(k)).transform(displ);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
