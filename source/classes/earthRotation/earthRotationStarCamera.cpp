/***********************************************/
/**
* @file earthRotationStarCamera.cpp
*
* @brief Earth rotation from quaternion (star camera) file.
* @see EarthRotation
*
* @author Andreas Kvas
* @date 2019-01-22
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/earthRotation/earthRotationStarCamera.h"

/***********************************************/

EarthRotationStarCamera::EarthRotationStarCamera(Config &config)
{
  try
  {
    FileName fileNameStarCamera;
    UInt polynomialDegree;
    readConfig(config, "inputfileStarCamera", fileNameStarCamera, Config::MUSTSET,  "",  "");
    readConfig(config, "interpolationDegree", polynomialDegree,   Config::DEFAULT,  "1", "degree of interpolation polynomial");
    if(isCreateSchema(config)) return;

    StarCameraArc arc = InstrumentFile::read(fileNameStarCamera);
    quaternions = arc.matrix().column(1, 4);
    polynomial.init(arc.times(), polynomialDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d EarthRotationStarCamera::rotaryMatrix(const Time &timeGPS) const
{
  try
  {
    Matrix q = polynomial.interpolate({timeGPS}, quaternions);
    return Rotary3d(1./norm(q)*q.trans());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
