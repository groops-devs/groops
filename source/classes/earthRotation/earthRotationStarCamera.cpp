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

    quaternions = Matrix(arc.size(), 4);
    times.reserve(arc.size());
    for(UInt k = 0; k<arc.size(); k++)
    {
      times.push_back(arc.at(k).time);
      copy(arc.at(k).data().trans(), quaternions.row(k));
    }

    polynomial.init(polynomialDegree);
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
    Matrix interpolatedQuaternion = polynomial.interpolate({timeGPS}, times, quaternions);

    return Rotary3d(1./norm(interpolatedQuaternion)*interpolatedQuaternion.trans());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
