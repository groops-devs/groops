/***********************************************/
/**
* @file fileEphemerides.h
*
* @brief Ephemerides of sun, moon, and planets.
*
* @author Torsten Mayer-Guerr
* @date 2020-08-11
*
*/
/***********************************************/

#ifndef __GROOPS_FILEEPHEMERIDES__
#define __GROOPS_FILEEPHEMERIDES__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Ephemerides
static const char *docstringEphemerides = R"(
Ephemerides of sun, moon, and planets stored as coefficients of Chebyshev polynomials.
Used in \configClass{Ephemerides:jpl}{ephemeridesType:jpl}.

See also \program{JplAscii2Ephemerides}.
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_EPHEMERIDES_TYPE    = "ephemerides";
constexpr UInt    FILE_EPHEMERIDES_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** CLASS ***********************************/

/** @brief Ephemerides of sun, moon, and planets. */
class InFileEphemerides
{
  InFileArchive     file;
  Double            earthMoonRatio;
  std::vector<Time> times;                // interval boundaries
  std::vector<UInt> subintervals;         // for each body: number of subintervals
  std::vector<UInt> components;           // for each body: number of components: 3 for x,y,z
  std::vector<UInt> degree;               // for each body: polynomial degree
  std::streampos    seekStart;
  std::streamoff    seekSize;
  UInt              idInterval;           // next interval to read
  std::vector<std::vector<Matrix>> coeff; // for each body, subinterval: Matrix(components, degree+1)

  void readInterval(UInt idInterval_);

public:
  /// planet identifier
  enum Planet {MERCURY             = 1,
               VENUS               = 2,
               EARTH               = 3,
               MARS                = 4,
               JUPITER             = 5,
               SATURN              = 6,
               URANUS              = 7,
               NEPTUNE             = 8,
               PLUTO               = 9,
               MOON                = 10,
               SUN                 = 11,
               SOLARBARYCENTER     = 12,
               EARTHMOONBARYCENTER = 13,
               NUTATION            = 14,
               LIBRATION           = 15,
               MOONFROMEARTH       = 99};

  InFileEphemerides() {}
  explicit InFileEphemerides(const FileName &name) {open(name);}     //!< Constructor.
  InFileEphemerides(const InFileEphemerides &) = delete;             //!< Disallow copy constructor
  InFileEphemerides &operator=(const InFileEphemerides &x) = delete; //!< Disallow copying
 ~InFileEphemerides() {close();}

  void open(const FileName &name);
  void close();

  /// ephemeris of @p planet relative to @p center.
  void ephemeris(const Time &timeGPS, Planet planet, Planet center, Vector3d &position, Vector3d &velocity);

  /// interpolate data (e.g. libration) to @p timeGPS (second column: derivative).
  Matrix interpolate(const Time &timeGPS, Planet planet);
};

/***** FUNCTIONS *******************************/

/** @brief Write into an ephemerides file. */
void writeFileEphemerides(const FileName &fileName, Double earthMoonRatio, const std::vector<Time> &times,
                          const std::vector<UInt> &subintervals, const std::vector<UInt> &components, const std::vector<UInt> &degree,
                          const std::vector<std::vector<std::vector<Matrix>>> &coefficients);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
