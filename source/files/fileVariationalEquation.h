/***********************************************/
/**
* @file fileVariationalEquation.h
*
* @brief Variational equation arcs.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-22
*
*/
/***********************************************/

#ifndef __GROOPS_FILEVARIATIONALEQUATION__
#define __GROOPS_FILEVARIATIONALEQUATION__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_VariationalEquation
static const char *docstringVariationalEquation = R"(
Arcs with reference orbit and state transition matrices.

The file contains a reference orbit (position and velocity),
the derivatives of the orbit with respect to the satellite state vector for each arc,
transformations (rotations) between the satellite, celestial, and terrestrial frame
and a satellite macro model (see \file{SatelliteModel}{satelliteModel}).

The reference orbit can be extracted with \program{Variational2Orbit}.

See also: \program{PreprocessingVariationalEquation}.

)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_VARIATIONALEQUATION_TYPE = "variationalEquation";

/***** CLASS ***********************************/

/** @brief Variational equation arc. */
class VariationalEquationArc
{
  public:
  std::vector<Time>     times;
  Matrix                pos0, PosState;
  Matrix                vel0, VelState;
  std::vector<Rotary3d> rotSat;   //!< Sat -> CRF
  std::vector<Rotary3d> rotEarth; //!< CRF -> TRF

  OrbitArc orbitArc() const;

  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***** CLASS ***********************************/

/** @brief Variational equation arcs. */
class FileVariationalEquation
{
  FileName          _fileName;
  InFileArchive     _file;
  SatelliteModelPtr _satellite;
  UInt              _arcCount;
  UInt              _index;

public:
  /// Default Constructor.
  FileVariationalEquation() : _arcCount(0) {}

  /// Constructor.
  FileVariationalEquation(const FileName &fileName) {open(fileName);}

  /// Destructor.
  ~FileVariationalEquation() {close();}

  /** @brief Open a new file.
  * Old open file is closed before. */
  void  open(const FileName &fileName);

  /** @brief Close the file. */
  void  close();

  /** @brief Number of Arc in file. */
  UInt arcCount() const {return _arcCount;}

  /** @brief Satellite model. */
  SatelliteModelPtr satellite() const {return _satellite;}

  /** @brief Read a single Arc.
  * The operation is faster, if the arcs is read in increasing order. */
  VariationalEquationArc readArc(UInt arcNo);
};

/***** FUNCTIONS *******************************/

/** @brief Write a variational equation file. */
void writeFileVariationalEquation(const FileName &fileName, SatelliteModelPtr satellite, VariationalEquationArc arc);

/** @brief Write a variational equation file. */
void writeFileVariationalEquation(const FileName &fileName, SatelliteModelPtr satellite, std::vector<VariationalEquationArc> arcs);

/***********************************************/

/// @}

#endif
