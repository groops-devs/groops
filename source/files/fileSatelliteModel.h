/***********************************************/
/**
* @file fileSatelliteModel.h
*
* @brief Satellite model.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-05-25
*
*/
/***********************************************/

#ifndef __GROOPS_FILESATELLITEMODEL__
#define __GROOPS_FILESATELLITEMODEL__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_SatelliteModel
static const char *docstringSatelliteModel = R"(
Properties of a satellite to model non-conservative forces (e.g. \configClass{miscAccelerations}{miscAccelerationsType}).
The file may contain surface properties, mass, drag coefficients, and antenna thrust values.

See \program{SatelliteModelCreate} and \program{SinexMetadata2SatelliteModel}.

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="satelliteModel" version="20190429">
   <satelliteCount>1</satelliteCount>
   <satellite>
       <satelliteName>GALILEO-2</satelliteName>
       <mass>7.00000000000000e+02</mass>
       <coefficientDrag>0.00000000000000e+00</coefficientDrag>
       <surfaceCount>15</surfaceCount>
       <surface>
           <type>0</type>
           <normal>
               <x>-1.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>0.00000000000000e+00</z>
           </normal>
           <area>4.40000000000000e-01</area>
           <reflexionVisible>0.00000000000000e+00</reflexionVisible>
           <diffusionVisible>7.00000000000000e-02</diffusionVisible>
           <absorptionVisible>9.30000000000000e-01</absorptionVisible>
           <reflexionInfrared>1.00000000000000e-01</reflexionInfrared>
           <diffusionInfrared>1.00000000000000e-01</diffusionInfrared>
           <absorptionInfrared>8.00000000000000e-01</absorptionInfrared>
           <hasThermalReemission>1</hasThermalReemission>
       </surface>
       ...
       <modulCount>2</modulCount>
       <modul>
           <type>1</type>
           <rotationAxis>
               <x>0.00000000000000e+00</x>
               <y>1.00000000000000e+00</y>
               <z>0.00000000000000e+00</z>
           </rotationAxis>
           <normal>
               <x>0.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>1.00000000000000e+00</z>
           </normal>
           <surface>
               <count>4</count>
               <cell>11</cell>
               <cell>12</cell>
               <cell>13</cell>
               <cell>14</cell>
           </surface>
       </modul>
       <modul>
           <type>2</type>
           <antennaThrust>
               <x>0.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>2.65000000000000e+02</z>
           </antennaThrust>
       </modul>
   </satellite>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"
#include "inputOutput/logging.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_SATELLITEMODEL_TYPE = "satelliteModel";

/***** TYPES ***********************************/

class SatelliteModel;
class SatelliteModelModule;
typedef std::shared_ptr<SatelliteModel> SatelliteModelPtr;
typedef std::shared_ptr<SatelliteModelModule> SatelliteModelModulePtr;

/***** CLASS ***********************************/

/** @brief Satellite model. */
class SatelliteModel
{
public:
  class Surface
  {
  public:
    enum Type {PLATE=0, SPHERE=1, CYLINDER=2};
    Type     type;
    Vector3d normal;
    Double   area;
    Double   absorptionVisible, absorptionInfrared;
    Double   diffusionVisible,  diffusionInfrared;
    Double   reflexionVisible,  reflexionInfrared;
    Bool     hasThermalReemission;

    Surface() : type(PLATE), area(0),
                absorptionVisible(0), absorptionInfrared(0),
                diffusionVisible(0),  diffusionInfrared(0),
                reflexionVisible(0),  reflexionInfrared(0),
                hasThermalReemission(FALSE) {}

    Double   sectionalArea(const Vector3d &direction) const;
    Vector3d accelerationDrag(const Vector3d &direction, Double velocity, Double density, Double temperature) const;
    void     accelerationPressure(const Vector3d &direction, Double visible, Double infrared, Vector3d &acc) const;
  };

  SatelliteModel(); //!< Constructor.
 ~SatelliteModel(); //!< Destructor.
  SatelliteModel &operator=(const SatelliteModel &x) = delete; //!< Disallow copying.
  SatelliteModel(const SatelliteModel &x) = delete;            //!< Disallow copy constructor.

  std::string satelliteName; //!< Name of satellite.
  Double      mass;
  Double      coefficientDrag;
  std::vector<Surface> surfaces;
  std::vector<SatelliteModelModulePtr> modules;

  /** @brief Compute new state of satellite.
  * E.g. Move parts of satellite surfaces, update mass.
  * @param time Time.
  * @param position in CRF [m].
  * @param velocity in CRF [m/s].
  * @param positionSun in CRF [m].
  * @param rotSat   Sat -> CRF
  * @param rotEarth CRF -> TRF */
  void changeState(const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);

  /** @brief Acceleration of satellite due to drag.
  * @param velocity of satellite relative to atmosphere (in satellite frame, SRF).
  * @param density of atmosphere [kg/m^3].
  * @param temperature of atmosphere [K], if zero drag is computed using simplifed coefficientDrag.
  * @return acceleration [m/s^2] in satellite frame (SRF). */
  Vector3d accelerationDrag(const Vector3d &velocity, Double density, Double temperature) const;

  /** @brief Acceleration of satellite due to solar radiation pressure/albedo.
  * @param direction unit vector of incoming radiation (in satellite frame, SRF)
  * @param visible magnitude of incident visible radiation [N/m^2]
  * @param infrared magnitude of incident infrared radiation [N/m^2]
  * @return acceleration [m/s^2] in satellite frame (SRF) */
  Vector3d accelerationPressure(const Vector3d &direction, Double visible, Double infrared) const;

  /** @brief Acceleration of satellite due to thrust. */
  Vector3d accelerationThrust() const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const SatelliteModelPtr &x);
template<> void load(InArchive  &ar, SatelliteModelPtr &x);

/** @brief Write into a SatelliteModel file. */
void writeFileSatelliteModel(const FileName &fileName, const SatelliteModelPtr &x);

/** @brief Write into a SatelliteModel file. */
void writeFileSatelliteModel(const FileName &fileName, const std::vector<SatelliteModelPtr> &x);

/** @brief Read from a SatelliteModel file. */
void readFileSatelliteModel(const FileName &fileName, SatelliteModelPtr &x);

/** @brief Read from a SatelliteModel file. */
void readFileSatelliteModel(const FileName &fileName, std::vector<SatelliteModelPtr> &x);

/// @}

/***** CLASS ***********************************/

// Internal class
class SatelliteModelModule
{
public:
  enum Type {SOLARPANEL = 1, ANTENNATHRUST = 2, MASSCHANGE = 3};

  virtual ~SatelliteModelModule() {}
  virtual Type type() const = 0;
  virtual void changeState(SatelliteModel &/*satellite*/, const Time &/*time*/, const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                           const Vector3d &/*positionSun*/, const Rotary3d &/*rotSat*/, const Rotary3d &/*rotEarth*/) {}
  virtual void accelerationDrag    (const SatelliteModel &/*satellite*/, const Vector3d &/*velocity*/,  Double /*density*/, Double /*temperature*/, Vector3d &/*a*/) const {}
  virtual void accelerationPressure(const SatelliteModel &/*satellite*/, const Vector3d &/*direction*/, Double /*visible*/, Double /*infrared*/,    Vector3d &/*a*/) const {}
  virtual void accelerationThrust  (const SatelliteModel &/*satellite*/, Vector3d &/*a*/) const {}

  /** @brief Create an epoch of given type (with new). */
  static SatelliteModelModulePtr create(Type type);

  virtual void save(OutArchive &ar) const = 0;
  virtual void load(InArchive  &ar) = 0;
};

/***********************************************/

class SatelliteModelModuleSolarPanel : public SatelliteModelModule
{
public:
  Vector3d          rotationAxis;
  Vector3d          normal;
  std::vector<UInt> indexSurface;

  Type type() const {return SOLARPANEL;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

class SatelliteModelModuleAntennaThrust : public SatelliteModelModule
{
public:
  Vector3d thrust;

  Type type() const {return ANTENNATHRUST;}
  void accelerationThrust(const SatelliteModel &satellite, Vector3d &a) const;
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

class SatelliteModelModuleMassChange : public SatelliteModelModule
{
public:
  std::vector<Time>   times;
  std::vector<Double> mass;

  Type type() const {return MASSCHANGE;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

#endif /* __GROOPS__ */
