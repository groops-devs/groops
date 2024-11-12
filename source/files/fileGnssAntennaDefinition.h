/***********************************************/
/**
* @file fileGnssAntennaDefinition.h
*
* @brief Antenna center variations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 201-11-22
*
*/
/***********************************************/

#ifndef __GROOPS_FILEGNSSANTENNADEFINITION__
#define __GROOPS_FILEGNSSANTENNADEFINITION__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GnssAntennaDefinition
static const char *docstringGnssAntennaDefinition = R"(
Contains a list of GNSS antennas which are identified by its
name (type), serial, and radome. Each antenna consists of
antenna center offsets (ACO) and antenna center variations (ACV)
for different signal \configClass{types}{gnssType} (code and phase).
The ACV values for each type are stored in an elevation and azimuth dependent grid.

%New \config{antenna} with \config{pattern}s for code (\config{type}=\verb|C**|) and phase
%(\config{type}=\verb|L**|).
%The standard deviation is expressed e.g. with \config{values}=\verb|0.001/cos(2*PI/180*zenith)|.

See also \program{GnssAntennaDefinitionCreate}, \program{GnssAntex2AntennaDefinition}.

\fig{!hb}{1.0}{fileFormatGnssAntennaDefinition}{fig:fileFormatGnssAntennaDefinition}{Antenna center variations of ASH701945D\_M for two frequencies of GPS and GLONASS}

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="antennaDefinition" version="20190429">
    <antennaCount>65</antennaCount>
    ...
    <antenna>
        <name>BLOCK IIIA</name>
        <serial>G074</serial>
        <radome>2018-109A</radome>
        <comment>PCO provided by the Aerospace Corporation, PV from estimations by ESA/CODE</comment>
        <pattern>
            <count>3</count>
            <cell>
                <type>*1*G**</type>
                <offset>
                    <x>-1.23333333333333e-03</x>
                    <y>4.33333333333333e-04</y>
                    <z>3.15200000000000e-01</z>
                </offset>
                <dZenit>1.00000000000000e+00</dZenit>
                <pattern>
                    <type>0</type>
                    <rows>1</rows>
                    <columns>18</columns>
                    <cell row="0" col="0">1.39000000000000e-02</cell>
                    <cell row="0" col="1">1.28000000000000e-02</cell>
                    <cell row="0" col="2">1.02000000000000e-02</cell>
                    <cell row="0" col="3">5.80000000000000e-03</cell>
                    <cell row="0" col="4">1.10000000000000e-03</cell>
                    <cell row="0" col="5">-4.50000000000000e-03</cell>
                    <cell row="0" col="6">-9.70000000000000e-03</cell>
                    <cell row="0" col="7">-1.28000000000000e-02</cell>
                    <cell row="0" col="8">-1.34000000000000e-02</cell>
                    <cell row="0" col="9">-1.18000000000000e-02</cell>
                    <cell row="0" col="10">-8.90000000000000e-03</cell>
                    <cell row="0" col="11">-4.50000000000000e-03</cell>
                    <cell row="0" col="12">1.20000000000000e-03</cell>
                    <cell row="0" col="13">7.20000000000000e-03</cell>
                    <cell row="0" col="14">1.33000000000000e-02</cell>
                    <cell row="0" col="15">1.33000000000000e-02</cell>
                    <cell row="0" col="16">1.33000000000000e-02</cell>
                    <cell row="0" col="17">1.33000000000000e-02</cell>
                </pattern>
            </cell>
            ...
        </pattern>
    </antenna>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"


/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSANTENNADEFINITION_TYPE    = "antennaDefinition";
constexpr UInt    FILE_GNSSANTENNADEFINITION_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class GnssAntennaDefinition;
class GnssAntennaPattern;
typedef std::shared_ptr<GnssAntennaDefinition> GnssAntennaDefinitionPtr;

/***** CLASS ***********************************/

/** @brief Antenna center variations.
* with different patterns for different signals. */
class GnssAntennaDefinition
{
  public:
  enum NoPatternFoundAction
  {
    IGNORE_OBSERVATION,
    USE_NEAREST_FREQUENCY,
    THROW_EXCEPTION
  };

  std::string name, serial;
  std::string radome;
  std::string comment;
  std::vector<GnssAntennaPattern> patterns;

  /** @brief Returns the separator between parts of the full antenna name. */
  static constexpr Char sep = '|';
  static std::string str(const std::string &name, const std::string &serial, const std::string &radome) {return name+sep+serial+sep+radome;}
  std::string str() const {return str(name, serial, radome);}

  Vector antennaVariations(Angle azimut, Angle elevation, const std::vector<GnssType> &types, NoPatternFoundAction noPatternFoundAction) const;

  /** @brief Returns id of pattern matching @p type or of nearest frequency pattern depending on @p noPatternFoundAction. Returns NULLINDEX if there are no patterns. */
  UInt findAntennaPattern(const GnssType &type, NoPatternFoundAction noPatternFoundAction) const;

  static GnssAntennaDefinitionPtr find(const std::vector<GnssAntennaDefinitionPtr> &antennaList, const std::string &name, const std::string &serial, const std::string radome);
};

/***** CLASS ***********************************/

/** @brief Antenna center variations.
* for one specfic observation type. */
class GnssAntennaPattern
{
  public:
  GnssType type;
  Vector3d offset;   // phase center relative to antenna reference point
  Angle    dZenit;
  Matrix   pattern;  // phase center variations (azimut(0..360) x zenit(0..dZenit*rows))

  // pattern estimation -> not written to file
  std::vector<std::vector<std::vector<Double>>> residuals;
  Matrix   ePe, redundancy;
  Matrix   sum, count;

  Double antennaVariations(Angle azimut, Angle elevation, Bool applyOffset=TRUE) const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GnssAntennaDefinition    &x);
template<> void save(OutArchive &ar, const GnssAntennaDefinitionPtr &x);
template<> void load(InArchive  &ar, GnssAntennaDefinition    &x);
template<> void load(InArchive  &ar, GnssAntennaDefinitionPtr &x);

template<> void save(OutArchive &ar, const GnssAntennaPattern &x);
template<> void load(InArchive  &ar, GnssAntennaPattern  &x);

/** @brief Write into a GnssAntennaDefinition file. */
void writeFileGnssAntennaDefinition(const FileName &fileName, const GnssAntennaDefinition &x);

/** @brief Write into a GnssAntennaDefinition file. */
void writeFileGnssAntennaDefinition(const FileName &fileName, const std::vector<GnssAntennaDefinitionPtr> &x);


/** @brief Read from a GnssAntennaDefinition file. */
void readFileGnssAntennaDefinition(const FileName &fileName, GnssAntennaDefinition &x);

/** @brief Read from a GnssAntennaDefinition file. */
void readFileGnssAntennaDefinition(const FileName &fileName, std::vector<GnssAntennaDefinitionPtr> &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
