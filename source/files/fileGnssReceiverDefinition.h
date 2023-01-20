/***********************************************/
/**
* @file fileGnssReceiverDefinition.h
*
* @brief GNSS receiver definition.
*
* @author Sebastian Strasser
* @date 2019-08-28
*/
/***********************************************/

#ifndef __GROOPS_FILEGNSSRECEIVERDEFINITION__
#define __GROOPS_FILEGNSSRECEIVERDEFINITION__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GnssReceiverDefinition
static const char *docstringGnssReceiverDefinition = R"(
Contains a list of GNSS receivers which are identified by its
name, serial, and version. Defines for each receiver a list of
signal \configClass{types}{gnssType} which can be observed.
Can also be used for GNSS transmitters to define a list of
transmitted signal types. For GLONASS satellites the frequency
number can be stored in the \emph{version} field.

See \program{GnssReceiverDefinitionCreate}.

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="receiverDefinition" version="20190429">
    <receiverCount>112</receiverCount>
    <receiver>
        <name>GLONASS</name>
        <serial>R779</serial>
        <version>2</version>
        <comment/>
        <types>
            <count>4</count>
            <cell>*1CR**J</cell>
            <cell>*1PR**J</cell>
            <cell>*2CR**J</cell>
            <cell>*2PR**J</cell>
        </types>
    </receiver>
    ...
    <receiver>
        <name>GLONASS-K1</name>
        <serial>R802</serial>
        <version>7</version>
        <comment/>
        <types>
            <count>10</count>
            <cell>*1CR**O</cell>
            <cell>*1PR**O</cell>
            <cell>*2CR**O</cell>
            <cell>*2PR**O</cell>
            <cell>*3IR**</cell>
            <cell>*3QR**</cell>
            <cell>*4AR**</cell>
            <cell>*4BR**</cell>
            <cell>*6AR**</cell>
            <cell>*6BR**</cell>
        </types>
    </receiver>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSRECEIVERDEFINITION_TYPE = "receiverDefinition";

/***** TYPES ***********************************/

class GnssReceiverDefinition;
typedef std::shared_ptr<GnssReceiverDefinition> GnssReceiverDefinitionPtr;

/***** CLASS ***********************************/

/** @brief GNSS receiver definition. */
class GnssReceiverDefinition
{
  public:
  std::string name;
  std::string serial;
  std::string version;
  std::string comment;
  std::vector<GnssType> types;

  /** @brief Returns the separator between parts of the full receiver name. */
  static constexpr Char sep = '|';
  static std::string str(const std::string &name, const std::string &serial, const std::string &version) {return name+sep+serial+sep+version;}
  std::string str() const {return str(name, serial, version);}

  static GnssReceiverDefinitionPtr find(const std::vector<GnssReceiverDefinitionPtr> &receivers, const std::string &name, const std::string &serial, const std::string version);
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GnssReceiverDefinition    &x);
template<> void save(OutArchive &ar, const GnssReceiverDefinitionPtr &x);
template<> void load(InArchive  &ar, GnssReceiverDefinition    &x);
template<> void load(InArchive  &ar, GnssReceiverDefinitionPtr &x);


/** @brief Write into a GnssReceiverDefinition file. */
void writeFileGnssReceiverDefinition(const FileName &fileName, const GnssReceiverDefinition &x);

/** @brief Write into a GnssReceiverDefinition file. */
void writeFileGnssReceiverDefinition(const FileName &fileName, const std::vector<GnssReceiverDefinitionPtr> &x);


/** @brief Read from a GnssReceiverDefinition file. */
void readFileGnssReceiverDefinition(const FileName &fileName, GnssReceiverDefinition &x);

/** @brief Read from a GnssReceiverDefinition file. */
void readFileGnssReceiverDefinition(const FileName &fileName, std::vector<GnssReceiverDefinitionPtr> &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
