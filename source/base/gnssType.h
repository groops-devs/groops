/***********************************************/
/**
* @file gnssType.h
*
* @brief Defines a GNSS observation type according to the RINEX 3 definition.
*
* Due to the strict weak ordering requirement, wildcard matching is not supported for e.g. sets or maps.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTYPE__
#define __GROOPS_GNSSTYPE__

// Latex documentation
#ifdef DOCSTRING_GnssType
static const char *docstringGnssType = R"(
\section{GnssType}\label{gnssType}
A GnssType string consists of six parts (type, frequency, attribute, system, PRN, frequency number)
represented by seven characters.
\begin{itemize}
\item The first three characters (representing type, frequency, and attribute) correspond to the observation codes of the
      \href{https://files.igs.org/pub/data/format/rinex305.pdf}{RINEX 3 definition}.
\item The satellite system character also follows the RINEX 3 definition:
      \begin{itemize}
        \item \verb|G| = GPS
        \item \verb|R| = GLONASS
        \item \verb|E| = Galileo
        \item \verb|C| = BeiDou
        \item \verb|S| = SBAS
        \item \verb|J| = QZSS
        \item \verb|I| = IRNSS
      \end{itemize}
\item PRN is a two-digit number identifying a satellite.
\item Frequency number is only used for GLONASS, where the range -7 to 14 is represented by letters starting with A.
\end{itemize}

Each part of a GnssType string can be replaced by a wildcard '\verb|*|', enabling the use of these strings as patterns,
for example to select a subset of observations (e.g. \verb|C**G**| matches all GPS code/range observations).
Trailing wildcards are optional, meaning \verb|L1*R| is automatically expanded to \verb|L1*R***|.
For some RINEX 2 types (e.g. Galileo L5) the RINEX 3 attribute is unknown/undefined and can be replaced by \verb|?|,
for example \verb|L5?E01|.

Examples:
\begin{itemize}
\item \verb|C1CG23| = code/range observation, L1 frequency, derived from C/A code, GPS, PRN 23
\item \verb|L2PR05B| = phase observation, G2 frequency, derived from P code, GLONASS, PRN 05, frequency number -6
\item \verb|*5*E**| = all observation types, E5a frequency, all attributes, Galileo, all PRNs
\end{itemize}
)";
#endif

/***********************************************/

#include "base/importStd.h"

/***** CLASS ***********************************/

/**
* @brief Defines a GNSS observation type according to the RINEX 3 definition.
* @ingroup base
*
* Due to the strict weak ordering requirement, wildcard matching is not supported for e.g. sets or maps. */
class GnssType
{
public:
  /// Bit masks
  static const GnssType PRN;        ///< satellite identification number (PRN, SBAS-100, ...)
  static const GnssType SYSTEM;     ///< satellite system (GPS, GLONASS, ...)
  static const GnssType FREQUENCY;  ///< frequency
  static const GnssType TYPE;       ///< phase, range, doppler, ...
  static const GnssType ATTRIBUTE;  ///< attribute
  static const GnssType FREQ_NO;    ///< GLONASS frequency number
  static const GnssType ALL;        ///< = TYPE + FREQUENCY + ATTRIBUTE + SYSTEM + PRN + GLO_FREQ_NO
  static const GnssType NOPRN;      ///< = TYPE + FREQUENCY + ATTRIBUTE + SYSTEM       + GLO_FREQ_NO
  /// system
  static const GnssType GPS;
  static const GnssType GLONASS;
  static const GnssType GALILEO;
  static const GnssType BDS;
  static const GnssType SBAS;
  static const GnssType QZSS;
  static const GnssType IRNSS;
  /// frequency
  static const GnssType L1;   ///< GPS
  static const GnssType L2;
  static const GnssType L5;
  static const GnssType E1;   ///< GALILEO
  static const GnssType E5a;
  static const GnssType E5b;
  static const GnssType E5;
  static const GnssType E6;
  static const GnssType G1;   ///< GLONASS
  static const GnssType G1a;
  static const GnssType G2;
  static const GnssType G2a;
  static const GnssType G3;
  static const GnssType B1;   ///< BDS
  static const GnssType B1C;
  static const GnssType B2a;
  static const GnssType B2b;
  static const GnssType B2;
  static const GnssType B3;
  static const GnssType L6;   ///< QZSS
  static const GnssType S9;   ///< IRNSS S
  /// observation type
  static const GnssType RANGE;
  static const GnssType PHASE;
  static const GnssType DOPPLER;
  static const GnssType SNR;
  static const GnssType IONODELAY;
  static const GnssType CHANNEL;
  static const GnssType AZIMUT;     ///< FREQUENCY L1: receiver, L2: transmitter
  static const GnssType ELEVATION;  ///< FREQUENCY L1: receiver, L2: transmitter
  static const GnssType ROTI;       ///< Rate of Tec Index
  static const GnssType IONOINDEX;  ///< Ionospheric index (sigma_phi)
  /// attributes
  static const GnssType C;          ///< C/A - Code
  static const GnssType W;          ///< P Z-tracking and similar (AS on)
  static const GnssType D;          ///< L1(C/A)+(P2-P1) (semi-codeless)
  static const GnssType X;
  static const GnssType S;
  static const GnssType L;
  static const GnssType I;
  static const GnssType Q;
  static const GnssType A;
  static const GnssType B;
  static const GnssType Z;
  static const GnssType P;         ///< military P (AS off) - Code
  static const GnssType Y;         ///< military
  static const GnssType M;         ///< military
  static const GnssType E;
  static const GnssType UNKNOWN_ATTRIBUTE;

  static const GnssType C1CG;
  static const GnssType C1SG;
  static const GnssType C1LG;
  static const GnssType C1XG;
  static const GnssType C1WG;
  static const GnssType C2CG;
  static const GnssType C2DG;
  static const GnssType C2SG;
  static const GnssType C2LG;
  static const GnssType C2XG;
  static const GnssType C2WG;
  static const GnssType C2UG;      ///< unknown attribute (RINEX 2)
  static const GnssType C5IG;
  static const GnssType C5QG;
  static const GnssType C5XG;
  static const GnssType C5UG;      ///< unknown attribute (RINEX 2)
  static const GnssType L1_G;
  static const GnssType L2_G;
  static const GnssType L5_G;

  static const GnssType C1CR;
  static const GnssType C1PR;
  static const GnssType C4AR;
  static const GnssType C4BR;
  static const GnssType C4XR;
  static const GnssType C2CR;
  static const GnssType C2PR;
  static const GnssType C6AR;
  static const GnssType C6BR;
  static const GnssType C6XR;
  static const GnssType C3IR;
  static const GnssType C3QR;
  static const GnssType C3XR;
  static const GnssType L1_R;
  static const GnssType L4_R;
  static const GnssType L2_R;
  static const GnssType L6_R;
  static const GnssType L3_R;

  static const GnssType C1AE;
  static const GnssType C1BE;
  static const GnssType C1CE;
  static const GnssType C1XE;
  static const GnssType C1ZE;
  static const GnssType C1UE;      ///< unknown attribute (RINEX 2)
  static const GnssType C5IE;
  static const GnssType C5QE;
  static const GnssType C5XE;
  static const GnssType C5UE;      ///< unknown attribute (RINEX 2)
  static const GnssType C7IE;
  static const GnssType C7QE;
  static const GnssType C7XE;
  static const GnssType C7UE;      ///< unknown attribute (RINEX 2)
  static const GnssType C8IE;
  static const GnssType C8QE;
  static const GnssType C8XE;
  static const GnssType C8UE;      ///< unknown attribute (RINEX 2)
  static const GnssType C6AE;
  static const GnssType C6BE;
  static const GnssType C6CE;
  static const GnssType C6XE;
  static const GnssType C6ZE;
  static const GnssType C6UE;      ///< unknown attribute (RINEX 2)
  static const GnssType L1_E;
  static const GnssType L5_E;
  static const GnssType L7_E;
  static const GnssType L8_E;
  static const GnssType L6_E;

  static const GnssType C2IC;
  static const GnssType C2QC;
  static const GnssType C2XC;
  static const GnssType C1DC;
  static const GnssType C1PC;
  static const GnssType C1XC;
  static const GnssType C1SC;
  static const GnssType C1LC;
  static const GnssType C1ZC;
  static const GnssType C5DC;
  static const GnssType C5PC;
  static const GnssType C5XC;
  static const GnssType C7IC;
  static const GnssType C7QC;
  static const GnssType C7XC;
  static const GnssType C7DC;
  static const GnssType C7PC;
  static const GnssType C7ZC;
  static const GnssType C8DC;
  static const GnssType C8PC;
  static const GnssType C8XC;
  static const GnssType C6IC;
  static const GnssType C6QC;
  static const GnssType C6XC;
  static const GnssType C6DC;
  static const GnssType C6PC;
  static const GnssType C6ZC;
  static const GnssType L2_C;
  static const GnssType L1_C;
  static const GnssType L5_C;
  static const GnssType L7_C;
  static const GnssType L8_C;
  static const GnssType L6_C;

  static const GnssType C1CJ;
  static const GnssType C1SJ;
  static const GnssType C1LJ;
  static const GnssType C1XJ;
  static const GnssType C1ZJ;
  static const GnssType C1BJ;
  static const GnssType C2SJ;
  static const GnssType C2LJ;
  static const GnssType C2XJ;
  static const GnssType C5IJ;
  static const GnssType C5QJ;
  static const GnssType C5XJ;
  static const GnssType C5DJ;
  static const GnssType C5PJ;
  static const GnssType C5ZJ;
  static const GnssType C6SJ;
  static const GnssType C6LJ;
  static const GnssType C6XJ;
  static const GnssType C6EJ;
  static const GnssType C6ZJ;
  static const GnssType L1_J;
  static const GnssType L2_J;
  static const GnssType L5_J;
  static const GnssType L6_J;

  UInt64 type;

  constexpr GnssType() : type(0) {}
  constexpr explicit GnssType(UInt64 t) : type(t) {}
  constexpr GnssType(const GnssType &t) : type(t.type) {}
  explicit GnssType(const std::string &str);
  GnssType &operator=(const GnssType &t) {type = t.type; return *this;}

  Double      frequency() const;
  Double      wavelength() const;
  std::string str() const;
  std::string prnStr() const {return str().substr(3,3);}
  UInt        prn() const    {return type & PRN.type;}
  Int         frequencyNumber() const; ///< GLONASS frequency number (9999 if not set).
  void        setFrequencyNumber(Int number);

  /** @brief First order STEC influence [m/TECU]. */
  Double ionosphericFactor() const;

  /** @brief Returns true if both vectors are of the same size and contain only the same types, independent of sorting. */
  static Bool allEqual(const std::vector<GnssType> &types1, const std::vector<GnssType> &types2, GnssType mask=GnssType::ALL);

  /** @brief Returns the index of types vector. If not found NULLINDEX is returned. */
  static UInt index(const std::vector<GnssType> &types, GnssType type);

  /** @brief Returns true if a wilcard (*) is used in the parts given by @a mask. */
  Bool hasWildcard(GnssType mask=GnssType::ALL) const;

  GnssType &operator+=(const GnssType &t);
  GnssType &operator&=(const GnssType &t);
  GnssType  operator+ (const GnssType &t) const {GnssType t1(*this); t1+=t; return t1;}
  GnssType  operator& (const GnssType &t) const {GnssType t1(*this); t1&=t; return t1;}
  GnssType  operator~ () const;

  /** @brief Returns true if all parts are either the same or match a wildcard. */
  Bool operator==(const GnssType &t) const;

  /** @brief Returns true if @a this is smaller than @a other or @a other is a wildcard. When used for sorting, wildcards are at the end. */
  Bool operator< (const GnssType &t) const;

  /** @brief Returns true if any part is neither the same nor matches a wildcard. */
  Bool operator!=(const GnssType &t) const { return !(*this==t); }
};

/***********************************************/

#endif /* __GROOPS___ */
