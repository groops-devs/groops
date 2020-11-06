/***********************************************/
/**
* @file gnssType.cpp
*
* @brief GNSS observation types.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-30
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/gnssType.h"

/***********************************************/

// Bit masks
const GnssType GnssType::PRN        = GnssType(static_cast<UInt64>(0xff));       // satellite identification number (PRN, SBAS-100, ...)
const GnssType GnssType::SYSTEM     = GnssType(static_cast<UInt64>(0x0f) << 8);  // satellite system (GPS, GLONASS, ...)
const GnssType GnssType::FREQUENCY  = GnssType(static_cast<UInt64>(0x0f) << 12); // frequency
const GnssType GnssType::TYPE       = GnssType(static_cast<UInt64>(0xff) << 16); // phase, range, doppler, ...
const GnssType GnssType::ATTRIBUTE  = GnssType(static_cast<UInt64>(0x7f) << 24); // attribute
const GnssType GnssType::FREQ_NO    = GnssType(static_cast<UInt64>(0xff) << 32); // GLONASS frequency number
const GnssType GnssType::ALL        = GnssType::TYPE + GnssType::FREQUENCY + GnssType::ATTRIBUTE + GnssType::SYSTEM + GnssType::PRN + GnssType::FREQ_NO;
const GnssType GnssType::NOPRN      = GnssType::TYPE + GnssType::FREQUENCY + GnssType::ATTRIBUTE + GnssType::SYSTEM                 + GnssType::FREQ_NO;
// system
const GnssType GnssType::GPS       = GnssType(static_cast<UInt64>(1) << 8);
const GnssType GnssType::GLONASS   = GnssType(static_cast<UInt64>(2) << 8);
const GnssType GnssType::SBAS      = GnssType(static_cast<UInt64>(3) << 8);
const GnssType GnssType::BDS       = GnssType(static_cast<UInt64>(4) << 8);
const GnssType GnssType::GALILEO   = GnssType(static_cast<UInt64>(5) << 8);
const GnssType GnssType::QZSS      = GnssType(static_cast<UInt64>(6) << 8);
const GnssType GnssType::IRNSS     = GnssType(static_cast<UInt64>(7) << 8);
// frequency
const GnssType GnssType::L1        = GnssType(static_cast<UInt64>(1) << 12); // GPS
const GnssType GnssType::L2        = GnssType(static_cast<UInt64>(2) << 12);
const GnssType GnssType::L5        = GnssType(static_cast<UInt64>(5) << 12);
const GnssType GnssType::E1        = GnssType(static_cast<UInt64>(1) << 12); // GALILEO
const GnssType GnssType::E5a       = GnssType(static_cast<UInt64>(5) << 12);
const GnssType GnssType::E5b       = GnssType(static_cast<UInt64>(7) << 12);
const GnssType GnssType::E5        = GnssType(static_cast<UInt64>(8) << 12);
const GnssType GnssType::E6        = GnssType(static_cast<UInt64>(6) << 12);
const GnssType GnssType::G1        = GnssType(static_cast<UInt64>(1) << 12); // GLONASS
const GnssType GnssType::G1a       = GnssType(static_cast<UInt64>(4) << 12);
const GnssType GnssType::G2        = GnssType(static_cast<UInt64>(2) << 12);
const GnssType GnssType::G2a       = GnssType(static_cast<UInt64>(6) << 12);
const GnssType GnssType::G3        = GnssType(static_cast<UInt64>(3) << 12);
const GnssType GnssType::B1        = GnssType(static_cast<UInt64>(1) << 12); // BDS
const GnssType GnssType::B1_2      = GnssType(static_cast<UInt64>(2) << 12);
const GnssType GnssType::B2        = GnssType(static_cast<UInt64>(8) << 12);
const GnssType GnssType::B2a       = GnssType(static_cast<UInt64>(5) << 12);
const GnssType GnssType::B2b       = GnssType(static_cast<UInt64>(7) << 12);
const GnssType GnssType::B3        = GnssType(static_cast<UInt64>(6) << 12);
const GnssType GnssType::L6        = GnssType(static_cast<UInt64>(6) << 12); // QZSS
const GnssType GnssType::S9        = GnssType(static_cast<UInt64>(9) << 12); // IRNSS S
// observation type
const GnssType GnssType::RANGE     = GnssType(static_cast<UInt64>(1) << 16);
const GnssType GnssType::PHASE     = GnssType(static_cast<UInt64>(2) << 16);
const GnssType GnssType::DOPPLER   = GnssType(static_cast<UInt64>(3) << 16);
const GnssType GnssType::SNR       = GnssType(static_cast<UInt64>(4) << 16);
const GnssType GnssType::IONODELAY = GnssType(static_cast<UInt64>(5) << 16);
const GnssType GnssType::AZIMUT    = GnssType(static_cast<UInt64>(6) << 16);
const GnssType GnssType::ELEVATION = GnssType(static_cast<UInt64>(7) << 16);
const GnssType GnssType::ROTI      = GnssType(static_cast<UInt64>(8) << 16);  // Rate of Tec Index
const GnssType GnssType::IONOINDEX = GnssType(static_cast<UInt64>(9) << 16);  // Ionospheric index (sigma_phi)
const GnssType GnssType::CHANNEL   = GnssType(static_cast<UInt64>(10)<< 16);
// attributes
const GnssType GnssType::C         = GnssType(static_cast<UInt64>(1) << 24);  // C/A - Code
const GnssType GnssType::W         = GnssType(static_cast<UInt64>(2) << 24);  // P Z-tracking and similar (AS on)
const GnssType GnssType::D         = GnssType(static_cast<UInt64>(3) << 24);  // L1(C/A)+(P2-P1) (semi-codeless)
const GnssType GnssType::X         = GnssType(static_cast<UInt64>(4) << 24);
const GnssType GnssType::S         = GnssType(static_cast<UInt64>(5) << 24);
const GnssType GnssType::L         = GnssType(static_cast<UInt64>(6) << 24);
const GnssType GnssType::I         = GnssType(static_cast<UInt64>(7) << 24);
const GnssType GnssType::Q         = GnssType(static_cast<UInt64>(8) << 24);
const GnssType GnssType::A         = GnssType(static_cast<UInt64>(9) << 24);
const GnssType GnssType::B         = GnssType(static_cast<UInt64>(10) << 24);
const GnssType GnssType::Z         = GnssType(static_cast<UInt64>(11) << 24);
const GnssType GnssType::P         = GnssType(static_cast<UInt64>(12) << 24);
const GnssType GnssType::Y         = GnssType(static_cast<UInt64>(13) << 24);
const GnssType GnssType::M         = GnssType(static_cast<UInt64>(14) << 24);
const GnssType GnssType::UNKNOWN_ATTRIBUTE = GnssType(static_cast<UInt64>(99) << 24);

// some codes
const GnssType GnssType::C1CG      = GnssType::RANGE + GnssType::L1  + GnssType::C + GnssType::GPS;
const GnssType GnssType::C1SG      = GnssType::RANGE + GnssType::L1  + GnssType::S + GnssType::GPS;
const GnssType GnssType::C1LG      = GnssType::RANGE + GnssType::L1  + GnssType::L + GnssType::GPS;
const GnssType GnssType::C1XG      = GnssType::RANGE + GnssType::L1  + GnssType::X + GnssType::GPS;
const GnssType GnssType::C1WG      = GnssType::RANGE + GnssType::L1  + GnssType::W + GnssType::GPS;
const GnssType GnssType::C2CG      = GnssType::RANGE + GnssType::L2  + GnssType::C + GnssType::GPS;
const GnssType GnssType::C2DG      = GnssType::RANGE + GnssType::L2  + GnssType::D + GnssType::GPS;
const GnssType GnssType::C2SG      = GnssType::RANGE + GnssType::L2  + GnssType::S + GnssType::GPS;
const GnssType GnssType::C2LG      = GnssType::RANGE + GnssType::L2  + GnssType::L + GnssType::GPS;
const GnssType GnssType::C2XG      = GnssType::RANGE + GnssType::L2  + GnssType::X + GnssType::GPS;
const GnssType GnssType::C2WG      = GnssType::RANGE + GnssType::L2  + GnssType::W + GnssType::GPS;
const GnssType GnssType::C2UG      = GnssType::RANGE + GnssType::L2  + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GPS;
const GnssType GnssType::C5IG      = GnssType::RANGE + GnssType::L5  + GnssType::I + GnssType::GPS;
const GnssType GnssType::C5QG      = GnssType::RANGE + GnssType::L5  + GnssType::Q + GnssType::GPS;
const GnssType GnssType::C5XG      = GnssType::RANGE + GnssType::L5  + GnssType::X + GnssType::GPS;
const GnssType GnssType::C5UG      = GnssType::RANGE + GnssType::L5  + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GPS;
const GnssType GnssType::L1_G      = GnssType::PHASE + GnssType::L1                + GnssType::GPS;
const GnssType GnssType::L2_G      = GnssType::PHASE + GnssType::L2                + GnssType::GPS;
const GnssType GnssType::L5_G      = GnssType::PHASE + GnssType::L5                + GnssType::GPS;

const GnssType GnssType::C1CR      = GnssType::RANGE + GnssType::G1  + GnssType::C + GnssType::GLONASS;
const GnssType GnssType::C1PR      = GnssType::RANGE + GnssType::G1  + GnssType::P + GnssType::GLONASS;
const GnssType GnssType::C4AR      = GnssType::RANGE + GnssType::G1a + GnssType::A + GnssType::GLONASS;
const GnssType GnssType::C4BR      = GnssType::RANGE + GnssType::G1a + GnssType::B + GnssType::GLONASS;
const GnssType GnssType::C4XR      = GnssType::RANGE + GnssType::G1a + GnssType::X + GnssType::GLONASS;
const GnssType GnssType::C2CR      = GnssType::RANGE + GnssType::G2  + GnssType::C + GnssType::GLONASS;
const GnssType GnssType::C2PR      = GnssType::RANGE + GnssType::G2  + GnssType::P + GnssType::GLONASS;
const GnssType GnssType::C6AR      = GnssType::RANGE + GnssType::G2a + GnssType::A + GnssType::GLONASS;
const GnssType GnssType::C6BR      = GnssType::RANGE + GnssType::G2a + GnssType::B + GnssType::GLONASS;
const GnssType GnssType::C6XR      = GnssType::RANGE + GnssType::G2a + GnssType::X + GnssType::GLONASS;
const GnssType GnssType::C3IR      = GnssType::RANGE + GnssType::G3  + GnssType::I + GnssType::GLONASS;
const GnssType GnssType::C3QR      = GnssType::RANGE + GnssType::G3  + GnssType::Q + GnssType::GLONASS;
const GnssType GnssType::C3XR      = GnssType::RANGE + GnssType::G3  + GnssType::X + GnssType::GLONASS;
const GnssType GnssType::L1_R      = GnssType::PHASE + GnssType::G1                + GnssType::GLONASS;
const GnssType GnssType::L4_R      = GnssType::PHASE + GnssType::G1a               + GnssType::GLONASS;
const GnssType GnssType::L2_R      = GnssType::PHASE + GnssType::G2                + GnssType::GLONASS;
const GnssType GnssType::L6_R      = GnssType::PHASE + GnssType::G2a               + GnssType::GLONASS;
const GnssType GnssType::L3_R      = GnssType::PHASE + GnssType::G3                + GnssType::GLONASS;

const GnssType GnssType::C1AE      = GnssType::RANGE + GnssType::E1  + GnssType::A + GnssType::GALILEO;
const GnssType GnssType::C1BE      = GnssType::RANGE + GnssType::E1  + GnssType::B + GnssType::GALILEO;
const GnssType GnssType::C1CE      = GnssType::RANGE + GnssType::E1  + GnssType::C + GnssType::GALILEO;
const GnssType GnssType::C1XE      = GnssType::RANGE + GnssType::E1  + GnssType::X + GnssType::GALILEO;
const GnssType GnssType::C1ZE      = GnssType::RANGE + GnssType::E1  + GnssType::Z + GnssType::GALILEO;
const GnssType GnssType::C1UE      = GnssType::RANGE + GnssType::E1  + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GALILEO;
const GnssType GnssType::C5IE      = GnssType::RANGE + GnssType::E5a + GnssType::I + GnssType::GALILEO;
const GnssType GnssType::C5QE      = GnssType::RANGE + GnssType::E5a + GnssType::Q + GnssType::GALILEO;
const GnssType GnssType::C5XE      = GnssType::RANGE + GnssType::E5a + GnssType::X + GnssType::GALILEO;
const GnssType GnssType::C5UE      = GnssType::RANGE + GnssType::E5a + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GALILEO;
const GnssType GnssType::C7IE      = GnssType::RANGE + GnssType::E5b + GnssType::I + GnssType::GALILEO;
const GnssType GnssType::C7QE      = GnssType::RANGE + GnssType::E5b + GnssType::Q + GnssType::GALILEO;
const GnssType GnssType::C7XE      = GnssType::RANGE + GnssType::E5b + GnssType::X + GnssType::GALILEO;
const GnssType GnssType::C7UE      = GnssType::RANGE + GnssType::E5b + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GALILEO;
const GnssType GnssType::C8IE      = GnssType::RANGE + GnssType::E5  + GnssType::I + GnssType::GALILEO;
const GnssType GnssType::C8QE      = GnssType::RANGE + GnssType::E5  + GnssType::Q + GnssType::GALILEO;
const GnssType GnssType::C8XE      = GnssType::RANGE + GnssType::E5  + GnssType::X + GnssType::GALILEO;
const GnssType GnssType::C8UE      = GnssType::RANGE + GnssType::E5  + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GALILEO;
const GnssType GnssType::C6AE      = GnssType::RANGE + GnssType::E6  + GnssType::A + GnssType::GALILEO;
const GnssType GnssType::C6BE      = GnssType::RANGE + GnssType::E6  + GnssType::B + GnssType::GALILEO;
const GnssType GnssType::C6CE      = GnssType::RANGE + GnssType::E6  + GnssType::C + GnssType::GALILEO;
const GnssType GnssType::C6XE      = GnssType::RANGE + GnssType::E6  + GnssType::X + GnssType::GALILEO;
const GnssType GnssType::C6ZE      = GnssType::RANGE + GnssType::E6  + GnssType::Z + GnssType::GALILEO;
const GnssType GnssType::C6UE      = GnssType::RANGE + GnssType::E6  + GnssType::UNKNOWN_ATTRIBUTE + GnssType::GALILEO;
const GnssType GnssType::L1_E      = GnssType::PHASE + GnssType::E1                + GnssType::GALILEO;
const GnssType GnssType::L5_E      = GnssType::PHASE + GnssType::E5a               + GnssType::GALILEO;
const GnssType GnssType::L7_E      = GnssType::PHASE + GnssType::E5b               + GnssType::GALILEO;
const GnssType GnssType::L8_E      = GnssType::PHASE + GnssType::E5                + GnssType::GALILEO;
const GnssType GnssType::L6_E      = GnssType::PHASE + GnssType::E6                + GnssType::GALILEO;

/***********************************************/

GnssType::GnssType(const std::string &str)
{
  try
  {
    type = 0;
    if(str.size()<=0)
      return;

    switch(str.at(0))
    {
      case 'C': type += RANGE.type;     break;
      case 'L': type += PHASE.type;     break;
      case 'D': type += DOPPLER.type;   break;
      case 'S': type += SNR.type;       break;
      case 'A': type += AZIMUT.type;    break;
      case 'E': type += ELEVATION.type; break;
      case 'I': type += IONODELAY.type; break;
      case 'X': type += CHANNEL.type;   break;
      case 'R': type += ROTI.type;      break;
      case 'P': type += IONOINDEX.type; break;
      case '*': break;
      default:
        throw(Exception("Unknown GnssType string: "+str));
    }
    if(str.size()<=1)
      return;

    switch(str.at(1))
    {
      case '1': type += L1.type;   break;
      case '2': type += L2.type;   break;
      case '3': type += G3.type;   break;
      case '4': type += G1a.type;  break;
      case '5': type += L5.type;   break;
      case '6': type += E6.type;   break;
      case '7': type += E5b.type;  break;
      case '8': type += E5.type;   break;
      case '9': type += S9.type;   break;
      case '*': break;
      default:
        throw(Exception("Unknown GnssType string: "+str));
    }
    if(str.size()<=2)
      return;

    switch(str.at(2))
    {
      case 'C': type += C.type; break;
      case 'W': type += W.type; break;
      case 'D': type += D.type; break;
      case 'X': type += X.type; break;
      case 'S': type += S.type; break;
      case 'L': type += L.type; break;
      case 'I': type += I.type; break;
      case 'Q': type += Q.type; break;
      case 'A': type += A.type; break;
      case 'B': type += B.type; break;
      case 'Z': type += Z.type; break;
      case 'P': type += P.type; break;
      case 'Y': type += Y.type; break;
      case 'M': type += M.type; break;
      case '?': type += UNKNOWN_ATTRIBUTE.type; break;
      case '*': break;
      default:
        throw(Exception("Unknown GnssType string: "+str));
    }
    if(str.size()<=3)
      return;

    switch(str.at(3))
    {
      case 'G': type += GPS.type;     break;
      case 'R': type += GLONASS.type; break;
      case 'E': type += GALILEO.type; break;
      case 'C': type += BDS.type;     break;
      case 'S': type += SBAS.type;    break;
      case 'J': type += QZSS.type;    break;
      case 'I': type += IRNSS.type;   break;
      case '*': break;
      default:
        throw(Exception("Unknown GnssType string: "+str));
    }
    if(str.size()<=4)
      return;

    if(str.at(4)!='*')
    {
      UInt prn;
      std::stringstream ss(str.substr(4,2));
      ss>>prn;
      type += prn;
    }

    // GLONASS frequency number
    if(str.size()>6 && 'A'<=str.at(6) && str.at(6)<='Z')
      setFrequencyNumber(str.at(6)-'A'-7);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssType::str() const
{
  try
  {
    std::stringstream ss;

    if     ((type & TYPE.type) == RANGE.type)     ss<<'C';
    else if((type & TYPE.type) == PHASE.type)     ss<<'L';
    else if((type & TYPE.type) == DOPPLER.type)   ss<<'D';
    else if((type & TYPE.type) == SNR.type)       ss<<'S';
    else if((type & TYPE.type) == AZIMUT.type)    ss<<'A';
    else if((type & TYPE.type) == ELEVATION.type) ss<<'E';
    else if((type & TYPE.type) == IONODELAY.type) ss<<'I';
    else if((type & TYPE.type) == CHANNEL.type)   ss<<'X';
    else if((type & TYPE.type) == ROTI.type)      ss<<'R';
    else if((type & TYPE.type) == IONOINDEX.type) ss<<'P';
    else if((type & TYPE.type) == 0)              ss<<'*';
    else ss<<'?';


    if     ((type & FREQUENCY.type) == L1.type)   ss<<'1';
    else if((type & FREQUENCY.type) == L2.type)   ss<<'2';
    else if((type & FREQUENCY.type) == L5.type)   ss<<'5';
    else if((type & FREQUENCY.type) == E5b.type)  ss<<'7';
    else if((type & FREQUENCY.type) == E5.type)   ss<<'8';
    else if((type & FREQUENCY.type) == E6.type)   ss<<'6';
    else if((type & FREQUENCY.type) == G3.type)   ss<<'3';
    else if((type & FREQUENCY.type) == G1a.type)  ss<<'4';
    else if((type & FREQUENCY.type) == S9.type)   ss<<'9';
    else if((type & FREQUENCY.type) == 0)         ss<<'*';
    else ss<<'?';

    if     ((type & ATTRIBUTE.type) == C.type) ss<<'C';
    else if((type & ATTRIBUTE.type) == W.type) ss<<'W';
    else if((type & ATTRIBUTE.type) == D.type) ss<<'D';
    else if((type & ATTRIBUTE.type) == X.type) ss<<'X';
    else if((type & ATTRIBUTE.type) == S.type) ss<<'S';
    else if((type & ATTRIBUTE.type) == L.type) ss<<'L';
    else if((type & ATTRIBUTE.type) == I.type) ss<<'I';
    else if((type & ATTRIBUTE.type) == Q.type) ss<<'Q';
    else if((type & ATTRIBUTE.type) == A.type) ss<<'A';
    else if((type & ATTRIBUTE.type) == B.type) ss<<'B';
    else if((type & ATTRIBUTE.type) == Z.type) ss<<'Z';
    else if((type & ATTRIBUTE.type) == P.type) ss<<'P';
    else if((type & ATTRIBUTE.type) == Y.type) ss<<'Y';
    else if((type & ATTRIBUTE.type) == M.type) ss<<'M';
    else if((type & ATTRIBUTE.type) == UNKNOWN_ATTRIBUTE.type) ss<<'?';
    else if((type & ATTRIBUTE.type) == 0)      ss<<'*';
    else ss<<'?';

    if     ((type & SYSTEM.type) == GPS.type)     ss<<'G';
    else if((type & SYSTEM.type) == GLONASS.type) ss<<'R';
    else if((type & SYSTEM.type) == GALILEO.type) ss<<'E';
    else if((type & SYSTEM.type) == BDS.type)     ss<<'C';
    else if((type & SYSTEM.type) == SBAS.type)    ss<<'S';
    else if((type & SYSTEM.type) == QZSS.type)    ss<<'J';
    else if((type & SYSTEM.type) == IRNSS.type)   ss<<'I';
    else if((type & SYSTEM.type) == 0)            ss<<'*';
    else ss<<'?';

    if((type & PRN.type) == 0)
      ss<<"**";
    else
      ss<<prn()%"%02i"s;

    if(frequencyNumber() != 9999)
      ss<<static_cast<Char>(frequencyNumber()+7+'A');
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Double GnssType::frequency() const
{
  try
  {
    if((*this==GPS) || (*this==SBAS) || (*this==GALILEO) || (*this==QZSS) || (*this==IRNSS))
    {
      if     (*this == L1)  return 1575.42e6;
      else if(*this == L2)  return 1227.60e6;
      else if(*this == L5)  return 1176.45e6;
      else if(*this == E5b) return 1207.140e6;
      else if(*this == E5)  return 1191.795e6;
      else if(*this == E6)  return 1278.75e6;
      else if(*this == S9)  return 2492.028e6;
    }
    else if(*this==GLONASS)
    {
      if((*this == G1 || *this == G2) && frequencyNumber() == 9999)
        throw(Exception(str()+": GLONASS frequency number not set"));
      if     (*this == G1)  return 1602e6 + frequencyNumber()*9e6/16;
      else if(*this == G1a) return 1600.995e6;
      else if(*this == G2)  return 1246e6 + frequencyNumber()*7e6/16;
      else if(*this == G2a) return 1248.06e6;
      else if(*this == G3)  return 1202.025e6;
    }
    else if(*this==BDS)
    {
      if     (*this == B1)   return 1575.42e6;
      else if(*this == B1_2) return 1561.098e6;
      else if(*this == B2)   return 1191.795e6;
      else if(*this == B2a)  return 1176.45e6;
      else if(*this == B2b)  return 1207.14e6;
      else if(*this == B3)   return 1268.52e6;
    }
    throw(Exception(str()+": frequency unknown"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssType::wavelength() const
{
  return LIGHT_VELOCITY/frequency();
}

/***********************************************/

Int GnssType::frequencyNumber() const
{
  if(!(type & FREQ_NO.type))
    return 9999;
  return static_cast<Int>((type & FREQ_NO.type)>>32)-8;
}

/***********************************************/

void GnssType::setFrequencyNumber(Int number)
{
  type &= ~FREQ_NO.type;
  if((-7 <= number) && (number <= 255-8))
    type += (static_cast<UInt64>(number+8)<<32) & FREQ_NO.type;
}

/***********************************************/

Double GnssType::ionosphericFactor() const
{
  try
  {
    if(*this == GnssType::PHASE)
      return -Ionosphere::Ap/std::pow(frequency(), 2);
    if(*this == GnssType::RANGE)
      return +Ionosphere::Ap/std::pow(frequency(), 2);
    throw(Exception("only defined for PHASE and RANGE"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssType::allEqual(const std::vector<GnssType> &types1, const std::vector<GnssType> &types2, GnssType mask)
{
  try
  {
    if(types1.size() != types2.size())
      return FALSE;
    for(UInt i=0; i<types1.size(); i++)
      if(index(types2, types1.at(i) & mask) == NULLINDEX)
        return FALSE;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssType::index(const std::vector<GnssType> &types, GnssType type)
{
  for(UInt i=0; i<types.size(); i++)
    if(type == types.at(i))
      return i;
  return NULLINDEX;
}

/***********************************************/

Bool GnssType::hasWildcard(GnssType mask) const
{
  return ((mask.type & TYPE.type)      && !(type & TYPE.type))      ||
         ((mask.type & FREQUENCY.type) && !(type & FREQUENCY.type)) ||
         ((mask.type & ATTRIBUTE.type) && !(type & ATTRIBUTE.type)) ||
         ((mask.type & SYSTEM.type)    && !(type & SYSTEM.type))    ||
         ((mask.type & PRN.type)       && !(type & PRN.type))       ||
         ((mask.type & FREQ_NO.type)   && !(type & FREQ_NO.type));
}

/***********************************************/

GnssType &GnssType::operator+=(const GnssType &t)
{
  try
  {
    if(((type & TYPE.type)      && (t.type & TYPE.type)      && (type & TYPE.type)       != (t.type & TYPE.type))       ||
       ((type & FREQUENCY.type) && (t.type & FREQUENCY.type) && ((type & FREQUENCY.type) != (t.type & FREQUENCY.type))) ||
       ((type & SYSTEM.type)    && (t.type & SYSTEM.type)    && ((type & SYSTEM.type)    != (t.type & SYSTEM.type)))    ||
       ((type & ATTRIBUTE.type) && (t.type & ATTRIBUTE.type) && ((type & ATTRIBUTE.type) != (t.type & ATTRIBUTE.type))) ||
       ((type & PRN.type)       && (t.type & PRN.type)       && ((type & PRN.type)       != (t.type & PRN.type)))       ||
       ((type & FREQ_NO.type)   && (t.type & FREQ_NO.type)   && ((type & FREQ_NO.type)   != (t.type & FREQ_NO.type))))
      throw(Exception("Incompatible addition: "+str()+" + "+t.str()));

    type |= t.type;
    return *this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssType &GnssType::operator&=(const GnssType &t)
{
  type &= t.type;
  return *this;
}

/***********************************************/

GnssType GnssType::operator~() const
{
  return GnssType(~type);
}

/***********************************************/

Bool GnssType::operator==(const GnssType &t) const
{
  if(type == t.type) return TRUE;
  if((type & SYSTEM.type)    && (t.type & SYSTEM.type)    && ((type & SYSTEM.type)    != (t.type & SYSTEM.type)))    return FALSE;
  if((type & FREQUENCY.type) && (t.type & FREQUENCY.type) && ((type & FREQUENCY.type) != (t.type & FREQUENCY.type))) return FALSE;
  if((type & TYPE.type)      && (t.type & TYPE.type)      && ((type & TYPE.type)      != (t.type & TYPE.type)))      return FALSE;
  if((type & ATTRIBUTE.type) && (t.type & ATTRIBUTE.type) && ((type & ATTRIBUTE.type) != (t.type & ATTRIBUTE.type))) return FALSE;
  if((type & PRN.type)       && (t.type & PRN.type)       && ((type & PRN.type)       != (t.type & PRN.type)))       return FALSE;
  if((type & FREQ_NO.type)   && (t.type & FREQ_NO.type)   && ((type & FREQ_NO.type)   != (t.type & FREQ_NO.type)))   return FALSE;
  return TRUE;
}

/***********************************************/

Bool GnssType::operator<(const GnssType &t) const
{
  if((type & SYSTEM.type)    != (t.type & SYSTEM.type))    return ((t.type & SYSTEM.type)    == 0) || (((type & SYSTEM.type)    != 0) && ((type & SYSTEM.type)    < (t.type & SYSTEM.type)));
  if((type & TYPE.type)      != (t.type & TYPE.type))      return ((t.type & TYPE.type)      == 0) || (((type & TYPE.type)      != 0) && ((type & TYPE.type)      < (t.type & TYPE.type)));
  if((type & FREQUENCY.type) != (t.type & FREQUENCY.type)) return ((t.type & FREQUENCY.type) == 0) || (((type & FREQUENCY.type) != 0) && ((type & FREQUENCY.type) < (t.type & FREQUENCY.type)));
  if((type & ATTRIBUTE.type) != (t.type & ATTRIBUTE.type)) return ((t.type & ATTRIBUTE.type) == 0) || (((type & ATTRIBUTE.type) != 0) && ((type & ATTRIBUTE.type) < (t.type & ATTRIBUTE.type)));
  if((type & PRN.type)       != (t.type & PRN.type))       return ((t.type & PRN.type)       == 0) || (((type & PRN.type)       != 0) && ((type & PRN.type)       < (t.type & PRN.type)));
  if((type & FREQ_NO.type)   != (t.type & FREQ_NO.type))   return ((t.type & FREQ_NO.type)   == 0) || (((type & FREQ_NO.type)   != 0) && ((type & FREQ_NO.type)   < (t.type & FREQ_NO.type)));
  return FALSE;
}

/***********************************************/
