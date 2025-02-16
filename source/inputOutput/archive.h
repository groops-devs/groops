/***********************************************/
/**
* @file archive.h
*
* @brief Interface for streaming of data in different formats.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#ifndef __GROOPS_ARCHIVE__
#define __GROOPS_ARCHIVE__

#include "base/importStd.h"

/** @defgroup archiveGroup Archive
* @brief Streaming of data in different formats.
* @ingroup inputOutputGroup */
/// @{

/***** TYPES ***********************************/

class Time;
class Angle;
class Vector;
class Matrix;
class MatrixSlice;
class const_MatrixSlice;
class SphericalHarmonics;
class Doodson;
class GnssType;
class Vector3d;
class Tensor3d;
class Rotary3d;
class Transform3d;
class Ellipsoid;

/***** CLASS ***********************************/

/** @brief Identifier for values in archives.
* Values need an identifier in xml format. */
template<typename T> class const_NameValue
{
  const char *_name;
  const T    *_value;
public:
  const_NameValue(const char *name, const T &value) : _name(name), _value(&value) {} //!< Constructor
  const char *name()  const {return _name;}   //!< name
  const T    &value() const {return *_value;} //!< value
};

/** @brief Adapter for const_NameValue. */
template<typename T>
inline const const_NameValue<T> nameValue(const char *name, const T &x) {return const_NameValue<T>(name,x);}

/***** CLASS ***********************************/

/** @brief Identifier for values in archives.
* Values need an identifier in xml format. */
template<typename T> class NameValue
{
  const char *_name;
  T          *_value;
public:
  NameValue(const char *name, T &value) : _name(name), _value(&value) {} //!< Constructor
  const char *name()  const {return _name;}   //!< name
  T          &value() const {return *_value;} //!< value
};

/** @brief Adapter for NameValue. */
template<typename T>
inline const NameValue<T> nameValue(const char *name, T &x) {return NameValue<T>(name,x);}


/***** CLASS ***********************************/

/** @brief Grouping values in archives.
* It is used in xml format. */
class BeginGroup
{
public:
  const char *name; //!< name
  BeginGroup(const char *_name) : name(_name) {} //!< Constructor
};

/** @brief Adapter for BeginGroup. */
inline const BeginGroup beginGroup(const char *name) {return BeginGroup(name);}

/***** CLASS ***********************************/

/** @brief Grouping values in archives.
* It is used in xml format. */
class EndGroup
{
public:
  const char *name; //!< name
  EndGroup(const char *_name) : name(_name) {} //!< Constructor
};

/** @brief Adapter for EndGroup. */
inline const EndGroup endGroup(const char *name) {return EndGroup(name);}

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/** @brief Abstract interface for output archives. */
class OutArchive
{
public:
  virtual ~OutArchive() {}

  enum ArchiveType {ASCII, XML, JSON, BINARY};
  virtual ArchiveType archiveType() const = 0;

  template<typename T> OutArchive &operator<<(const NameValue<T> &nv);
  template<typename T> OutArchive &operator<<(const const_NameValue<T> &nv);
  OutArchive &operator<<(const BeginGroup   &gr);
  OutArchive &operator<<(const EndGroup     &gr);

  virtual void save(const std::string &x) = 0;
  virtual void save(const Int      &x) = 0;
  virtual void save(const UInt     &x) = 0;
  virtual void save(const Double   &x) = 0;
  virtual void save(const Bool     &x) = 0;
  virtual void save(const Time     &x) = 0;
  virtual void save(const Angle    &x) = 0;
  virtual void save(const Vector   &x) = 0;
  virtual void save(const const_MatrixSlice  &x) = 0;
  virtual void save(const SphericalHarmonics &x) = 0;
  virtual void save(const Doodson  &x) = 0;
  virtual void save(const GnssType &x) = 0;

  virtual void startTag(const std::string &/*name*/) {}
  virtual void endTag  (const std::string &/*name*/) {}
  virtual void endLine() {}
  virtual void comment(const std::string &/*text*/) {}
};

/***** CLASS ***********************************/

/** @brief Abstract interface for input archives. */
class InArchive
{
protected:
  static UInt versionStr2version(const std::string &str);

public:
  virtual ~InArchive() {}

  template<typename T>
  InArchive &operator>>(const NameValue<T> &nv);
  InArchive &operator>>(const BeginGroup   &gr);
  InArchive &operator>>(const EndGroup     &gr);

  enum ArchiveType {ASCII, XML, JSON, BINARY};
  virtual ArchiveType archiveType() const = 0;
  virtual std::string type()        const = 0;
  virtual UInt        version()     const = 0;

  virtual void startTag(const std::string &/*name*/) {}
  virtual void endTag  (const std::string &/*name*/) {}

  virtual void load(std::string &x) = 0;
  virtual void load(Int      &x) = 0;
  virtual void load(UInt     &x) = 0;
  virtual void load(Double   &x) = 0;
  virtual void load(Bool     &x) = 0;
  virtual void load(Time     &x) = 0;
  virtual void load(Angle    &x) = 0;
  virtual void load(Vector   &x) = 0;
  virtual void load(Matrix   &x) = 0;
  virtual void load(SphericalHarmonics &x) = 0;
  virtual void load(Doodson  &x) = 0;
  virtual void load(GnssType &x) = 0;
};

/***** FUNCTIONS *******************************/

template<typename T> void save(OutArchive &ar, const std::vector<T> &x);
template<typename T> inline void save(OutArchive &ar, const T &x)        {x.save(ar);}
template<> inline void save(OutArchive &ar, const std::string &x)        {ar.save(x);}
template<> inline void save(OutArchive &ar, const Int     &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const UInt    &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Double  &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Bool    &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Time    &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Angle   &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Vector  &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const Matrix  &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const MatrixSlice  &x)       {ar.save(x);}
template<> inline void save(OutArchive &ar, const const_MatrixSlice  &x) {ar.save(x);}
template<> inline void save(OutArchive &ar, const SphericalHarmonics &x) {ar.save(x);}
template<> inline void save(OutArchive &ar, const Doodson &x)            {ar.save(x);}
template<> inline void save(OutArchive &ar, const GnssType &x)           {ar.save(x);}
template<> void save(OutArchive &ar, const Vector3d &x);
template<> void save(OutArchive &ar, const Tensor3d &x);
template<> void save(OutArchive &ar, const Rotary3d &x);
template<> void save(OutArchive &ar, const Transform3d &x);
template<> void save(OutArchive &ar, const Ellipsoid &x);

template<typename T> void load(InArchive  &ar, std::vector<T> &x);
template<typename T> inline void load(InArchive &ar, T &x)               {x.load(ar);}
template<> inline void load(InArchive  &ar, std::string &x)              {ar.load(x);}
template<> inline void load(InArchive  &ar, Int     &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, UInt    &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Double  &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Bool    &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Time    &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Angle   &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Vector  &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, Matrix  &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, SphericalHarmonics &x)       {ar.load(x);}
template<> inline void load(InArchive  &ar, Doodson &x)                  {ar.load(x);}
template<> inline void load(InArchive  &ar, GnssType &x)                 {ar.load(x);}
template<> void load(InArchive  &ar, Vector3d  &x);
template<> void load(InArchive  &ar, Tensor3d  &x);
template<> void load(InArchive  &ar, Rotary3d  &x);
template<> void load(InArchive  &ar, Transform3d &x);
template<> void load(InArchive  &ar, Ellipsoid &x);

/// @}

/***********************************************/
/***** INLINES ***********************************/
/***********************************************/

template<typename T> void save(OutArchive &ar, const std::vector<T> &x)
{
  UInt size = x.size();
  ar<<nameValue("count", size);
  ar.endLine();
  for(UInt i=0; i<size; i++)
  {
    ar<<nameValue("cell", const_cast<std::vector<T>&>(x).at(i));
    ar.endLine();
  }
}

/***********************************************/

template<typename T> void load(InArchive  &ar, std::vector<T> &x)
{
  UInt size;
  ar>>nameValue("count", size);
  x.resize(size);
  for(UInt i=0; i<size; i++)
    ar>>nameValue("cell", x.at(i));
}

/***********************************************/
/***********************************************/

inline OutArchive &OutArchive::operator<<(const BeginGroup &gr)
{
  startTag(gr.name);
  endLine();
  return *this;
}

/***********************************************/

inline InArchive &InArchive::operator>>(const BeginGroup &gr)
{
  startTag(gr.name);
  return *this;
}

/***********************************************/

inline OutArchive &OutArchive::operator<<(const EndGroup &gr)
{
  endTag(gr.name);
  endLine();
  return *this;
}

/***********************************************/

inline InArchive &InArchive::operator>>(const EndGroup &gr)
{
  endTag(gr.name);
  return *this;
}

/***********************************************/

template<typename T>
inline OutArchive &OutArchive::operator<<(const const_NameValue<T> &nv)
{
  startTag(nv.name());
  ::save(*this,nv.value());
  endTag(nv.name());
  return *this;
}

/***********************************************/

template<typename T>
inline OutArchive &OutArchive::operator<<(const NameValue<T> &nv)
{
  startTag(nv.name());
  ::save(*this,nv.value());
  endTag(nv.name());
  return *this;
}

/***********************************************/

template<typename T>
inline InArchive &InArchive::operator>>(const NameValue<T> &nv)
{
  startTag(nv.name());
  ::load(*this,nv.value());
  endTag(nv.name());
  return *this;
}

/***********************************************/
/***********************************************/

#endif
