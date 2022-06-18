/***********************************************/
/**
* @file fileInstrument.h
*
* @brief Satellites instrument data organized in arcs.
*
* @author Torsten Mayer-Guerr
* @date 2004-11-29
*
*/
/***********************************************/

#ifndef __GROOPS_FILEINSTRUMENT__
#define __GROOPS_FILEINSTRUMENT__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Instrument
static const char *docstringInstrument = R"(
This template file format can store different observations in a epoch wise manner. Each epoch consists of a time and
additional data, e.g orbits, accelerometer data, star camera quaternions (see \configClass{InstrumentType}{instrumentTypeType}).
The time series can be divided in several arcs (see \program{InstrumentSynchronize}).

Also a \file{matrix}{matrix} file is allowed as one single arc. The first column must contain times [MJD]. Without any extra column
the instrument type is INSTRUMENTTIME, with one additional column the type is MISCVALUE, and for more columns the type
MISCVALUES is used.

\begin{verbatim}
groops instrument version=20200123
# SATELLITETRACKING
         -9         60  # instrument type, number of arcs
# Time [MJD]               data0: range [m]          data1: range-rate [m/s]   data2: range-acc [m/s^2]
# =====================================================================================================
         12 # number of epochs of 1. arc
 54588.000000000000000000 -5.074649470097549492e+05  5.755440207134928654e-01  1.877605261528093308e-03
 54588.000057870370255841 -5.074620458130163024e+05  5.849357691551860805e-01  1.878948916234051596e-03
 54588.000115740740966430 -5.074590976427756250e+05  5.943331739937073310e-01  1.879937220634776869e-03
 54588.000173611111222272 -5.074561024756557308e+05  6.037340169611068452e-01  1.880370529387525701e-03
 54588.000231481481478113 -5.074530602992626373e+05  6.131368121270999172e-01  1.880680632122925426e-03
 54588.000289351851733954 -5.074499711071007187e+05  6.225398878861636565e-01  1.880495369480403561e-03
 54588.000347222222444543 -5.074468349029610981e+05  6.319414138081351773e-01  1.880073731783055927e-03
 54588.000405092592700385 -5.074436516971451929e+05  6.413404243585696385e-01  1.879464843086203459e-03
 54588.000462962962956226 -5.074404215058300761e+05  6.507353310092597320e-01  1.878578987216372124e-03
 54588.000520833333212067 -5.074371443491023383e+05  6.601267978060636477e-01  1.877878184949659246e-03
 54588.000578703703922656 -5.074338202460713219e+05  6.695136489207137442e-01  1.876962042758626532e-03
 54588.000636574074178498 -5.074304492190054734e+05  6.788964444122400632e-01  1.876091925462087043e-03
         12 # number of epochs of 2. arc
 54588.000694444444434339 -5.074270312892858055e+05  6.882748400534359767e-01  1.875376456928801432e-03
 54588.000752314814690180 -5.074235664742725785e+05  6.976508178537534910e-01  1.874929898412159559e-03
 54588.000810185185400769 -5.074200547868391732e+05  7.070236200716006891e-01  1.874312324351668077e-03
 54588.000868055555656611 -5.074164962409950094e+05  7.163943828291452487e-01  1.873924188388115340e-03
 54588.000925925925912452 -5.074128908454515622e+05  7.257639682023964145e-01  1.874025826380292404e-03
 54588.000983796296168293 -5.074092386012640782e+05  7.351333608427884636e-01  1.873680487441316657e-03
 54588.001041666666878882 -5.074055395130896359e+05  7.445020815182646912e-01  1.873849502509668122e-03
 54588.001099537037134724 -5.074017935789784533e+05  7.538716732272922050e-01  1.873971633320137753e-03
 54588.001157407407390565 -5.073980007962241652e+05  7.632414098560330595e-01  1.873984767500571974e-03
 54588.001215277777646406 -5.073941611626467784e+05  7.726123093411200182e-01  1.874295246964456478e-03
 54588.001273148148356995 -5.073902746728868224e+05  7.819835205798950639e-01  1.874226146744964808e-03
 54588.001331018518612836 -5.073863413272026228e+05  7.913547196412918927e-01  1.874173804634685515e-03
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/gnssType.h"
#include "inputOutput/fileArchive.h"

/**
* @defgroup fileInstrumentGroup FileInstrument
* @brief Satellite Instrument data.
* @ingroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_INSTRUMENT_TYPE = "instrument";

/***** CLASS ***********************************/

/** @brief Epoch with satellite instrument data.
* An Epoch contain instrument data given at a specific @a time.
* This is an abstract class. */
class Epoch
{
public:
  /// Instrument data type (type > 0: MISCVALUES where type is the number of data columns).
  enum Type : Int  {EMPTY             =  0,
                    MISCVALUESOLD     = -9999,
                    MISCVALUES        = -1, // MISCVALUES with unknown or zero size
                    INSTRUMENTTIME    = -2,
                    MISCVALUE         = -3,
                    VECTOR3D          = -4,
                    COVARIANCE3D      = -5,
                    ORBIT             = -6,
                    STARCAMERA        = -7,
                    ACCELEROMETER     = -8,
                    SATELLITETRACKING = -9,
                    GRADIOMETER       = -10,
                    GNSSRECEIVER      = -11,
                    OBSERVATIONSIGMA  = -12,
                    MASS              = -16,
                    THRUSTER          = -17,
                    MAGNETOMETER      = -18,
                    ACCHOUSEKEEPING   = -19,
                    CLOCK             = -20,
                    STARCAMERA1A      = -21,
                    ACCELEROMETER1A   = -22};

  Time time; //!< Time of Epoch.

  Epoch()                         = default; //!< Constructor
  Epoch(const Epoch &)            = default; //!< Copy constructor
  Epoch(Epoch &&)                 = default; //!< Move constructor
  virtual ~Epoch()                = default; //!< Destructor.
  Epoch &operator=(const Epoch &) = default; //!< Assignement
  Epoch &operator=(Epoch &&)      = default; //!< Move assignment

  /** @brief Create an epoch of given type (with new). */
  static Epoch *create(Type type);

  /** @brief Data type (e.g. ORBIT, ACCELEROMETER). */
  virtual Type getType() const = 0;

  /** @brief Name of data type (e.g. ORBIT, ACCELEROMETER). */
  std::string getTypeName() const {return getTypeName(getType());}

  /** @brief Name of data type (e.g. ORBIT, ACCELEROMETER). */
  static std::string getTypeName(Type type);

  /** @brief Description of  data type (e.g. ORBIT, ACCELEROMETER). */
  static std::string fileFormatString(Type type);

  /** @brief Number of data columns
  * If undefined NULLINDEX is returned
  * and if additionally @p mustDefined==TRUE, an Exception is thrown. */
  static UInt dataCount(Type type, Bool mustDefined=FALSE);

  /** @brief The data as as list of Doubles.
  * Without time. */
  virtual Vector data() const = 0; // data without time

  /** @brief set data from a list of Doubles.
  * Without time. */
  virtual void setData(const Vector &x) = 0; // data without time

  /** @brief A copy is created with new.
  * The user has to delete the new Epoch. */
  virtual Epoch *clone() const = 0;

  virtual void save(OutArchive &oa) const = 0;
  virtual void load(InArchive  &ia) = 0;
};

/***** CLASS ***********************************/

/** @brief Arc with satellite instrument data.
* An Arc consists of a list of Epochs of an arbitrary instrument. */
class Arc
{
protected:
  std::vector<std::unique_ptr<Epoch>> epoch;

  virtual void checkType(Epoch::Type type)
  {
    if((type != Epoch::EMPTY) && (getType() != Epoch::EMPTY) && (type != getType()))
      throw(Exception("In Arc<"+Epoch::getTypeName(getType())+">: assignment of wrong data type: "+Epoch::getTypeName(type)));
  }

  /** @brief Random access iterator with implicit double dereferencing to directly access epochs of type @a T.
   * Constness of iterator is defined via template parameter @a isConstIterator. */
  template<typename T, bool isConstIterator>
  class Iterator : public std::conditional<isConstIterator, std::vector<std::unique_ptr<Epoch>>::const_iterator, std::vector<std::unique_ptr<Epoch>>::iterator>::type
  {
    using BaseIterator  = typename std::conditional<isConstIterator, std::vector<std::unique_ptr<Epoch>>::const_iterator, std::vector<std::unique_ptr<Epoch>>::iterator>::type;
    using PointerType   = typename std::conditional<isConstIterator, const T*, T*>::type;
    using ReferenceType = typename std::conditional<isConstIterator, const T&, T&>::type;

  public:
    using iterator_category = typename BaseIterator::iterator_category;
    using difference_type   = typename BaseIterator::difference_type;
    using value_type        = T;
    using pointer           = PointerType;
    using reference         = ReferenceType;

    Iterator() : BaseIterator() {}                           /// Default constructor
    Iterator(const BaseIterator &it) : BaseIterator(it) {}   /// Base constructor

    PointerType operator->() const { return dynamic_cast<PointerType>(BaseIterator::operator*().get()); }                           // double dereferencing
    ReferenceType operator*() const { return *this->operator->(); }                                                                 // double dereferencing
    ReferenceType operator[](difference_type rhs) const { return *dynamic_cast<PointerType>(BaseIterator::operator[](rhs).get()); } // double dereferencing
    Iterator& operator+=(difference_type rhs)           { BaseIterator::operator+=(rhs); return *this; }
    Iterator& operator-=(difference_type rhs)           { BaseIterator::operator-=(rhs); return *this; }
    Iterator& operator++()                              { BaseIterator::operator++(); return *this; }
    Iterator& operator--()                              { BaseIterator::operator--(); return *this; }
    Iterator operator++(int rhs)                        { return BaseIterator::operator++(rhs); }
    Iterator operator--(int rhs)                        { return BaseIterator::operator--(rhs); }
    Iterator operator+(difference_type rhs) const       { return BaseIterator::operator+(rhs); }
    Iterator operator-(difference_type rhs) const       { return BaseIterator::operator-(rhs); }
    friend Iterator operator+(difference_type lhs, const Iterator &rhs) { return lhs+BaseIterator(rhs); }
  };

public:
  using iterator               = Iterator<Epoch, false>;
  using const_iterator         = Iterator<Epoch, true>;
  using reverse_iterator       = typename std::reverse_iterator<iterator>;
  using const_reverse_iterator = typename std::reverse_iterator<const_iterator>;

  iterator begin()                       { return epoch.begin(); }
  iterator end()                         { return epoch.end(); }
  const_iterator begin()           const { return epoch.cbegin(); }
  const_iterator end()             const { return epoch.cend(); }
  const_iterator cbegin()          const { return epoch.cbegin(); }
  const_iterator cend()            const { return epoch.cend(); }
  reverse_iterator rbegin()              { return reverse_iterator(this->end()); }
  reverse_iterator rend()                { return reverse_iterator(this->begin()); }
  const_reverse_iterator rbegin()  const { return const_reverse_iterator(this->cend()); }
  const_reverse_iterator rend()    const { return const_reverse_iterator(this->cbegin()); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(this->cend()); }
  const_reverse_iterator crend()   const { return const_reverse_iterator(this->cbegin()); }

public:
  Arc() = default;          //!< Default Constructor
  Arc(const Arc &arc);      //!< Copy constructor
  Arc(Arc &&) = default;    //!< Move constructor
  virtual ~Arc() = default; //!< Destructor

  /// Constructor from matrix (first column is MJD)
  explicit Arc(const_MatrixSliceRef A, Epoch::Type type=Epoch::EMPTY);

  /** @brief Constructor from matrix (first column is MJD)
  First columns of @a A is ignored, instead @a times is used. */
  explicit Arc(const std::vector<Time> &times, const_MatrixSliceRef A, Epoch::Type type=Epoch::EMPTY);

  Arc &operator=(const Arc &x);     //!< Assignement
  Arc &operator=(Arc &&) = default; //!< Move assignment

  /** @brief Epoch at index @a i. */
  const Epoch &at(UInt i) const {return *epoch.at(i);}

  /** @brief Writeable reference to epoch at index @a i. */
  Epoch &at(UInt i) {return *epoch.at(i);}

  /** @brief First epoch in arc */
  const Epoch &front() const {return at(0);}
        Epoch &front()       {return at(0);}

  /** @brief Last epoch in arc */
  const Epoch &back() const {return at(size()-1);}
        Epoch &back()       {return at(size()-1);}

  /** @brief Number of epochs. */
  UInt size() const {return epoch.size();}

  /** @brief Data type (e.g. ORBIT, ACCELEROMETER). */
  Epoch::Type getType() const {return epoch.size() ? front().getType() : Epoch::EMPTY;}

  /** @brief Name of data type (e.g. ORBIT, ACCELEROMETER). */
  std::string getTypeName() const {return Epoch::getTypeName(getType());}

  /** @brief Sort epochs in ascending temporal order. */
  void sort();

  /** @brief Remove consecutive epochs with equal time stamps. */
  void removeDuplicateEpochs(Bool keepFirst = TRUE, Double margin=1e-5);

  /** @brief Creates a new arc with @a len epochs starting at index @a start. */
  Arc subArc(UInt start, UInt len) const;

  /** @brief Removes @a len epochs starting at index @a index. */
  void remove(UInt index, UInt len=1);

  /** @brief Append another arc. */
  void append(const Arc &arc);

  /** @brief Synchronize all epochs with a given time list.
  * All epochs are removed which are not in the time list.
  * @param time list of valid points in time (all other epochs will be removed)
  * @param margin time interval around the given @a time in which epochs are valid [seconds] */
  void synchronize(const std::vector<Time> &time, Double margin=1e-5);

  /** @brief Divide arc at time gaps and append new sub arcs to @a arcList.
  * @param minGap minimal time span for a gap.
  * @param minArcLen arcs shorter than minArcLen are dropped
  * @param arcList divided arcs are appended to this list. */
  void divide(const Time &minGap, UInt minArcLen, std::vector<Arc> &arcList) const;

  /** @brief Insert a copy of @a epoch at @a index. */
  void insert(UInt index, const Epoch &epoch);

  /** @brief Append a copy of @a epoch to the arc. */
  void push_back(const Epoch &epoch);

  /** @brief Time series of epochs. */
  std::vector<Time> times() const;

  /** @brief Time series of data as matrix (first column is MJD). */
  Matrix matrix() const;

  void load(InArchive  &ia);
  void save(OutArchive &oa) const;

  template<typename T>
  static Epoch::Type getType(const T &arcList)
  {
    // determine type
    Epoch::Type type = Epoch::EMPTY;
    for(const Arc &arc : arcList)
      if(arc.size())
      {
        if((type != Epoch::EMPTY) && (type != arc.getType()))
          throw(Exception("arcList contain different instruments types "+Epoch::getTypeName(type)+", "+arc.getTypeName()));
        type = arc.getType();
      }
    return type;
  }

  /** @brief Test of synchronicity of multiple arcs.
  * If the arcs are not synchronous an expection is thrown.
  * Empty arcs are ignored. */
  static void checkSynchronized(const std::vector<std::reference_wrapper<const Arc>> &arc);

  /** @brief Log information about arc size, number of epochs and so on. */
  static void printStatistics(const Arc &arc);

  /** @brief Log information about arc size, number of epochs and so on. */
  static void printStatistics(const std::vector<Arc> &arcList);
};

/***** CLASS ***********************************/

/** @brief Arc with specific satellite instrument data.
* An Arc consists of a list of Epoch of specific instrument data. */
template<class EpochType>
class ArcTemplate : public Arc
{
  void checkType(Epoch::Type type) override
  {
    Arc::checkType(type);
    if((type != Epoch::EMPTY) && (type != EpochType::type) && !((EpochType::type == Epoch::MISCVALUES) && (static_cast<Int>(type) > 0)))
      throw(Exception("In Arc<"+Epoch::getTypeName(EpochType::type)+">: assignment of wrong data type: "+Epoch::getTypeName(type)));
  }

public:
  using iterator               = Iterator<EpochType, false>;
  using const_iterator         = Iterator<EpochType, true>;
  using reverse_iterator       = typename std::reverse_iterator<iterator>;
  using const_reverse_iterator = typename std::reverse_iterator<const_iterator>;

  iterator begin()                       { return epoch.begin(); }
  iterator end()                         { return epoch.end(); }
  const_iterator begin()           const { return epoch.cbegin(); }
  const_iterator end()             const { return epoch.cend(); }
  const_iterator cbegin()          const { return epoch.cbegin(); }
  const_iterator cend()            const { return epoch.cend(); }
  reverse_iterator rbegin()              { return reverse_iterator(this->end()); }
  reverse_iterator rend()                { return reverse_iterator(this->begin()); }
  const_reverse_iterator rbegin()  const { return const_reverse_iterator(this->cend()); }
  const_reverse_iterator rend()    const { return const_reverse_iterator(this->cbegin()); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(this->cend()); }
  const_reverse_iterator crend()   const { return const_reverse_iterator(this->cbegin()); }

public:
  ArcTemplate() = default;                                             //!< Default Constructor
  ArcTemplate(const ArcTemplate<EpochType> &arc) = default;            //!< Copy constructor
  ArcTemplate(const Arc &arc) : Arc(arc) {checkType(getType());}       //!< Copy constructor
  ArcTemplate(ArcTemplate<EpochType> &&arc) = default;                 //!< Move constructor
  ArcTemplate(Arc &&arc) : Arc(std::move(arc)) {checkType(getType());} //!< Move constructor

  ArcTemplate &operator=(const ArcTemplate<EpochType> &arc) {Arc::operator=(arc); return *this;}                                  //!< Assignement
  ArcTemplate &operator=(const Arc &arc)                    {Arc::operator=(arc); checkType(getType()); return *this;}            //!< Assignement
  ArcTemplate &operator=(ArcTemplate<EpochType> &&arc)      {Arc::operator=(std::move(arc)); return *this;}                       //!< Move assignment
  ArcTemplate &operator=(Arc &&arc)                         {Arc::operator=(std::move(arc)); checkType(getType()); return *this;} //!< Move assignment

  /** @brief Instrument specific epoch at index @a i. */
  const EpochType &at(UInt i) const  {return *dynamic_cast<EpochType*>(epoch.at(i).get());}

  /** @brief Instrument specific writable epoch at index @a i. */
  EpochType &at(UInt i) {return *dynamic_cast<EpochType*>(epoch.at(i).get());}

  /** @brief First epoch in arc */
  const EpochType &front() const {return at(0);}
        EpochType &front()       {return at(0);}

  /** @brief Last epoch in arc */
  const EpochType &back() const {return at(size()-1);}
        EpochType &back()       {return at(size()-1);}
};

/***** TYPES ***********************************/

class InstrumentFile;
typedef std::shared_ptr<InstrumentFile> InstrumentFilePtr;

/***** CLASS ***********************************/

/** @brief File with arbitrary satellite instrument data.
* A File consists of a list of arcs (Arc).
* Each arc consists of a list of Epoch of a specific instrument
* (Only one instrument type per file is allowed).
* The file is not read at once, but only read arc by arc.
*
* @code
* InstrumentFile file(FileName("grace_satelliteTracking.dat"));
* UInt arcCount = file.arcCount();
* for(UInt arcNo=0; arcNo<arcCount; arcNo++)
* {
*   SatelliteTrackingrArc sstArc = file.readArc(arcNo); // read arc
*   for(UInt i=0; i<sstArc.size(); i++)
*     logInfo<<arc.at(i).time<<": "<<arc.at(i).range<<Log::endl; // access to Epoch with .at(i)
* }
* @endcode
*/
class InstrumentFile
{
  InFileArchive file;
  FileName      fileName;
  Epoch::Type   type;
  UInt          arcCount_;
  UInt          index;
  Matrix        A; // if a matrix file is open

public:
  InstrumentFile() : type(Epoch::EMPTY), arcCount_(0) {}       //!< Default constructor.
  explicit InstrumentFile(const FileName &name) {open(name);}  //!< Constructor.
  InstrumentFile(const InstrumentFile &) = delete;             //!< Disallow copy constructor
  InstrumentFile &operator=(const InstrumentFile &x) = delete; //!< Disallow copying
 ~InstrumentFile() {close();}                                  //!< Destructor.

  /** @brief Open a new file.
  * Old open file is closed before. */
  void  open(const FileName &name);

  /** @brief Close the file. */
  void  close();

  /** @brief Number of Arc in file. */
  UInt arcCount() const {return arcCount_;}

  /** @brief Instrument type. */
  Epoch::Type getType()  const {return type;}

  /** @brief Name of data type (e.g. ORBIT, ACCELEROMETER). */
  std::string getTypeName() const {return Epoch::getTypeName(getType());}

  /** @brief Number of data columns.
  * If undefined NULLINDEX is returned
  * and if additionally @p mustDefined==TRUE, an Exception is thrown. */
  UInt dataCount(Bool mustDefined=FALSE) const {return Epoch::dataCount(getType(), mustDefined);}

  /** @brief Read a single Arc.
  * The operation is faster, if the arcs in read in increasing order.
  * If the file is not open, a empty Arc is returned. */
  Arc readArc(UInt arcNo);

  /** @brief Test number of arcs of multiple files.
  * Test whether files are divided into the same number of arcs otherwise an expection is thrown.
  * Files which are not open are ignored. */
  static void checkArcCount(const std::vector<std::reference_wrapper<const InstrumentFile>> &fileList);

  /** @brief Read all arcs and concatenate to one arc. */
  static Arc read(const FileName &name);

  /** @brief Write an Arc to file. */
  static void write(const FileName &name, const Arc &arc) {write(name, std::vector<Arc>(1, arc));}

  /** @brief Write an Arc to file. */
  template<typename T>
  static void write(const FileName &name, const ArcTemplate<T> &arc) {write(name, std::vector<Arc>(1, arc));}

  /** @brief Write a list of Arc to file. */
  template<typename Container>
  static void write(const FileName &name, const Container &arcList)
  {
    try
    {
      OutFileArchive file(name, FILE_INSTRUMENT_TYPE);
      Epoch::Type type = Arc::getType(arcList);
      file.comment(Epoch::getTypeName(type));
      file<<nameValue("satelliteType", static_cast<Int>(type));
      file<<nameValue("arcCount",      arcList.size());
      const std::string comment = Epoch::fileFormatString(type);
      file.comment(comment);
      file.comment(std::string(comment.size(), '='));
      for(const Arc &arc : arcList)
      {
        file<<beginGroup("arc");
        file<<nameValue("pointCount", arc.size());
        for(UInt i=0; i<arc.size(); i++)
        {
          file<<beginGroup("epoch");
          arc.at(i).save(file.outArchive());
          file<<endGroup("epoch");
        }
        file<<endGroup("arc");
      }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }

  /** @brief Factory for Instrument file. */
  static InstrumentFilePtr newFile(const std::string &name="") {return std::make_shared<InstrumentFile>(name);}
};

/***********************************************/

/** @brief Write into an Instrument file. */
inline void writeFileInstrument(const FileName &name, const Arc &arc) {InstrumentFile::write(name, arc);}
template<typename T> inline void writeFileInstrument(const FileName &name, ArcTemplate<T> &arc) {InstrumentFile::write(name, arc);}
template<typename Container> inline void writeFileInstrument(const FileName &name, const Container &arcList) {InstrumentFile::write(name, arcList);}

/***********************************************/
/***** Derived data types **********************/
/***********************************************/
// implementation of a new data
// 1. define new enum in Epoch (see top)
// 2. expand Epoch::create(Type type) in fileSatellite.cpp
// 3. expand Epoch::getTypeName(Type type) in fileSatellite.cpp
// 4. expand Epoch::dataCount(Type type) in fileSatellite.cpp
// 5. implement load, save, data, and setData in fileSatellite.cpp
/***********************************************/

/** @brief Epoch with a vector of values. */
class MiscValuesEpoch : public Epoch
{
public:
  Vector values;

  MiscValuesEpoch(UInt size) : values(size) {}

  static constexpr Type type = MISCVALUES;
  virtual Type   getType() const {return values.size() ? static_cast<Type>(values.size()) : MISCVALUES;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new MiscValuesEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<MiscValuesEpoch> MiscValuesArc;

/***********************************************/

/** @brief Epoch with time data only. */
class InstrumentTimeEpoch : public Epoch
{
public:
  static constexpr Type type = INSTRUMENTTIME;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new InstrumentTimeEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<InstrumentTimeEpoch> InstrumentTimeArc;

/***********************************************/

/** @brief Epoch with a single value. */
class MiscValueEpoch : public Epoch
{
public:
  Double value;

  MiscValueEpoch() : value(0) {}

  static constexpr Type type = MISCVALUE;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new MiscValueEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<MiscValueEpoch> MiscValueArc;

/***********************************************/

/** @brief Epoch with a vector of values. */
class Vector3dEpoch : public Epoch
{
public:
  Vector3d vector3d;

  static constexpr Type type = VECTOR3D;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new Vector3dEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<Vector3dEpoch> Vector3dArc;

/***********************************************/

/** @brief Epoch with 3x3 covariance matrix of Vector3d data. */
class Covariance3dEpoch : public Epoch
{
public:
  Tensor3d covariance;

  static constexpr Type type = COVARIANCE3D;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new Covariance3dEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<Covariance3dEpoch> Covariance3dArc;

/***********************************************/

/** @brief Epoch with orbit data. */
class OrbitEpoch : public Epoch
{
public:
  Vector3d position;
  Vector3d velocity;
  Vector3d acceleration;

  static constexpr Type type = ORBIT;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new OrbitEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<OrbitEpoch> OrbitArc;

/***********************************************/

/** @brief Epoch with star camera data (rotation). */
class StarCameraEpoch : public Epoch
{
public:
  Rotary3d rotary;

  static constexpr Type type = STARCAMERA;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new StarCameraEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<StarCameraEpoch> StarCameraArc;

/***********************************************/

/** @brief Epoch with accelerometer data (accelerations). */
class AccelerometerEpoch : public Epoch
{
public:
  Vector3d acceleration;

  static constexpr Type type = ACCELEROMETER;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new AccelerometerEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<AccelerometerEpoch> AccelerometerArc;

/***********************************************/

/** @brief Epoch with (low-low) satellite tracking data. */
class SatelliteTrackingEpoch : public Epoch
{
public:
  Double range;
  Double rangeRate;
  Double rangeAcceleration;

  SatelliteTrackingEpoch() : range(0), rangeRate(0), rangeAcceleration(0) {}

  static constexpr Type type = SATELLITETRACKING;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new SatelliteTrackingEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<SatelliteTrackingEpoch> SatelliteTrackingArc;

/***********************************************/

/** @brief Epoch with gradiometer data (gravity gradients). */
class GradiometerEpoch : public Epoch
{
public:
  Tensor3d gravityGradient;

  static constexpr Type type = GRADIOMETER;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new GradiometerEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<GradiometerEpoch> GradiometerArc;

/***********************************************/

/** @brief Epoch with GNSS receiver observations. */
class GnssReceiverEpoch : public Epoch
{
public:
  std::vector<GnssType> obsType;
  std::vector<GnssType> satellite;
  std::vector<Double>   observation;
  Double                clockError;

  GnssReceiverEpoch() : clockError(0.) {}

  static constexpr Type type = GNSSRECEIVER;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new GnssReceiverEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<GnssReceiverEpoch> GnssReceiverArc;

/***********************************************/
/***********************************************/

/** @brief Epoch with observations sigma. */
class ObservationSigmaEpoch : public Epoch
{
public:
  Double sigma;

  ObservationSigmaEpoch() : sigma(0) {}

  static constexpr Type type = OBSERVATIONSIGMA;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new ObservationSigmaEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<ObservationSigmaEpoch> ObservationSigmaArc;

/***********************************************/

/** @brief Epoch with mass data. */
class MassEpoch : public Epoch
{
public:
  Double massThr;
  Double massTank;

  static constexpr Type type = MASS;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new MassEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<MassEpoch> MassArc;

/***********************************************/

/** @brief Epoch with thruster data. */
class ThrusterEpoch : public Epoch
{
public:
  UInt onTime1, onTime2, onTime3,  onTime4,  onTime5,  onTime6,  onTime7,
         onTime8, onTime9, onTime10, onTime11, onTime12, onTime13, onTime14;

  static constexpr Type type = THRUSTER;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new ThrusterEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<ThrusterEpoch> ThrusterArc;

/***********************************************/

/** @brief Epoch with magnetometer data. */
class MagnetometerEpoch : public Epoch
{
public:
  Vector3d magneticField;
  Vector3d torquerA;
  Vector3d torquerB;
  Vector3d magneticFieldCalibration;
  Double   torquerCalibration;

  static constexpr Type type = MAGNETOMETER;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new MagnetometerEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<MagnetometerEpoch> MagnetometerArc;

/***********************************************/

/** @brief Epoch with IPU housekeeping data. */
class AccHousekeepingEpoch : public Epoch
{
public:
  Double   biasVoltage; // proof mass bias voltage (averaged) [V]
  Double   vd;          // amplitude of the AC voltages that operates the position sensors [Vrms]
  Vector3d xOut;        // displacement of capacitive sensor X1, X2, X3 [m]
  Vector3d yOut;        // displacement of capacitive sensor Y2, Y2, Z1 [m]
  Double   tempSU;      // temperature of SU electronics          [°C]
  Double   tempICU;     // temperature of ICU power supply board  [°C]
  Double   tempCore;    // temperature of internal core           [°C]
  Double   tempICUConv; // temperature of ICU A/D converter board [°C]
  UInt     blkNrICU;    // ICU block number

  static constexpr Type type = ACCHOUSEKEEPING;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new AccHousekeepingEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<AccHousekeepingEpoch> AccHousekeepingArc;

/***********************************************/

/** @brief Epoch with clock data. */
class ClockEpoch : public Epoch
{
public:
  Int    rcvTime;
  Double epsTime;
  Double epsError, epsDrift, driftError;
  UInt   qualityFlag;

  static constexpr Type type = CLOCK;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new ClockEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<ClockEpoch> ClockArc;

/***********************************************/

/** @brief Epoch with Level-1A star camera data (quaternions). */
class StarCamera1AEpoch : public Epoch
{
public:
  Int    rcvTime;
  Double epsTime;
  UInt   scaDesign; // 1 = primary SCA, 2 = secondary SCA
  Double q0, q1, q2, q3;
  UInt   nLocks;
  UInt   nStars;

  static constexpr Type type = STARCAMERA1A;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new StarCamera1AEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<StarCamera1AEpoch> StarCamera1AArc;

/***********************************************/

/** @brief Epoch with Level-1A accelerometer data (accelerations). */
class Accelerometer1AEpoch : public Epoch
{
public:
  Int      rcvTimeInt;
  Double   rcvTimeFrac;
  Vector3d acceleration;

  static constexpr Type type = ACCELEROMETER1A;
  virtual Type   getType() const {return type;}
  virtual Vector data()    const; // data without time
  virtual void   setData(const Vector &x);
  virtual Epoch *clone()   const {return new Accelerometer1AEpoch(*this);}
  virtual void   save(OutArchive &oa) const;
  virtual void   load(InArchive  &ia);
};

typedef ArcTemplate<Accelerometer1AEpoch> Accelerometer1AArc;

/***********************************************/

/// @}

#endif // __GROOPS_FILESATELLITE__
