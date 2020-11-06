/***********************************************/
/**
* @file fileNetCdf.h
*
* @brief Input/output of netCDF files.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-05-15
*
*/
/***********************************************/

#ifndef __GROOPS_FILENETCDF__
#define __GROOPS_FILENETCDF__

#include "base/import.h"
#include "inputOutput/fileName.h"

/** @brief GROOPS relevant subset of NetCDF data model */
namespace NetCdf
{
  enum DataType {BYTE, CHAR, SHORT, INT, FLOAT, DOUBLE};

  class Attribute;
  class Dimension;
  class Variable;

  /***** CLASS ***********************************/

  /** @brief netCDF Group representation. */
  class Group
  {
  protected:
    Int groupId; // netCDF Group ID

  public:
    /// Constructor from netCDF group ID
    Group(Int groupId=-1) : groupId(groupId) {}

    /// Adds a new dimension
    Dimension addDimension(const std::string &name, UInt length = MAX_UINT);

    /// Return all dimensions available in this group.
    std::vector<Dimension> dimensions() const;

    // -------------------

    /// Adds a new variable
    Variable addVariable(const std::string &name, DataType dtype, const std::vector<Dimension> &dims);

    /// Return all variables of this groups
    std::vector<Variable> variables() const;

    /// Returns the variable with name @a name.
    Variable variable(const std::string &name) const;

    // -------------------

    /// Add string attribute @a name
    Attribute addAttribute(const std::string &name, const std::string &val);

    /// Add value attribute @a name
    Attribute addAttribute(const std::string &name, DataType dtype, const Vector &val);

    /// Return all attributes of this group.
    std::vector<Attribute> attributes() const;
  };

  /***** CLASS ***********************************/

  /** @brief netCDF Dimension representation. */
  class Dimension
  {
    Int groupId; // containing group ID
    Int dimId;   // dimension ID

  public:
    /// Constructor from dimension id and group id
    Dimension(Int groupId=-1, Int dimId=-1) : groupId(groupId), dimId(dimId) {}

    /// Return the name of the dimension
    std::string name() const;

    /// Return the current length of the dimension
    UInt length() const;

    Bool operator==(const Dimension &x) const {return (groupId == x.groupId) && (dimId == x.dimId);}
    Bool operator!=(const Dimension &x) const {return !(*this == x);}

    friend class Group;
  };

  /***** CLASS ***********************************/

  /** @brief netCDF variable representation. */
  class Variable
  {
    Int groupId; // ID of containing group
    Int varId;   // variable ID

  public:
    /// Constructor from variable ID and group ID
    Variable(Int groupId=-1, Int varId=-1) : groupId(groupId), varId(varId) {}

    /// Return the name (string identifier) of the variable
    std::string name() const;

    /// Return all dimensions of the variable
    std::vector<Dimension> dimensions() const;

    // -------------------

    /// Add string attribute @a name
    Attribute addAttribute(const std::string &name, const std::string &text);

    /// Add value attribute @a name
    Attribute addAttribute(const std::string &name, DataType dtype, const Vector &val);

    /// Return all attributes of the variable
    std::vector<Attribute> attributes() const;

    /// Return the attribute with string identifier name
    Attribute attribute(const std::string &name) const;

    // -------------------

    /** @brief Write a slice of data values to the variable.
    * Double values will be cast to the data type of the variable.
    * @param start start index of the data point (must be of size @a dimensionCount)
    * @param count number of elements to be retrieved along each dimension
    * @param val data. */
    void setValues(const std::vector<UInt> &start, const std::vector<UInt> &count, const Vector &val) const;

    /** @brief Write a slice of data values to a 1D variable.
    * Double values will be cast to the data type of the variable.
    * @param val to a @a Vector which will hold the resulting data slice a contiguous chunk in memory. */
    void setValues(const Vector &val) const;

    /** @brief Retrieve a slice of data values from the variable.
    * All numeric data types will be cast to Double.
    * @param start start index of the data point (must be of size @a dimensionCount)
    * @param count number of elements to be retrieved along each dimension
    * @return Vector which will hold the resulting data slice a contiguous chunk in memory. */
    Vector values(const std::vector<std::size_t> &start, const std::vector<std::size_t> &count) const;

    /// Convenience function for 1D/2D variables which retrieves the whole data (useful for dimensions)
    Vector values() const;
  };

  /***** CLASS ***********************************/

  /** @brief netCDF Attribute representation. */
  class Attribute
  {
    Int groupId;        // ID of containing group
    Int varId;          // id of the associated variable
    std::string name_;  // name of the attribute

  public:
    /// Constructor from attribute name, group id and variable id
    Attribute(Int groupId, Int varId, const std::string &name) : groupId(groupId), varId(varId), name_(name) {}

    /// Return a readable reference to the name of the attribute
    const std::string &name() const {return name_;}

    /// Convenience function for string attributes
    std::string value() const;
  };

  /***** CLASS ***********************************/

  /** @brief File input entry point for netCDF data model. */
  class InFile : public Group
  {
  public:
    /// Constructor from file name. The specified netCDF file is immediately opened.
    InFile(const FileName &fileName);

    /// Destructor
   ~InFile();
  };

  /***** CLASS ***********************************/

  /** @brief File output entry point for netCDF data model. */
  class OutFile : public Group
  {
  public:
    /// Constructor from file name. The specified netCDF file is immediately opened.
    OutFile(const FileName &fileName);

    /// Destructor
   ~OutFile();
  };


  /***** FUNCTIONS ***********************************/

  /// convert a vector of generic dates/times to Time objects
  std::vector<Time> convertTimes(const_MatrixSliceRef values, const std::string &unitString);

  /// convert a vector of angles in arbitrary units to Angle objects
  std::vector<Angle> convertAngles(const_MatrixSliceRef values);
}

/***********************************************/

#endif
