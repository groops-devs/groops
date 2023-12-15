/***********************************************/
/**
* @file fileNetCdf.cpp
*
* @brief Input/output of netCDF files.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-05-15
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "fileNetCdf.h"

#ifdef GROOPS_DISABLE_NETCDF
#else
#include <netcdf.h>

/***********************************************/

template<typename T> static void convert(const Vector &values, std::vector<T> &tmp)
{
  tmp.resize(values.rows());
  for(UInt i=0; i<values.rows(); i++)
    tmp.at(i) = static_cast<T>(values(i));
}

/***********************************************/

template<typename T> static Vector convert(const std::vector<T> &tmp)
{
  Vector values(tmp.size());
  for(UInt i=0; i<values.rows(); i++)
    values(i) = static_cast<Double>(tmp.at(i));
  return values;
}

/***********************************************/
/***********************************************/

NetCdf::Dimension NetCdf::Group::addDimension(const std::string &name, UInt length)
{
  if(length == MAX_UINT)
    length = NC_UNLIMITED;

  Int dimId = 0;
  nc_def_dim(groupId, name.c_str(), length, &dimId);

  return Dimension(groupId, dimId);
}

/***********************************************/

std::vector<NetCdf::Dimension> NetCdf::Group::dimensions() const
{
  Int countDimensions = 0;
  nc_inq_ndims(groupId, &countDimensions);
  std::vector<Dimension> dims;
  for(Int dimId=0; dimId<countDimensions; dimId++)
    dims.push_back(NetCdf::Dimension(groupId, dimId));
  return dims;
}

/***********************************************/

NetCdf::Variable NetCdf::Group::addVariable(const std::string &name, DataType dtype, const std::vector<NetCdf::Dimension> &dims)
{
  std::vector<Int> dimIds(dims.size());
  for(UInt i=0; i<dims.size(); i++)
    dimIds.at(i) = dims.at(i).dimId;

  nc_type ncType = NC_DOUBLE;
  switch(dtype)
  {
    case BYTE:   ncType = NC_BYTE;   break;
    case CHAR:   ncType = NC_CHAR;   break;
    case SHORT:  ncType = NC_SHORT;  break;
    case INT:    ncType = NC_INT;    break;
    case FLOAT:  ncType = NC_FLOAT;  break;
    case DOUBLE: ncType = NC_DOUBLE; break;
  }

  Int varId = 0;
  nc_def_var(groupId, name.c_str(), ncType, dimIds.size(), dimIds.data(), &varId);

  return Variable(groupId, varId);
}

/***********************************************/

std::vector<NetCdf::Variable> NetCdf::Group::variables() const
{
  std::vector<Variable> vars;

  Int countVariables = 0;
  nc_inq_nvars(groupId, &countVariables);
  std::vector<Int> varIds(countVariables);
  nc_inq_varids(groupId, &countVariables, varIds.data());

  for(Int varId : varIds)
    vars.push_back(Variable(groupId, varId));

  return vars;
}

/***********************************************/

Bool NetCdf::Group::hasVariable(const std::string &name) const
{
  Int varId;
  return (nc_inq_varid(groupId, name.c_str(), &varId) == NC_NOERR);
}

/***********************************************/

NetCdf::Variable NetCdf::Group::variable(const std::string &name) const
{
  try
  {
    Int varId;
    if(nc_inq_varid(groupId, name.c_str(), &varId) != NC_NOERR)
      throw(Exception("Variable <"+name+"> not found in group"));
    return Variable(groupId, varId);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

NetCdf::Attribute NetCdf::Group::addAttribute(const std::string &name, const std::string &value)
{
  nc_put_att_text(groupId, NC_GLOBAL, name.c_str(), value.size(), value.c_str());
  return Attribute(groupId, NC_GLOBAL, name);
}

/***********************************************/

NetCdf::Attribute NetCdf::Group::addAttribute(const std::string &name, DataType dtype, const Vector &val)
{
  switch(dtype)
  {
    case BYTE:   {std::vector<Byte>  tmp; convert(val, tmp); nc_put_att      (groupId, NC_GLOBAL, name.c_str(), NC_BYTE,  tmp.size(), reinterpret_cast<void*>(tmp.data())); break;}
    case CHAR:   {std::vector<char>  tmp; convert(val, tmp); nc_put_att_text (groupId, NC_GLOBAL, name.c_str(), tmp.size(), tmp.data());  break;}
    case SHORT:  {std::vector<short> tmp; convert(val, tmp); nc_put_att_short(groupId, NC_GLOBAL, name.c_str(), NC_SHORT, tmp.size(), tmp.data());  break;}
    case INT:    {std::vector<int>   tmp; convert(val, tmp); nc_put_att_int  (groupId, NC_GLOBAL, name.c_str(), NC_INT,   tmp.size(), tmp.data());  break;}
    case FLOAT:  {std::vector<float> tmp; convert(val, tmp); nc_put_att_float(groupId, NC_GLOBAL, name.c_str(), NC_FLOAT, tmp.size(), tmp.data());  break;}
    case DOUBLE: {nc_put_att_double(groupId, NC_GLOBAL, name.c_str(), NC_DOUBLE, val.size(), val.field()); break;}
  }
  return Attribute(groupId, NC_GLOBAL, name);
}

/***********************************************/

std::vector<NetCdf::Attribute> NetCdf::Group::attributes() const
{
  Int countAttributes = 0;
  nc_inq_natts(groupId, &countAttributes);

  std::vector<Attribute> atts;
  for(Int attrId=0; attrId<countAttributes; attrId++)
  {
    char name[NC_MAX_NAME+1];
    nc_inq_attname(groupId, NC_GLOBAL, attrId, name);
    atts.push_back(Attribute(groupId, NC_GLOBAL, std::string(name)));
  }

  return atts;
}

/***********************************************/
/***********************************************/

std::string NetCdf::Dimension::name() const
{
  char dimName[NC_MAX_NAME+1];
  nc_inq_dimname(groupId, dimId, dimName);
  return std::string(dimName);
}

/***********************************************/

UInt NetCdf::Dimension::length() const
{
  UInt len = 0;
  nc_inq_dimlen(groupId, dimId, &len);
  return len;
}

/***********************************************/
/***********************************************/

std::string NetCdf::Variable::name() const
{
  char varName[NC_MAX_NAME+1];
  nc_inq_varname(groupId, varId, varName);
  return std::string(varName);
}

/***********************************************/

std::vector<NetCdf::Dimension> NetCdf::Variable::dimensions() const
{
  Int ndims = 0;
  nc_inq_varndims(groupId, varId, &ndims);

  std::vector<Int> dimIds(ndims);
  nc_inq_vardimid(groupId, varId, dimIds.data());

  std::vector<Dimension> dims;
  for(Int dimId : dimIds)
    dims.push_back(NetCdf::Dimension(groupId, dimId));
  return dims;
}

/***********************************************/

NetCdf::Attribute NetCdf::Variable::addAttribute(const std::string &name, const std::string &value)
{
  nc_put_att_text(groupId, varId, name.c_str(), value.size(), value.c_str());
  return Attribute(groupId, varId, name);
}

/***********************************************/

NetCdf::Attribute NetCdf::Variable::addAttribute(const std::string &name, DataType dtype, const Vector &val)
{
  switch(dtype)
  {
    case BYTE:   {std::vector<Byte>  tmp; convert(val, tmp); nc_put_att      (groupId, varId, name.c_str(), NC_BYTE,  tmp.size(), reinterpret_cast<void*>(tmp.data())); break;}
    case CHAR:   {std::vector<char>  tmp; convert(val, tmp); nc_put_att_text (groupId, varId, name.c_str(),           tmp.size(), tmp.data());  break;}
    case SHORT:  {std::vector<short> tmp; convert(val, tmp); nc_put_att_short(groupId, varId, name.c_str(), NC_SHORT, tmp.size(), tmp.data());  break;}
    case INT:    {std::vector<int>   tmp; convert(val, tmp); nc_put_att_int  (groupId, varId, name.c_str(), NC_INT,   tmp.size(), tmp.data());  break;}
    case FLOAT:  {std::vector<float> tmp; convert(val, tmp); nc_put_att_float(groupId, varId, name.c_str(), NC_FLOAT, tmp.size(), tmp.data());  break;}
    case DOUBLE: {nc_put_att_double(groupId, varId, name.c_str(), NC_DOUBLE, val.size(), val.field()); break;}
  }
  return Attribute(groupId, varId, name);
}

/***********************************************/

std::vector<NetCdf::Attribute> NetCdf::Variable::attributes() const
{
  Int countAttributes = 0;
  nc_inq_varnatts(groupId, varId, &countAttributes);

  std::vector<Attribute> atts;
  for(Int attrId=0; attrId<countAttributes; attrId++)
  {
    char name[NC_MAX_NAME+1];
    nc_inq_attname(groupId, varId, attrId, name);
    atts.push_back(Attribute(groupId, varId, std::string(name)));
  }

  return atts;
}

/***********************************************/

NetCdf::Attribute NetCdf::Variable::attribute(const std::string &name) const
{
  try
  {
    std::vector<Attribute> attr = attributes();
    auto iter = std::find_if(attr.begin(), attr.end(), [&name](const Attribute &attr) {return attr.name() == name;});
    if(iter == attr.end())
      throw(Exception("Attribute <"+name+"> not found in variable <"+this->name()+">"));
    return *iter;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NetCdf::Variable::setValues(const std::vector<UInt> &start, const std::vector<UInt> &count, const Vector &val) const
{
  nc_type ncType = 0;
  nc_inq_vartype(groupId, varId, &ncType);

  switch(ncType)
  {
    case NC_BYTE:   {std::vector<Byte>  tmp; convert(val, tmp); nc_put_vara      (groupId, varId, start.data(), count.data(), reinterpret_cast<void*>(tmp.data()));  break;}
    case NC_CHAR:   {std::vector<char>  tmp; convert(val, tmp); nc_put_vara_text (groupId, varId, start.data(), count.data(), tmp.data());  break;}
    case NC_SHORT:  {std::vector<short> tmp; convert(val, tmp); nc_put_vara_short(groupId, varId, start.data(), count.data(), tmp.data());  break;}
    case NC_INT:    {std::vector<int>   tmp; convert(val, tmp); nc_put_vara_int  (groupId, varId, start.data(), count.data(), tmp.data());  break;}
    case NC_FLOAT:  {std::vector<float> tmp; convert(val, tmp); nc_put_vara_float(groupId, varId, start.data(), count.data(), tmp.data());  break;}
    case NC_DOUBLE: {nc_put_vara_double(groupId, varId, start.data(), count.data(), val.field()); break;}
    default: throw(Exception("Unsupported data type."));
  }
}

/***********************************************/

void NetCdf::Variable::setValues(const Vector &val) const
{
  std::vector<Dimension> dims = dimensions();
  std::vector<UInt> count;
  for(auto &dim : dims)
    count.push_back(dim.length());
  setValues(std::vector<UInt>(dims.size(), 0), count, val);
}

/***********************************************/

Vector NetCdf::Variable::values(const std::vector<UInt> &start, const std::vector<UInt> &count) const
{
  try
  {
    UInt countElements = 1;
    for(UInt c : count)
      countElements *= c;

    nc_type ncType = 0;
    nc_inq_vartype(groupId, varId, &ncType);

    switch(ncType)
    {
      case NC_BYTE:   {std::vector<Byte>  tmp(countElements); nc_get_vara       (groupId, varId, start.data(), count.data(), reinterpret_cast<void*>(tmp.data())); return convert(tmp);}
      case NC_CHAR:   {std::vector<char>  tmp(countElements); nc_get_vara_text  (groupId, varId, start.data(), count.data(), tmp.data());  return convert(tmp);}
      case NC_SHORT:  {std::vector<short> tmp(countElements); nc_get_vara_short (groupId, varId, start.data(), count.data(), tmp.data());  return convert(tmp);}
      case NC_INT:    {std::vector<int>   tmp(countElements); nc_get_vara_int   (groupId, varId, start.data(), count.data(), tmp.data());  return convert(tmp);}
      case NC_FLOAT:  {std::vector<float> tmp(countElements); nc_get_vara_float (groupId, varId, start.data(), count.data(), tmp.data());  return convert(tmp);}
      case NC_DOUBLE: {Vector             tmp(countElements); nc_get_vara_double(groupId, varId, start.data(), count.data(), tmp.field()); return tmp;}
      default: throw(Exception("Unsupported data type."));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NetCdf::Variable::values() const
{
  std::vector<Dimension> dims = dimensions();
  std::vector<UInt> count;
  for(auto &dim : dims)
    count.push_back(dim.length());
  return values(std::vector<UInt>(dims.size(), 0), count);
}

/***********************************************/
/***********************************************/

std::string NetCdf::Attribute::value() const
{
  try
  {
    UInt len = 0;
    nc_inq_attlen(groupId, varId, name_.c_str(), &len);
    if(len == 0)
      return std::string();

    nc_type ncType = 0;
    nc_inq_atttype(groupId, varId, name_.c_str(), &ncType);

    Vector data;
    switch(ncType)
    {
      case NC_CHAR:   {std::vector<char>   tmp(len); nc_get_att_text  (groupId, varId, name_.c_str(), tmp.data()); return std::string(tmp.begin(), tmp.end());}
      case NC_SHORT:  {std::vector<short>  tmp(len); nc_get_att_short (groupId, varId, name_.c_str(), tmp.data()); data = convert(tmp); break;}
      case NC_INT:    {std::vector<int>    tmp(len); nc_get_att_int   (groupId, varId, name_.c_str(), tmp.data()); data = convert(tmp); break;}
      case NC_FLOAT:  {std::vector<float>  tmp(len); nc_get_att_float (groupId, varId, name_.c_str(), tmp.data()); data = convert(tmp); break;}
      case NC_DOUBLE: {std::vector<double> tmp(len); nc_get_att_double(groupId, varId, name_.c_str(), tmp.data()); data = convert(tmp); break;}
      case NC_STRING:
      {
        std::size_t attlen = 0;
        nc_inq_attlen(groupId, varId, name_.c_str(), &attlen);
        std::vector<char*> strings(attlen, nullptr);
        nc_get_att_string(groupId, varId, name_.c_str(), strings.data());
        std::stringstream ss;
        ss<<strings.at(0);
        for(UInt i=1; i<data.rows(); i++)
          ss<<", "<<strings.at(i);
        nc_free_string(attlen, strings.data());
        return ss.str();
      }
      default:
        return "Unsupported data type.";
    }

    std::stringstream ss;
    ss<<data(0);
    for(UInt i=1; i<data.rows(); i++)
      ss<<", "<<data(i);
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

NetCdf::InFile::InFile(const FileName &fileName)
{
  if(nc_open(fileName.c_str(), NC_NOWRITE, &groupId) != NC_NOERR)
    throw(Exception("Error opening NetCDF file <"+fileName.str()+">"));
}

/***********************************************/

NetCdf::InFile::~InFile()
{
  nc_close(groupId);
}

/***********************************************/

NetCdf::OutFile::OutFile(const FileName &fileName)
{
  if(nc_create(fileName.c_str(), NC_NETCDF4 | NC_SHARE, &groupId) != NC_NOERR)
    throw(Exception("Error opening NetCDF file <"+fileName.str()+">"));
}

/***********************************************/

NetCdf::OutFile::~OutFile()
{
  nc_close(groupId);
}

/***********************************************/
/***********************************************/

std::vector<Time> NetCdf::convertTimes(const_MatrixSliceRef values, const std::string &unitString)
{
  try
  {
// if you want to use uduinits2 library instead:
//   ut_set_error_message_handler(ut_ignore);
//   ut_system    *usys       = ut_read_xml(nullptr);
//   ut_unit      *timeUnit   = ut_parse(usys, unitString.c_str(), UT_ASCII);
//   ut_unit      *targetUnit = ut_parse(usys, "days since 1858-11-17 00:00:00", UT_ASCII);
//   cv_converter *converter  = ut_get_converter(timeUnit, targetUnit);
//   std::vector<Time> times;
//   times.reserve(values.size());
//   for(UInt i=0; i<values.size(); i++)
//     times.push_back(mjd2time(cv_convert_double(converter, values.at(i))));
//   cv_free(converter);
//   ut_free(targetUnit);
//   ut_free(timeUnit);
//   ut_free_system(usys);
//   return times;

    std::vector<std::string> parts;
    for(const std::string &part : String::split(unitString, " \t\n"))
      if(!part.empty())
        parts.push_back(String::lowerCase(part));

    if((parts.size() < 3) || ((parts.at(1) != "since") && (parts.at(1) != "from") && (parts.at(1) != "ref")))
      throw(Exception("Cannot interpret time units: '"+unitString+"'"));

    // time unit
    // ---------
    Double factor = 1.;
    if(String::startsWith(parts.at(0), "ms"))      factor = 1e-3/(24*60*60);          // milliseconds
    else if(String::startsWith(parts.at(0), "mo")) factor = 3.15569259747e7/86400/12; // month
    else if(String::startsWith(parts.at(0), "y"))  factor = 3.15569259747e7/86400;    // year
    else if(String::startsWith(parts.at(0), "s"))  factor = 1./(24*60*60);            // seconds
    else if(String::startsWith(parts.at(0), "m"))  factor = 1./(24*60);               // minutes
    else if(String::startsWith(parts.at(0), "h"))  factor = 1./24.;                   // hours
    else if(String::startsWith(parts.at(0), "d"))  factor = 1.;                       // days
    else if(String::startsWith(parts.at(0), "w"))  factor = 7.;                       // weeks
    else throw(Exception("Cannot interpret units: '"+unitString+"'"));

    // reference time
    // --------------
    UInt   year=0, month=1, day=1, hour=0, minute=0;
    Double second=0;

    std::vector<std::string> dateParts = String::split(parts.at(2), 't'); // format yyyy-MM-ddTHH:mm:ss.SSS
    if(dateParts.size() > 1)
      parts.insert(parts.begin()+3, dateParts.at(1));
    // date
    if(dateParts.at(0).find('-') == std::string::npos)
    {
      year  = static_cast<UInt>(String::toInt(dateParts.at(0).substr(0, 4)));
      if(dateParts.at(0).size() > 4) month = static_cast<UInt>(String::toInt(dateParts.at(0).substr(4, 2)));
      if(dateParts.at(0).size() > 6) day   = static_cast<UInt>(String::toInt(dateParts.at(0).substr(6, 2)));
    }
    else
    {
      Char c;
      std::stringstream ss(parts.at(2));
      ss>>year>>c>>month>>c>>day;
    }
    // time
    if(parts.size() > 3)
    {
      if(parts.at(3).find(':') == std::string::npos)
      {
        hour   = static_cast<UInt>(String::toInt(parts.at(3).substr(0, 4)));
        if(parts.at(3).size() > 4) minute = static_cast<UInt>(String::toInt(parts.at(3).substr(4, 2)));
        if(parts.at(3).size() > 6) second = String::toDouble(parts.at(3).substr(6));
      }
      else
      {
        Char c;
        std::stringstream ss(parts.at(3));
        ss>>hour>>c>>minute>>c>>second;
      }
    }
    Time time0 = date2time(year, month, day, hour, minute, second);

    std::vector<Time> times;
    times.reserve(values.size());
    for(UInt i=0; i<values.size(); i++)
      times.push_back(mjd2time(factor*values(i,0)) + time0);

    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Angle> NetCdf::convertAngles(const_MatrixSliceRef values)
{
  std::vector<Angle> lonLat;
  lonLat.reserve(values.size());
  for(UInt i=0; i<values.size(); i++)
    lonLat.push_back(Angle(DEG2RAD * values(i,0)));
  return lonLat;
}

/***********************************************/
/***********************************************/

#endif
