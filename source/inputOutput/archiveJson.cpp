/***********************************************/
/**
* @file archiveJson.cpp
*
* @brief Read/write archive files in JSON format.
*
* @author Torsten Mayer-Guerr
* @date 2025-01-22
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/string.h"
#include "base/doodson.h"
#include "base/sphericalHarmonics.h"
#include "base/gnssType.h"
#include "parser/xml.h"
#include "archive.h"
#include "archiveJson.h"

/***********************************************/

static std::string quoted(const std::string &str)
{
  if(str.empty())
    return str;
  std::stringstream ss;
  ss<<std::quoted(str);
  return ss.str();
}

/***********************************************/

static void writeJsonData(std::ostream &stream, XmlNodePtr xmlNode, UInt depth)
{
  // has children -> struct
  if(xmlNode->hasChildren())
  {
    stream<<"{"<<std::endl;
    while(xmlNode->hasChildren())
    {
      XmlNodePtr child = xmlNode->getNextChild();
      stream<<std::string(depth, '\t')<<std::quoted(child->getName())<<": ";
      // is array?
      if(xmlNode->findChild(child->getName()))
      {
        stream<<"[";
        while(child)
        {
          writeJsonData(stream, child, depth+1);
          child = xmlNode->getChild(child->getName());
          if(child)
            stream<<", ";
        }
        stream<<"]";
      }
      else
        writeJsonData(stream, child, depth+1);
      if(xmlNode->hasChildren())
        stream<<","<<std::endl;
    }
    stream<<std::endl<<std::string(depth-1, '\t')<<"}";
  }
  else if(!xmlNode->getText().empty()) // value
  {
    stream<<xmlNode->getText();
  }
  else // empty
    stream<<"null";
}

/***********************************************/

OutArchiveJson::OutArchiveJson(std::ostream &_stream, const std::string &type, UInt version) : stream(_stream)
{
  XmlNodePtr xmlNode = XmlNode::create("groops");
  writeXml(xmlNode, "fileType",    quoted(type));
  writeXml(xmlNode, "fileVersion", version);
  stack.push(xmlNode);
}

/***********************************************/

OutArchiveJson::~OutArchiveJson()
{
  XmlNodePtr xmlNode = stack.top();
  writeJsonData(stream, xmlNode, 1);
}

/***********************************************/

void OutArchiveJson::startTag(const std::string &name)
{
  stack.push(createXmlNode(stack.top(), name));
}

/***********************************************/

void OutArchiveJson::endTag(const std::string &/*name*/)
{
  stack.top() = XmlNodePtr();
  stack.pop();
}

/***********************************************/
/***********************************************/

static void readObject(std::istream &stream, XmlNodePtr parent);

static void readData(std::istream &stream, XmlNodePtr xmlNode)
{
  char c;
  stream>>c;
  if(c == '{')
  {
    do
    {
      readObject(stream, xmlNode);
      stream>>c;
    }
    while(c == ',');
    if(c != '}')
      throw(Exception("json parser error in '"+xmlNode->getName()+"': '}' expected"));
    return;
  }

  stream.putback(c);
  std::string value;
  if(c == '"')
  {
    stream>>std::quoted(value);
    xmlNode->setValue(value);
  }
  else
  {
    // read until delimiter (,]} or white space)
    while(std::string(",]} \n\t").find_first_of((c = stream.get())) == std::string::npos)
      value.push_back(c);
    stream.putback(c);
    if(value == "true")
      xmlNode->setValue(1);
    else if(value == "false")
      xmlNode->setValue(0);
    else if(value != "null")
      xmlNode->setValue(value);
  }
}

/***********************************************/

static void readObject(std::istream &stream, XmlNodePtr parent)
{
  std::string name;
  stream>>std::quoted(name);
  XmlNodePtr xmlNode = createXmlNode(parent, name);

  char c;
  stream>>c; // ":"
  if(c != ':')
    throw(Exception("json parser error in '"+name+"': ':' expected"));

  // content
  stream>>c;
  if(c == '[') // array
  {
    for(;;)
    {
      readData(stream, xmlNode);
      stream>>c;
      if(c != ',')
        break;
      xmlNode = createXmlNode(parent, name);
    }
    if(c != ']')
      throw(Exception("json parser error in '"+name+"': ']' expected"));
    return;
  }

  stream.putback(c);
  readData(stream, xmlNode);
}

/***********************************************/

InArchiveJson::InArchiveJson(std::istream &stream) : _version(0)
{
  XmlNodePtr xmlNode = XmlNode::create("groops");
  stack.push(xmlNode);
  readData(stream, xmlNode);
  readXml(xmlNode, "fileType",    typeStr);
  readXml(xmlNode, "fileVersion", _version);
}

/***********************************************/

InArchiveJson::~InArchiveJson()
{
  while(!stack.empty())
    stack.pop();
}

/***********************************************/

void InArchiveJson::startTag(const std::string &name)
{
  stack.push(getChild(stack.top(), name, TRUE));
}

/***********************************************/

void InArchiveJson::endTag(const std::string &/*name*/)
{
  stack.pop();
}

/***********************************************/
/***********************************************/

void OutArchiveJson::save(const Int         &x) {stack.top()->setValue(x);}
void OutArchiveJson::save(const UInt        &x) {stack.top()->setValue(x);}
void OutArchiveJson::save(const Double      &x) {stack.top()->setValue(x);}
void OutArchiveJson::save(const Bool        &x) {stack.top()->setValue(x ? "true" : "false");}
void OutArchiveJson::save(const std::string &x) {stack.top()->setValue(quoted(x));}
void OutArchiveJson::save(const Angle       &x) {stack.top()->setValue(x);}
void OutArchiveJson::save(const Doodson     &x) {stack.top()->setValue(quoted(x.code()));}
void OutArchiveJson::save(const GnssType    &x) {stack.top()->setValue(quoted(x.str()));}

/***********************************************/

void InArchiveJson::load(Int         &x) {stack.top()->getValue(x);}
void InArchiveJson::load(UInt        &x) {stack.top()->getValue(x);}
void InArchiveJson::load(Double      &x) {stack.top()->getValue(x);}
void InArchiveJson::load(Bool        &x) {stack.top()->getValue(x);}
void InArchiveJson::load(std::string &x) {stack.top()->getValue(x);}
void InArchiveJson::load(Angle       &x) {stack.top()->getValue(x);}
void InArchiveJson::load(Doodson     &x) {std::string str; stack.top()->getValue(str); x = Doodson(str);}
void InArchiveJson::load(GnssType    &x) {std::string str; stack.top()->getValue(str); x = GnssType(str);}

/***********************************************/
/***********************************************/

void OutArchiveJson::save(const Time &x)
{
  try
  {
    LongDouble mjd = x.mjdMod();
    mjd += x.mjdInt();

    std::stringstream stream_;
    stream_.setf(std::ios::fixed,std::ios::floatfield);
    stream_.precision(16);
    stream_<<mjd;

    stack.top()->setValue(stream_.str());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InArchiveJson::load(Time &x)
{
  LongDouble mjd;
  stack.top()->getValue(mjd);
  Int mjdInt = static_cast<Int>(std::floor(mjd));
  x = Time(mjdInt, static_cast<Double>(mjd-mjdInt));
}

/***********************************************/
/***********************************************/

void OutArchiveJson::save(const Vector &x)
{
  Matrix A = x;
  save(A);
}

/***********************************************/

void InArchiveJson::load(Vector &x)
{
  Matrix A;
  load(A);
  x = A;
}

/***********************************************/
/***********************************************/

void OutArchiveJson::save(const const_MatrixSlice &A)
{
  XmlNodePtr xmlNode = stack.top();

  std::string type;
  if(A.getType()==Matrix::GENERAL)
    type = "general";
  else
  {
    type = (A.isUpper() ? "upper" : "lower");
    if(A.getType()==Matrix::TRIANGULAR)
      type += "TriangularMatrix";
    else if(A.getType()==Matrix::SYMMETRIC)
      type += "SymmetricMatrix";
  }

  std::stringstream ss;
  ss.setf(std::ios::scientific,std::ios::floatfield);
  ss.precision(18);
  UInt minCol = 0, maxCol = A.columns();
  ss<<std::endl<<"[";
  for(UInt i=0; i<A.rows(); i++)
  {
    if((A.getType() != Matrix::GENERAL) &&  A.isUpper()) minCol = i;
    if((A.getType() != Matrix::GENERAL) && !A.isUpper()) maxCol = i+1;
    ss<<"[";
    for(UInt k=minCol; k<maxCol; k++)
      ss<<A(i,k)<<((k<maxCol-1) ? ", " : "]");
    ss<<((i < A.rows()-1) ? ",\n " : "]");
  }

  writeXml(xmlNode, "type",    quoted(type));
  writeXml(xmlNode, "rows",    A.rows());
  writeXml(xmlNode, "columns", A.columns());
  if(A.size())
    writeXml(xmlNode, "data",  ss.str());
  else
    writeXml(xmlNode, "data",  "[]");
}

/***********************************************/

void InArchiveJson::load(Matrix &A)
{
  XmlNodePtr xmlNode = stack.top();

  std::string type;
  UInt rows, columns = 1;
  readXml(xmlNode, "type",    type, TRUE);
  readXml(xmlNode, "rows",    rows, TRUE);
  readXml(xmlNode, "columns", columns);

  if(type == "general")
    A = Matrix(rows, columns, Matrix::NOFILL);
  else if(type == "upperTriangularMatrix")
    A = Matrix(rows, Matrix::TRIANGULAR, Matrix::UPPER);
  else if(type == "lowerTriangularMatrix")
    A = Matrix(rows, Matrix::TRIANGULAR, Matrix::LOWER);
  else if(type == "upperSymmetricMatrix")
    A = Matrix(rows, Matrix::SYMMETRIC, Matrix::UPPER);
  else if(type == "lowerSymmetricMatrix")
    A = Matrix(rows, Matrix::SYMMETRIC, Matrix::LOWER);
  else
    throw(Exception("unknown matrix type: "+type));

  UInt minCol = 0, maxCol = A.columns();
  for(UInt i=0; i<A.rows(); i++)
  {
    if((A.getType() != Matrix::GENERAL) &&  A.isUpper()) minCol = i;
    if((A.getType() != Matrix::GENERAL) && !A.isUpper()) maxCol = i+1;
    for(UInt k=minCol; k<maxCol; k++)
      readXml(xmlNode, "data", A(i,k), TRUE);
  }
}

/***********************************************/
/***********************************************/

void OutArchiveJson::save(const SphericalHarmonics &harm)
{
  XmlNodePtr xmlNode = stack.top();
  writeXml(xmlNode, "GM",        harm.GM());
  writeXml(xmlNode, "R",         harm.R());
  writeXml(xmlNode, "maxDegree", harm.maxDegree());
  for(UInt n=0; n<=harm.maxDegree(); n++)
    for(UInt m=0; m<=n; m++)
    {
      writeXml(xmlNode, "cnm", harm.cnm()(n,m));
      writeXml(xmlNode, "snm", harm.snm()(n,m));
      if(harm.sigma2cnm().size())
      {
        writeXml(xmlNode, "sigmacnm", std::sqrt(harm.sigma2cnm()(n,m)));
        writeXml(xmlNode, "sigmasnm", std::sqrt(harm.sigma2snm()(n,m)));
      }
    }
}

/***********************************************/

void InArchiveJson::load(SphericalHarmonics &harm)
{
  XmlNodePtr xmlNode = stack.top();
  UInt   maxDegree;
  Double GM = DEFAULT_GM;
  Double R  = DEFAULT_R;

  readXml(xmlNode, "GM", GM);
  readXml(xmlNode, "R",  R);
  readXml(xmlNode, "maxDegree", maxDegree, TRUE);

  Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
  Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

  for(UInt n=0; n<=maxDegree; n++)
    for(UInt m=0; m<=n; m++)
    {
      readXml(xmlNode, "cnm", cnm(n,m), TRUE);
      readXml(xmlNode, "snm", snm(n,m), TRUE);
      readXml(xmlNode, "sigmacnm", sigma2cnm(n,m));
      readXml(xmlNode, "sigmasnm", sigma2snm(n,m));
      sigma2cnm(n,m) *= sigma2cnm(n,m);
      sigma2snm(n,m) *= sigma2snm(n,m);
    }

  if(isStrictlyZero(sigma2cnm))
    harm = SphericalHarmonics(GM, R, cnm, snm);
  else
    harm = SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm);
}

/***********************************************/
