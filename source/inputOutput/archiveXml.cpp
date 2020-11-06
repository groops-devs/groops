/***********************************************/
/**
* @file archiveXml.cpp
*
* @brief Read/write archive files in XML format.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-01
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/doodson.h"
#include "base/sphericalHarmonics.h"
#include "base/gnssType.h"
#include "parser/xml.h"
#include "archive.h"
#include "archiveXml.h"

/***********************************************/

OutArchiveXml::OutArchiveXml(std::ostream &_stream, const std::string &type, UInt version) : stream(_stream)
{
  XmlNodePtr xmlNode = XmlNode::create("groops");
  if(!type.empty())
    writeAttribute(xmlNode, "type", type);
  if(version)
    writeAttribute(xmlNode, "version", version);
  stack.push(xmlNode);
}

/***********************************************/

OutArchiveXml::~OutArchiveXml()
{
  XmlNodePtr xmlNode = stack.top();
  XmlNode::write(stream, xmlNode);
}

/***********************************************/

void OutArchiveXml::startTag(const std::string &name)
{
  stack.push(createXmlNode(stack.top(), name));
}

/***********************************************/

void OutArchiveXml::endTag(const std::string &/*name*/)
{
  stack.top() = XmlNodePtr();
  stack.pop();
}

/***********************************************/
/***********************************************/

InArchiveXml::InArchiveXml(std::istream &stream) : _version(0)
{
  XmlNodePtr xmlNode = XmlNode::read(stream);

  XmlAttrPtr attr = xmlNode->getAttribute("type");
  if(attr)
    typeStr = attr->getText();

  attr = xmlNode->getAttribute("version");
  if(attr)
    _version = versionStr2version(attr->getText());

  stack.push(xmlNode);
}

/***********************************************/

InArchiveXml::~InArchiveXml()
{
  while(!stack.empty())
    stack.pop();
}

/***********************************************/

void InArchiveXml::startTag(const std::string &name)
{
  stack.push(getChild(stack.top(), name, TRUE));
}

/***********************************************/

void InArchiveXml::endTag(const std::string &/*name*/)
{
  stack.pop();
}

/***********************************************/
/***********************************************/

void OutArchiveXml::save(const Int         &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const UInt        &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const Double      &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const Bool        &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const std::string &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const Angle       &x) {stack.top()->setValue(x);}
void OutArchiveXml::save(const Doodson     &x) {stack.top()->setValue(x.code());}
void OutArchiveXml::save(const GnssType    &x) {stack.top()->setValue(x.str());}

/***********************************************/

void InArchiveXml::load(Int         &x) {stack.top()->getValue(x);}
void InArchiveXml::load(UInt        &x) {stack.top()->getValue(x);}
void InArchiveXml::load(Double      &x) {stack.top()->getValue(x);}
void InArchiveXml::load(Bool        &x) {stack.top()->getValue(x);}
void InArchiveXml::load(std::string &x) {stack.top()->getValue(x);}
void InArchiveXml::load(Angle       &x) {stack.top()->getValue(x);}
void InArchiveXml::load(Doodson     &x) {std::string str; stack.top()->getValue(str); x = Doodson(str);}
void InArchiveXml::load(GnssType    &x) {std::string str; stack.top()->getValue(str); x = GnssType(str);}

/***********************************************/
/***********************************************/

void OutArchiveXml::save(const Time &x)
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

void InArchiveXml::load(Time &x)
{
  LongDouble mjd;
  stack.top()->getValue(mjd);
  Int mjdInt = static_cast<Int>(std::floor(mjd));
  x = Time(mjdInt, static_cast<Double>(mjd-mjdInt));
}

/***********************************************/
/***********************************************/

void OutArchiveXml::save(const Vector &x)
{
  Matrix A = x;
  save(A);
}

/***********************************************/

void InArchiveXml::load(Vector &x)
{
  Matrix A;
  load(A);
  x = A;
}

/***********************************************/
/***********************************************/

void OutArchiveXml::save(const const_MatrixSlice &x)
{
  XmlNodePtr xmlNode = stack.top();
  writeXml(xmlNode, "type", x.getType());

  if(x.getType()==Matrix::GENERAL)
  {
    writeXml(xmlNode, "rows",    x.rows());
    writeXml(xmlNode, "columns", x.columns());
    for(UInt s=0; s<x.columns(); s++)
      for(UInt z=0; z<x.rows(); z++)
      {
        if(x(z,s)==0.0) continue;
        XmlNodePtr xmlNode2 = writeXml(xmlNode, "cell", x(z,s));
        writeAttribute(xmlNode2, "row", z);
        writeAttribute(xmlNode2, "col", s);
      }
  }
  else
  {
    writeXml(xmlNode, "upper", x.isUpper());
    writeXml(xmlNode, "dimension", x.rows());
    if(x.isUpper())
    {
      for(UInt s=0; s<x.columns(); s++)
        for(UInt z=0; z<=s; z++)
        {
          if(x(z,s)==0.0) continue;
          XmlNodePtr xmlNode2 = writeXml(xmlNode, "cell",x(z,s));
          writeAttribute(xmlNode2, "row", z);
          writeAttribute(xmlNode2, "col", s);
        }
    }
    else
    {
      for(UInt s=0; s<x.columns(); s++)
        for(UInt z=s; z<x.rows(); z++)
        {
          if(x(z,s)==0.0) continue;
          XmlNodePtr xmlNode2 = writeXml(xmlNode, "cell", x(z,s));
          writeAttribute(xmlNode2, "row", z);
          writeAttribute(xmlNode2, "col", s);
        }
    }
  }
}

/***********************************************/

void InArchiveXml::load(Matrix &x)
{
  XmlNodePtr xmlNode = stack.top();

  Int  type = Matrix::GENERAL;
  readXml(xmlNode, "type", type);

  if(type==Matrix::GENERAL)
  {
    UInt rows, columns = 1;
    readXml(xmlNode, "rows",    rows, TRUE);
    readXml(xmlNode, "columns", columns);
    x = Matrix(rows, columns);
  }
  else
  {
    UInt dimension;
    Bool upper;
    readXml(xmlNode, "upper", upper);
    readXml(xmlNode, "dimension", dimension);
    x = Matrix(dimension, static_cast<Matrix::Type>(type), ((upper) ? Matrix::UPPER : Matrix::LOWER));
  }

  UInt  count = childCount(xmlNode, "cell");
  for(UInt i=0; i<count; i++)
  {
    Double a;
    UInt n, m = 0;
    XmlNodePtr xmlNode2 = readXml(xmlNode, "cell", a, TRUE);
    readAttribute(xmlNode2, "row", n, TRUE);
    readAttribute(xmlNode2, "col", m, FALSE);
    x(n,m) = a;
  }
}

/***********************************************/
/***********************************************/

void OutArchiveXml::save(const SphericalHarmonics &harm)
{
  XmlNodePtr xmlNode = stack.top();
  writeXml(xmlNode, "GM", harm.GM());
  writeXml(xmlNode, "R", harm.R());
  writeXml(xmlNode, "maxDegree", harm.maxDegree());

  for(UInt n=0; n<=harm.maxDegree(); n++)
    for(UInt m=0; m<=n; m++)
    {
      if(harm.cnm()(n,m)!=0.0)
      {
        XmlNodePtr xmlNode2 = writeXml(xmlNode, "cnm", harm.cnm()(n,m));
        writeAttribute(xmlNode2, "degree", n);
        writeAttribute(xmlNode2, "order",  m);
      }
      if(harm.snm()(n,m)!=0.0)
      {
        XmlNodePtr xmlNode2 = writeXml(xmlNode, "snm", harm.snm()(n,m));
        writeAttribute(xmlNode2, "degree", n);
        writeAttribute(xmlNode2, "order",  m);
      }
      if(harm.sigma2cnm().size() && harm.sigma2cnm()(n,m)!=0.0)
      {
        XmlNodePtr xmlNode2 = writeXml(xmlNode, "sigmacnm", sqrt(harm.sigma2cnm()(n,m)));
        writeAttribute(xmlNode2, "degree", n);
        writeAttribute(xmlNode2, "order",  m);
      }
      if(harm.sigma2snm().size() && harm.sigma2snm()(n,m)!=0.0)
      {
        XmlNodePtr xmlNode2 = writeXml(xmlNode, "sigmasnm", sqrt(harm.sigma2snm()(n,m)));
        writeAttribute(xmlNode2, "degree", n);
        writeAttribute(xmlNode2, "order",  m);
      }
    }
}

/***********************************************/

void InArchiveXml::load(SphericalHarmonics &harm)
{
  XmlNodePtr xmlNode = stack.top();
  UInt   n,m;
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
  Bool hasSigma = FALSE;

  std::string name;
  while(xmlNode->hasChildren())
  {
    XmlNodePtr xmlNode2 = getNextChild(xmlNode, name);
    if(name=="cnm")
    {
      readAttribute(xmlNode2, "degree", n, TRUE);
      readAttribute(xmlNode2, "order",  m, TRUE);
      xmlNode2->getValue(cnm(n,m));
    }
    else if(name=="snm")
    {
      readAttribute(xmlNode2, "degree", n, TRUE);
      readAttribute(xmlNode2, "order",  m, TRUE);
      xmlNode2->getValue(snm(n,m));
    }
    else if(name=="sigmacnm")
    {
      readAttribute(xmlNode2, "degree", n, TRUE);
      readAttribute(xmlNode2, "order",  m, TRUE);
      xmlNode2->getValue(sigma2cnm(n,m));
      sigma2cnm(n,m) *= sigma2cnm(n,m);
      hasSigma = TRUE;
    }
    else if(name=="sigmasnm")
    {
      readAttribute(xmlNode2, "degree", n, TRUE);
      readAttribute(xmlNode2, "order",  m, TRUE);
      xmlNode2->getValue(sigma2snm(n,m));
      sigma2snm(n,m) *= sigma2snm(n,m);
      hasSigma = TRUE;
    }
  }

  if(hasSigma)
    harm = SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm);
  else
    harm = SphericalHarmonics(GM, R, cnm, snm);
}

/***********************************************/
