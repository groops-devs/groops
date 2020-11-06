/***********************************************/
/**
* @file archiveAscii.cpp
*
* @brief Write/read archives in ASCII format.
*
* @author Torsten Mayer-Guerr
* @date 2009-11-12
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/doodson.h"
#include "base/sphericalHarmonics.h"
#include "base/gnssType.h"
#include "archive.h"
#include "archiveAscii.h"

/***********************************************/

OutArchiveAscii::OutArchiveAscii(std::ostream &_stream, const std::string &type, UInt version) : stream(_stream)
{
  i_width = 10;
  isNewLine = TRUE;

  stream<<"groops";
  if(!type.empty())
    stream<<" "<<type;
  if(version)
    stream<<" version="<<version;
  stream<<std::endl;
}

/***********************************************/

void OutArchiveAscii::endLine()
{
  if(!isNewLine)
    stream<<std::endl;
  isNewLine = TRUE;
}

/***********************************************/

void OutArchiveAscii::comment(const std::string &text)
{
  if(text.empty())
    return;
  endLine();
  stream<<"# "<<text<<std::endl;
}

/***********************************************/

void OutArchiveAscii::saveDouble(Double x, Int width, Int precision, Bool science)
{
  if(science)
    stream.setf(std::ios::scientific,std::ios::floatfield);
  else
    stream.setf(std::ios::fixed,std::ios::floatfield);
  stream.width(width);
  stream.precision(precision);
  stream<<x;
}

/***********************************************/
/***********************************************/

InArchiveAscii::InArchiveAscii(std::istream &_stream) : stream(_stream), _version(0)
{
  stripComments();

  char c;
  stream>>c;
  stream.putback(c);
  if(c!='g')
    return;

  // headerline
  std::string line;
  std::getline(stream, line);
  if(stream.fail() || line.empty())
    return;
  std::stringstream ss(line);
  std::string text;
  ss>>text;
  if(text!="groops")
    throw(Exception("in InArchiveAscii: 'groops' expected"));

  // type
  ss>>text;
  if(ss.fail())
    return;
  if(text.find("version=")==std::string::npos)
  {
    typeStr = text;
    ss>>text;
    if(ss.fail())
      return;
  }

  // text contains now version
  if(text.find("version=")==std::string::npos)
    throw(Exception("in InArchiveAscii: 'version=' expected:"+text));
  _version = versionStr2version(text.substr(8));
}

/***********************************************/

void InArchiveAscii::stripComments()
{
  try
  {
    char c;
    stream>>c;
    stream.putback(c);
    if(c!='#')
      return;
    std::string s;
    std::getline(stream, s); // skip rest of line
    stripComments();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double InArchiveAscii::readDouble(std::istream &stream_)
{
  try
  {
    std::string dummy;
    stream_>>dummy;

    auto dpos = dummy.find_first_of("Dd");
    if(dpos != std::string::npos)
      dummy[dpos] = 'e';

    try
    {
      return std::stod(dummy);
    }
    catch(std::invalid_argument &)
    {
      throw(Exception("cannot read number: "+dummy));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void OutArchiveAscii::save(const Int      &x) {stream<<' '<<std::setw(i_width)<<x; isNewLine = FALSE;}
void OutArchiveAscii::save(const UInt     &x) {stream<<' '<<std::setw(i_width)<<x; isNewLine = FALSE;}
void OutArchiveAscii::save(const Double   &x) {saveDouble(x);                      isNewLine = FALSE;}
void OutArchiveAscii::save(const Bool     &x) {stream<<' '<<x;                     isNewLine = FALSE;}
void OutArchiveAscii::save(const Angle    &x) {saveDouble(x*RAD2DEG);              isNewLine = FALSE;}
void OutArchiveAscii::save(const Doodson  &x) {stream<<' '<<x.code();              isNewLine = FALSE;}
void OutArchiveAscii::save(const GnssType &x) {stream<<' '<<x.str();               isNewLine = FALSE;}

/***********************************************/

void InArchiveAscii::load(Int      &x) {stripComments(); stream>>x;}
void InArchiveAscii::load(UInt     &x) {stripComments(); stream>>x;}
void InArchiveAscii::load(Double   &x) {stripComments(); x = readDouble(stream);}
void InArchiveAscii::load(Bool     &x) {stripComments(); stream>>x;}
void InArchiveAscii::load(Angle    &x) {stripComments(); Double w = readDouble(stream); x = Angle(w*DEG2RAD);}
void InArchiveAscii::load(Doodson  &x) {stripComments(); std::string s; stream>>s; x = Doodson(s);}
void InArchiveAscii::load(GnssType &x) {stripComments(); std::string s; stream>>s; x = GnssType(s);}

/***********************************************/
/***********************************************/

void OutArchiveAscii::save(const std::string &x)
{
  try
  {
    stream<<' ';
    if((!x.empty()) && (x.find_first_of(" \t\n#") == std::string::npos))
      stream<<x; // string without special characters
    else
      stream<<'"'<<x<<'"';
    isNewLine = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void InArchiveAscii::load(std::string &x)
{
  try
  {
    stripComments();
    stream>>x;
    if(x.at(0)=='"')
    {
      if(x.size()==1)
        x += stream.get();
      while(x.at(x.size()-1)!='"')
        x += stream.get();
      x = x.substr(1, x.size()-2);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void OutArchiveAscii::save(const Time &x)
{
  try
  {
    isNewLine = FALSE;

    LongDouble mjd = x.mjdMod();
    mjd += x.mjdInt();

    stream<<' ';
    stream.setf(std::ios::fixed,std::ios::floatfield);
    stream.width(5+1+18);
    stream.precision(18);
    stream<<mjd;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InArchiveAscii::load(Time &x)
{
  try
  {
    stripComments();
    LongDouble mjd;
    stream>>mjd;

    Int mjdInt = static_cast<Int>(std::floor(mjd));
    x = Time(mjdInt, static_cast<Double>(mjd-mjdInt));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void OutArchiveAscii::save(const Vector &x)
{
  Matrix A = x;
  save(A);
}

/***********************************************/

void InArchiveAscii::load(Vector &x)
{
  Matrix A;
  load(A);
  x = A;
}

/***********************************************/
/***********************************************/

void OutArchiveAscii::save(const const_MatrixSlice &A)
{
  try
  {
    endLine();

    if(A.getType()==Matrix::GENERAL)
    {
      stream<<"Matrix( "<<A.rows()<<" x "<<A.columns()<<" )"<<std::endl;
      for(UInt i=0; i<A.rows(); i++)
      {
        for(UInt k=0; k<A.columns(); k++)
          saveDouble(A(i,k));
        stream<<std::endl;
      }
      stream<<std::endl;
      return;
    }

    if(A.isUpper())
      stream<<"Upper";
    else
      stream<<"Lower";
    if(A.getType()==Matrix::TRIANGULAR)
      stream<<"TriangularMatrix";
    else if(A.getType()==Matrix::SYMMETRIC)
      stream<<"SymmetricMatrix";
    stream<<"( "<<A.rows()<<" x "<<A.columns()<<" )"<<std::endl;

    if(A.isUpper())
    {
      for(UInt i=0; i<A.rows(); i++)
      {
        for(UInt k=i; k<A.columns(); k++)
          saveDouble(A(i,k));
        stream<<std::endl;
      }
    }
    else
    {
      for(UInt i=0; i<A.rows(); i++)
      {
        for(UInt k=0; k<=i; k++)
          saveDouble(A(i,k));
        stream<<std::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InArchiveAscii::load(Matrix &A)
{
  try
  {
    stripComments();

    char c;
    stream>>c;
    stream.putback(c);
    if(!isalpha(c) || (c == 'n') || (c == 'N')) // check for Nan
    {
      std::vector< std::vector<Double> > values;
      for(UInt i=0; ; i++)
      {
        std::string line;
        try
        {
          stripComments();
          std::getline(stream, line);
        }
        catch(std::exception &e)
        {
          if(stream.eof())
            break;
          GROOPS_RETHROW(e)
        }
        if(line.empty())
          break;
        values.resize(i+1);
        std::stringstream ss(line);
        for(;;)
        {
          std::string dummy;
          ss>>dummy;
          if(ss.fail() || ss.bad() || (dummy.at(0) == '#'))
            break;
          values.at(i).push_back(std::stod(dummy));
        }
      }
      A = Matrix(values.size(), values.at(0).size());
      for(UInt i=0; i<A.rows(); i++)
        for(UInt k=0; k<A.columns(); k++)
          A(i,k) = values.at(i).at(k);
      return;
    }

    std::string type;
    UInt rows, columns;
    stream>>type>>rows>>c>>columns>>c;
    if(type=="Matrix(")
    {
      A = Matrix(rows,columns);
      for(UInt i=0; i<A.rows(); i++)
        for(UInt k=0; k<A.columns(); k++)
          A(i,k) = readDouble(stream);
      return;
    }

    if(type=="UpperTriangularMatrix(")
      A = Matrix(rows, Matrix::TRIANGULAR, Matrix::UPPER);
    else if(type=="LowerTriangularMatrix(")
      A = Matrix(rows, Matrix::TRIANGULAR, Matrix::LOWER);
    else if(type=="UpperSymmetricMatrix(")
      A = Matrix(rows, Matrix::SYMMETRIC, Matrix::UPPER);
    else if(type=="LowerSymmetricMatrix(")
      A = Matrix(rows, Matrix::SYMMETRIC, Matrix::LOWER);
    else
      throw(Exception("unknown matrix type: "+type));

    if(A.isUpper())
      for(UInt i=0; i<A.rows(); i++)
        for(UInt k=i; k<A.columns(); k++)
          A(i,k) = readDouble(stream);
    else
      for(UInt i=0; i<A.rows(); i++)
        for(UInt k=0; k<=i; k++)
          A(i,k) = readDouble(stream);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void OutArchiveAscii::save(const SphericalHarmonics &harm)
{
  try
  {
    endLine();

    const Bool hasErrors = harm.sigma2cnm().size();

    stream<<"modelname              GROOPS"<<std::endl;
    stream<<"product_type           gravity_field"<<std::endl;
    stream<<"earth_gravity_constant "; saveDouble(harm.GM(),16,10,TRUE); stream<<std::endl;
    stream<<"radius                 "; saveDouble(harm.R(), 16,10,TRUE); stream<<std::endl;
    stream<<"max_degree             "<<std::setw(0)<<harm.maxDegree()<<std::endl;
    stream<<"norm                   fully_normalized"<<std::endl;
    //stream<<"tide_system            zero_tide"<<std::endl;
    stream<<"errors                 "<<(hasErrors ? "formal" : "no")<<std::endl;
    stream<<std::endl;
    if(harm.sigma2cnm().size())
    {
      stream<<"key   L   M       C                   S                   sigma C             sigma S"<<std::endl;
      stream<<"end_of_head ================================================================================"<<std::endl;
    }
    else
    {
      stream<<"key   L   M       C                   S"<<std::endl;
      stream<<"end_of_head =========================================="<<std::endl;
    }
    for(UInt n=0; n<=harm.maxDegree(); n++)
      for(UInt m=0; m<=n; m++)
      {
        // RR: changed setw(4) to setw(5)
        stream<<"gfc "<<std::setw(5)<<n<<std::setw(5)<<m;
        saveDouble(harm.cnm()(n,m),20,12,TRUE);
        saveDouble(harm.snm()(n,m),20,12,TRUE);
        if(hasErrors) saveDouble(sqrt(harm.sigma2cnm()(n,m)),20,12,TRUE);
        if(hasErrors) saveDouble(sqrt(harm.sigma2snm()(n,m)),20,12,TRUE);
        stream<<std::endl;
      }
    stream<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InArchiveAscii::load(SphericalHarmonics &harm)
{
  try
  {
    stripComments();

    Double GM=0, R=0;
    UInt   maxDegree=0;
    Bool   hasErrors = FALSE;

    // Header einlesen
    std::string line;
    for(;;)
    {
      std::getline(stream, line);
      std::stringstream ss(line);
      std::string tag;
      ss>>tag;

      if(tag == "earth_gravity_constant")
        GM = readDouble(ss);
      if(tag == "radius")
        R = readDouble(ss);
      if(tag == "max_degree")
        ss>>maxDegree;
      if(tag == "errors")
      {
        std::string no;
        ss>>no;
        if(no != "no")
          hasErrors = TRUE;
      }
      // Header fertig ?
      if(tag == "end_of_head")
        break;
    }

    Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    for(;;)
    {
      try
      {
        std::getline(stream, line);
      }
      catch(std::exception &e)
      {
        if(stream.eof())
          break;
        GROOPS_RETHROW(e)
      }
      if(line.empty())
        break;
      std::stringstream ss(line);
      ss.exceptions(std::ios::badbit|std::ios::failbit);
      std::string tag;
      ss>>tag;
      if((tag != "gfc")&&(tag != "gfct"))
        continue;
      UInt n,m;
      ss>>n>>m;
      cnm(n,m) = readDouble(ss);
      snm(n,m) = readDouble(ss);
      if(hasErrors)
      {
        sigma2cnm(n,m) = std::pow(readDouble(ss), 2);
        sigma2snm(n,m) = std::pow(readDouble(ss), 2);
      }
    }

    if(hasErrors)
      harm = SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm);
    else
      harm = SphericalHarmonics(GM, R, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
