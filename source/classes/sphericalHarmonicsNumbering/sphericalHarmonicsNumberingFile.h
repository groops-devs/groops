/***********************************************/
/**
* @file sphericalHarmonicsNumberingFile.h
*
* @brief Numbering schema of spherical harmonics coefficients.
* Numbering as specified in the chosen file.
* E.g. the kite numbering scheme.
* @see SphericalHarmonicsNumbering
*
* @author Torsten Mayer-Guerr
* @date 2017-12-20
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSNUMBERINGFILE__
#define __GROOPS_SPHERICALHARMONICSNUMBERINGFILE__

#ifdef DOCSTRING_SphericalHarmonicsNumbering
static const char *docstringSphericalHarmonicsNumberingFile = R"(
\subsection{File}
Numbering as specified in the chosen file.
The \configFile{inputfile}{matrix} is a matrix with the first column indicating cnm/snm with 0 or 1.
The second and third column specify degree and order.
)";
#endif

/***********************************************/

#include "files/fileMatrix.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Numbering schema of spherical harmonics coefficients.
* Numbering as specified in the chosen file.
* E.g. the kite numbering scheme.
* @see SphericalHarmonicsNumbering */
class SphericalHarmonicsNumberingFile : public SphericalHarmonicsNumbering
{
private:
  std::vector<UInt> degree, order, cs;
  UInt maxDegree, minDegree;

public:
  SphericalHarmonicsNumberingFile(Config &config);

  virtual UInt parameterCount(UInt maxDegree, UInt minDegree) const override;
  virtual void numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const override;
};

/***********************************************/

inline SphericalHarmonicsNumberingFile::SphericalHarmonicsNumberingFile(Config &config)
{
  try
  {
    FileName inName;
    readConfig(config, "inputfile", inName, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    Matrix A;
    readFileMatrix(inName, A);

    cs.resize(A.rows());
    degree.resize(A.rows());
    order.resize(A.rows());
    for(UInt i=0; i<A.rows(); i++)
    {
      if(!((cs.at(i) == 0) || (cs.at(i) == 1)))
        throw(Exception("First column must contain 0 or 1 for cnm/snm"));
      if(order.at(i) > degree.at(i))
        throw(Exception("order must not greater than degree in column "+i%"%i"s));
      cs.at(i)     = static_cast<UInt>(A(i,0));
      degree.at(i) = static_cast<UInt>(A(i,1));
      order.at(i)  = static_cast<UInt>(A(i,2));
    }

    minDegree = *std::min_element(degree.begin(), degree.end());
    maxDegree = *std::max_element(degree.begin(), degree.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt SphericalHarmonicsNumberingFile::parameterCount(UInt maxDegree, UInt minDegree) const
{
  return std::count_if(degree.begin(), degree.end(), [&](UInt n) {return (n <= maxDegree) && (n >= minDegree);});
}

/***********************************************/

void SphericalHarmonicsNumberingFile::numbering(UInt maxDegree, UInt minDegree, std::vector<std::vector<UInt>> &Cnm, std::vector<std::vector<UInt>> &Snm) const
{
  try
  {
    Cnm.clear(); Cnm.resize(maxDegree+1);
    Snm.clear(); Snm.resize(maxDegree+1);
    for(UInt n=0; n<=maxDegree; n++)
    {
      Cnm[n].resize(n+1, NULLINDEX);
      Snm[n].resize(n+1, NULLINDEX);
    }

    UInt idx = 0;
    for(UInt i=0; i<degree.size(); i++)
      if((degree[i] <= maxDegree) && (degree[i] >= minDegree))
      {
        if(cs[i] == 0)
          Cnm[degree[i]][order[i]] = idx++;
        else
          Snm[degree[i]][order[i]] = idx++;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
