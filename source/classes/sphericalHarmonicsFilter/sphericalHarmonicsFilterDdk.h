/***********************************************/
/**
* @file sphericalHarmonicsFilterDdk.h
*
* @brief Smoothing by a DDK filter (order wise).
* @see SphericalHarmonicsFilter
*
* @author Andreas Kvas
* @date 2020-04-16
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSFILTERDDK__
#define __GROOPS_SPHERICALHARMONICSFILTERDDK__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsFilter
static const char *docstringSphericalHarmonicsFilterDdk = R"(
\subsection{DDK}
Orderwise filtering with the DDK filter by Kusche et al. 2009.
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Smoothing by a DDK filter (order wise).
* @ingroup sphericalHarmonicsFilterGroup
* @see SphericalHarmonicsFilter */
class SphericalHarmonicsFilterDdk : public SphericalHarmonicsFilterBase
{
  std::vector<Matrix> matrix;

public:
  SphericalHarmonicsFilterDdk(Config &config);

  SphericalHarmonics filter(const SphericalHarmonics &harm) const override;
};

/***********************************************/

inline SphericalHarmonicsFilterDdk::SphericalHarmonicsFilterDdk(Config &config)
{
  try
  {
    FileName inName;
    UInt     level;

    renameDeprecatedConfig(config, "inputfileNormalequation", "inputfileNormalEquation", date2time(2020, 6, 3));

    readConfig(config, "level",                   level,  Config::MUSTSET, "", "DDK filter level (1, 2, 3, ...)");
    readConfig(config, "inputfileNormalEquation", inName, Config::MUSTSET, "{groopsDataDir}/sphericalHarmonicsFilter/DDK/normalsKuscheGfzBlock_n2-120_orderwiseNonAlternating.dat.gz", "");
    if(isCreateSchema(config)) return;

    const Double factor = std::pow(10, 15-level);
    const Double power  = 4.0;

    Matrix N, n;
    NormalEquationInfo info;
    readFileNormalEquation(inName, info, N, n);
    fillSymmetric(N);

    std::vector<std::vector<UInt>> idxC, idxS;
    std::vector<UInt>              cs, orders, degrees;
    for(UInt idx=0; idx<info.parameterName.size(); idx++)
    {
      std::string parameter = info.parameterName.at(idx).type;
      auto pos = parameter.find('.');
      std::string baseName = parameter.substr(0, pos);
      if(baseName != "sphericalHarmonics")
        throw(Exception("Non-spherical harmonic parameter found."));

      std::istringstream is(parameter.substr(pos+1));
      std::vector<std::string> tokens;
      std::string part;
      while(std::getline(is, part, '_'))
        tokens.push_back(part);

      cs.push_back((tokens.at(0) == "c") ? 0 : 1);
      degrees.push_back(std::stoul(tokens.at(1)));
      orders.push_back(std::stoul(tokens.at(2)));
    }

    const UInt maxDegree = *std::max_element(degrees.begin(), degrees.end());
    idxC.resize(maxDegree+1);
    idxS.resize(maxDegree+1);
    for(UInt n=0; n<=maxDegree; n++)
    {
      idxC.at(n).resize(n+1, NULLINDEX);
      idxS.at(n).resize(n+1, NULLINDEX);
    }

    for(UInt idx=0; idx<cs.size(); idx++)
    {
      const UInt n = degrees[idx];
      const UInt m = orders[idx];
      if(cs[idx] == 0)
        idxC.at(n).at(m) = idx;
      else
        idxS.at(n).at(m) = idx;
    }

    matrix.resize(2*maxDegree+1);
    matrix.at(0) = Matrix(maxDegree+1, maxDegree+1);
    for(UInt m=1; m<=maxDegree; m++)
      matrix.at(2*m-1) = matrix.at(2*m-0) = Matrix(maxDegree+1-m, maxDegree+1-m);

    // order 0
    for(UInt n1=0; n1<=maxDegree; n1++) // rows
      for(UInt n2=0; n2<=maxDegree; n2++) // columns
        if((idxC[n1][0]!=NULLINDEX)&&(idxC[n2][0]!=NULLINDEX))
          matrix.at(0)(n1,n2) = N(idxC[n1][0], idxC[n2][0]);

    // order >0
    for(UInt m=1; m<=maxDegree; m++)
    {
      UInt idx1 = 0;
      for(UInt n1=m; n1<=maxDegree; n1++) // rows
      {
          UInt idx2 = 0;
          for(UInt n2=m; n2<=maxDegree; n2++) // columns
          {
            if((idxC[n1][m]!=NULLINDEX)&&(idxC[n2][m]!=NULLINDEX))
              matrix.at(2*m-1)(idx1, idx2) = N(idxC[n1][m], idxC[n2][m]);
            if((idxS[n1][m]!=NULLINDEX)&&(idxS[n2][m]!=NULLINDEX))
              matrix.at(2*m-0)(idx1, idx2) = N(idxS[n1][m], idxS[n2][m]);
            idx2++;
          } // for(n2)
        idx1++;
      } // for(n1)
    } // for(m)

    // order 0
    Matrix matrix2 = matrix.at(0);
    matrix2(0,0) = 1.0;
    for(UInt n=1; n<=maxDegree; n++)
      matrix2(n,n) += factor*pow(static_cast<Double>(n),power);
    matrix2.setType(Matrix::SYMMETRIC);
    solveInPlace(matrix2,matrix.at(0));

    // order >0
    for(UInt m=1; m<=maxDegree; m++)
    {
      Matrix matrix2c = matrix.at(2*m-1);
      Matrix matrix2s = matrix.at(2*m-0);
      UInt idx = 0;
      for(UInt n=m; n<=maxDegree; n++)
      {
        matrix2c(idx,idx) += factor*pow(static_cast<Double>(n),power);
        matrix2s(idx,idx) += factor*pow(static_cast<Double>(n),power);
        idx++;
      }
      matrix2c.setType(Matrix::SYMMETRIC);
      matrix2s.setType(Matrix::SYMMETRIC);
      solveInPlace(matrix2c,matrix.at(2*m-1));
      solveInPlace(matrix2s,matrix.at(2*m-0));
    }

    const UInt minDegreeFilter = 2;
    for(UInt n=0; n<minDegreeFilter; n++)
      matrix.at(0)(n,n) = 1.0;

    for(UInt m=1; m<=maxDegree; m++)
    {
      UInt idx = 0;
      for(UInt n=m; n<minDegreeFilter; n++)
      {
        matrix.at(2*m-1)(idx, idx) = 1.0;
        matrix.at(2*m-0)(idx, idx) = 1.0;
        idx++;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SphericalHarmonics SphericalHarmonicsFilterDdk::filter(const SphericalHarmonics &harm) const
{
  try
  {
    const UInt maxDegree = std::min(harm.maxDegree(), matrix.size()-1);

    Matrix cnm = harm.cnm();
    Matrix snm = harm.snm();

    Matrix cnmNew(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmNew(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    // filter order 0
    Vector x(maxDegree+1);
    for(UInt n=0; n<=maxDegree; n++)
      x(n) = cnm(n,0);
    x = matrix.at(0).slice(0,0,x.rows(),x.rows()) * x;
    for(UInt n=0; n<=maxDegree; n++)
      cnmNew(n,0) = x(n);

    // filter order >0
    for(UInt m=1; m<=maxDegree; m++)
    {
      Vector x(maxDegree+1-m);

      // cnm
      UInt idx = 0;
      for(UInt n=m; n<=maxDegree; n++)
        x(idx++) = cnm(n,m);
      x = matrix.at(2*m-1).slice(0,0,x.rows(),x.rows()) * x;
      idx = 0;
      for(UInt n=m; n<=maxDegree; n++)
        cnmNew(n,m) = x(idx++);

      // snm
      idx = 0;
      for(UInt n=m; n<=maxDegree; n++)
        x(idx++) = snm(n,m);
      x = matrix.at(2*m-0).slice(0,0,x.rows(),x.rows()) * x;
      idx = 0;
      for(UInt n=m; n<=maxDegree; n++)
        snmNew(n,m) = x(idx++);
    }

    return SphericalHarmonics(harm.GM(), harm.R(), cnmNew, snmNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
