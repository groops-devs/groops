/***********************************************/
/**
* @file slrParametrizationGroup.h
*
* @brief Group.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONGROUP__
#define __GROOPS_SLRPARAMETRIZATIONGROUP__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationGroup = R"(
\subsection{Group}\label{slrParametrizationType:group}
Groups a set of parameters. This class can be used to structure complex parametrizations
and has no further effect itself.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Group.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationGroup : public SlrParametrizationBase
{
  SlrParametrizationPtr base;

public:
  SlrParametrizationGroup(Config &config);
 ~SlrParametrizationGroup() {}

  void   getParametrizationGravity(std::vector<const SlrParametrizationGravityField*> &paramGravityField) const override;
  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   observationCorrections(SlrObservationEquation &eqn) const override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  void   constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   updateCovariance(const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

inline SlrParametrizationGroup::SlrParametrizationGroup(Config &config)
{
  try
  {
    readConfig(config, "parametrization", base, Config::MUSTSET, "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::getParametrizationGravity(std::vector<const SlrParametrizationGravityField*> &paramGravityField) const
{
  std::vector<const SlrParametrizationGravityField*> params = base->getParametrizationGravity();
  paramGravityField.insert(paramGravityField.end(), params.begin(), params.end());
}

/***********************************************/

inline void SlrParametrizationGroup::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField)
{
  try
  {
    base->init(slr, paramGravityField);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    base->initParameter(normalEquationInfo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    base->observationCorrections(eqn);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    x0 += base->aprioriParameter(normalEquationInfo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    base->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    base->constraints(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double SlrParametrizationGroup::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    return base->updateParameter(normalEquationInfo, x, Wz);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::updateCovariance(const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    base->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationGroup::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    base->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
