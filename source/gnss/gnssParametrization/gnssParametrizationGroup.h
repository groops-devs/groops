/***********************************************/
/**
* @file gnssParametrizationGroup.h
*
* @brief Group.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONGROUP__
#define __GROOPS_GNSSPARAMETRIZATIONGROUP__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationGroup = R"(
\subsection{Group}\label{gnssParametrizationType:group}
Groups a set of parameters. This class can be used to structure complex parametrizations
and has no further effect itself.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Group.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationGroup : public GnssParametrizationBase
{
  GnssParametrizationPtr base;

public:
  GnssParametrizationGroup(Config &config);
 ~GnssParametrizationGroup() {}

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                      std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   observationCorrections(GnssObservationEquation &eqn) const override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  void   constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                          const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger) override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

inline GnssParametrizationGroup::GnssParametrizationGroup(Config &config)
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

inline void GnssParametrizationGroup::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    base->init(gnss, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationGroup::requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                                                   std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch)
{
  try
  {
    base->requirements(normalEquationInfo, transCount, transCountEpoch, recvCount, recvCountEpoch);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationGroup::initParameter(GnssNormalEquationInfo &normalEquationInfo)
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

inline void GnssParametrizationGroup::observationCorrections(GnssObservationEquation &eqn) const
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

inline void GnssParametrizationGroup::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
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

inline void GnssParametrizationGroup::designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
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

inline void GnssParametrizationGroup::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    base->constraintsEpoch(normalEquationInfo, idEpoch, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationGroup::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
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

inline Double GnssParametrizationGroup::ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                                                         const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger)
{
  try
  {
    return base->ambiguityResolve(normalEquationInfo, normals, n, lPl, obsCount, searchInteger);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double GnssParametrizationGroup::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
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

inline void GnssParametrizationGroup::updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
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

inline void GnssParametrizationGroup::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
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
