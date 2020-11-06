/***********************************************/
/**
* @file gnssParametrizationIonosphere.cpp
*
* @brief GNSS ionosphere representation.
*
* @author Torsten Mayer-Guerr
* @date 2018-11-20
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationIonosphere

#include "base/import.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "inputOutput/logging.h"
#include "config/config.h"
#include "config/configRegister.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationIonosphere.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationIonosphere, "gnssParametrizationIonosphereType")
GROOPS_READCONFIG_CLASS(GnssParametrizationIonosphere, "gnssParametrizationIonosphereType")

/***********************************************/

GnssParametrizationIonosphere::GnssParametrizationIonosphere(Config &config, const std::string &name)
{
  try
  {
    sigmaSTEC = 0;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "outputfileVTEC",          fileNameVTEC,      Config::OPTIONAL,  "output/vtec_{loopTime:%D}.{station}.dat", "");
    readConfig(config, "magnetosphere",           magnetosphere,     Config::MUSTSET,   "",        "");
    readConfig(config, "apply2ndOrderCorrection", apply2ndOrder,     Config::DEFAULT,   "1",       "apply ionospheric correction");
    readConfig(config, "apply3rdOrderCorrection", apply3rdOrder,     Config::DEFAULT,   "1",       "apply ionospheric correction");
    readConfig(config, "applyBendingCorrection",  applyBending,      Config::DEFAULT,   "1",       "apply ionospheric correction");
//     readConfig(config, "sigmaSTEC",               sigmaSTEC,         Config::DEFAULT,   "50",      "(0 = unconstrained) sigma [TECU] for STEC constraint");
    readConfig(config, "singleLayerHeight",       singleLayerHeight, Config::DEFAULT,   "450e3",   "[m] height of the ionospheric pierce point");
    readConfig(config, "mapR",                    mapR,              Config::DEFAULT,   "6371e3",  "constant of MSLM mapping function");
    readConfig(config, "mapH",                    mapH,              Config::DEFAULT,   "506.7e3", "constant of MSLM mapping function");
    readConfig(config, "mapAlpha",                mapAlpha,          Config::DEFAULT,   "0.9782",  "constant of MSLM mapping function");
    endSequence(config);
    if(isCreateSchema(config)) return;

    apply1stOrder = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::initIntervalEarly(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    VTEC = Matrix();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// intersection point in ionosphere height
// classic approach is to use ~450/500 km, not applicable for satellites flying possibly higher
// Hence, satellite position + 50 km, no evidence if true, but LEOs are in principle flying
// IN the ionosphere and so pierce point is set slightly higher
Vector3d GnssParametrizationIonosphere::intersection(const Vector3d &posRecv, const Vector3d &posTrans) const
{
  try
  {
    const Vector3d k  = normalize(posRecv - posTrans); // line of sight from transmitter to receiver
    const Double   r  = posRecv.r();
    const Double   H  = std::max(DEFAULT_R + singleLayerHeight, r + 50e3);  // single layer ionosphere in 350 km or 50 km above satellite
    const Double   rk = inner(k, posRecv);
    return posRecv - (rk + std::sqrt(rk*rk + 2*r*H + H*H)) * k;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Mapping function
Double GnssParametrizationIonosphere::mapping(const Vector3d &/*posRecv*/, Angle elevation) const
{
  return 1./std::cos(std::asin(mapR/(mapR+mapH) * std::sin(mapAlpha*(PI/2-elevation))));
}

/***********************************************/
/***********************************************/

void GnssParametrizationIonosphere::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameterVTEC.clear();
    indexParameterVTEC.resize(gnss().receiver.size(), std::vector<Gnss::ParameterIndex>(gnss().times.size()));

    if(!(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_VTEC))
      return;

    if((VTEC.rows() != gnss().receiver.size()) || (VTEC.columns() != gnss().times.size()))
      VTEC = Matrix(gnss().receiver.size(), gnss().times.size());

    // distribute usable receiver epochs
    Matrix useableEpoch(gnss().receiver.size(), gnss().times.size());
    for(UInt idRecv=0; idRecv<gnss().receiver.size(); idRecv++)
      if(normalEquationInfo.estimateReceiver.at(gnss().receiver.at(idRecv)->idRecv()) && gnss().receiver.at(idRecv)->useable())
        for(UInt idEpoch=0; idEpoch<gnss().times.size(); idEpoch++)
          useableEpoch(idRecv, idEpoch) = gnss().receiver.at(idRecv)->isEpochEstimable(normalEquationInfo.analysisType, idEpoch);
    Parallel::reduceSum(useableEpoch, 0, normalEquationInfo.comm);
    Parallel::broadCast(useableEpoch, 0, normalEquationInfo.comm);

    // VTEC
    // ----
    for(UInt idRecv=0; idRecv<gnss().receiver.size(); idRecv++)
      if(sum(useableEpoch.row(idRecv)))
      {
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(useableEpoch(idRecv, idEpoch))
          {
            const ParameterName parameterName(gnss().receiver.at(idRecv)->name(), "VTEC", "", gnss().times.at(idEpoch));
            indexParameterVTEC.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, gnss().receiver.at(idRecv)->idRecv(), {parameterName});
          } // for(idRecv, idEpoch)
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::aprioriParameter(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, MatrixSliceRef /*x0*/) const
{
  try
  {
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::eliminateTecParameter(const Gnss::NormalEquationInfo &normalEquationInfo, Gnss::ObservationEquation &eqn) const
{
  try
  {
    if(!eqn.B.size())
      return;

    if(sigmaSTEC && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC))
    {
      // extend l, A, and B by one row at beginning
      for(Matrix &A : std::vector<std::reference_wrapper<Matrix>>{eqn.l, eqn.A, eqn.B})
      {
        const Matrix A0 = A;
        A = Matrix(1+A0.rows(), A0.columns());
        copy(A0, A.row(1, A0.rows()));
      }

      // constrain STEC;
      Matrix B = eqn.B;
      eqn.l(0) = -1./sigmaSTEC * eqn.dSTEC; // constrain towards zero (0-x0)
      B(0, 0)  =  1./sigmaSTEC;             // in TECU

      // eliminate STEC;
      eliminationParameter(B, {eqn.l, eqn.A, eqn.B});
    }
    else if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_VTEC)
    {
      eqn.l += eqn.B * eqn.dSTEC; // add back
    }
    else
    {
      // eliminate ionosphere parameter
      eliminationParameter(eqn.B, {eqn.A, eqn.l});
      eqn.B = Matrix();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::updateAndEliminateTecParameter(const Gnss::NormalEquationInfo &normalEquationInfo, Gnss::ObservationEquation &eqn,
                                                    Vector &We, Matrix &AWz, Double &maxChangeTec, std::string &infoMaxChangeTec) const
{
  try
  {
    if(!eqn.B.size())
      return;

    Matrix B = eqn.B;
    eqn.B = Matrix();
    if(sigmaSTEC && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC))
    {
      // extend We, AWz, and B by one row at beginning
      for(Matrix &A : std::vector<std::reference_wrapper<Matrix>>{We, AWz, B})
      {
        const Matrix A0 = A;
        A = Matrix(1+A0.rows(), A0.columns());
        copy(A0, A.row(1, A0.rows()));
      }

      // constrain STEC;
      We(0)   = -1./sigmaSTEC * eqn.dSTEC; // constrain towards zero (0-x0)
      B(0, 0) =  1./sigmaSTEC;             // in TECU
    }

    // estimate and eliminate STEC parameter
    Vector tau = QR_decomposition(B);
    QTransMult(B, tau, We);
    triangularSolve(1., B.row(0,B.columns()), We.row(0,B.columns()));
    eqn.receiver->observation(eqn.transmitter->idTrans(), eqn.idEpoch)->updateParameter(We.row(0,B.columns())); // ionosphere parameter

    if((std::fabs(We(0)) > maxChangeTec) && (norm(eqn.sigma-eqn.sigma0) < 1e-8)) // without outlier
    {
      maxChangeTec     = std::fabs(We(0));
      infoMaxChangeTec = "  "+eqn.receiver->name()+"."+eqn.transmitter->name()+":  "+eqn.timeRecv.dateTimeStr()+" tecChange = "
                        + We(0)%"%6.2f TECU = "s+(We(0)*Ionosphere::Ap/pow((GnssType::L2+GnssType::GPS).frequency(),2)*1000)%"%7.1f mm = "s
                        + (We(0)/LIGHT_VELOCITY*1e9)%"%8.3f ns"s;
    }

    // remove STEC from residuals
    We.row(0,B.columns()).setNull();
    QMult(B, tau, We);

    // influence of B parameters = B(B'B)^(-1)B'
    QTransMult(B, tau, AWz);
    AWz.row(0, B.columns()).setNull();
    for(UInt i=0; i<B.columns(); i++)
      AWz(i,i) = 1.0;
    QMult(B, tau, AWz);

    if(sigmaSTEC && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC))
    {
      // remove first row
      We  = We.row (1, We.rows() -1);
      AWz = AWz.row(1, AWz.rows()-1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationIonosphere::isDesignMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt /*idTrans*/, UInt idEpoch) const
{
  try
  {
    if(indexParameterVTEC.at(idRecv).at(idEpoch))
      return TRUE;

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::designMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    const UInt idRecv = eqn.receiver->idRecv();

    // VTEC at station per epoch
    // -------------------------
    if(indexParameterVTEC.at(idRecv).at(eqn.idEpoch))
      axpy(mapping(eqn.posRecv, eqn.elevationRecvLocal), eqn.B.column(0), A.column(indexParameterVTEC.at(idRecv).at(eqn.idEpoch)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIonosphere::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool printStatistics)
{
  try
  {
    Double minVTEC   = 1e+99;
    Double maxVTEC   = 1e-99;
    Double meanVTEC  = 0;
    Double stdVTEC   = 0;
    UInt   countVTEC = 0;

    Double      maxChangeVTEC = 0;
    std::string infoMaxChangeVTEC;
    for(UInt idRecv=0; idRecv<indexParameterVTEC.size(); idRecv++)
      for(UInt idEpoch=0; idEpoch<indexParameterVTEC.at(idRecv).size(); idEpoch++)
        if(indexParameterVTEC.at(idRecv).at(idEpoch))
        {
          const Double dVTEC = x(normalEquationInfo.index(indexParameterVTEC.at(idRecv).at(idEpoch)), 0);
          VTEC(idRecv, idEpoch) += dVTEC;
          if(std::fabs(dVTEC) > maxChangeVTEC)
          {
            maxChangeVTEC = std::fabs(dVTEC);
            infoMaxChangeVTEC = "  "+gnss().receiver.at(idRecv)->name()+": "+gnss().times.at(idEpoch).dateTimeStr()+" vtecChange     = "+dVTEC%"%6.2f TECU"s;
          }

          // VTEC statistics
          minVTEC   = std::min(VTEC(idRecv, idEpoch), minVTEC);
          maxVTEC   = std::max(VTEC(idRecv, idEpoch), maxVTEC);
          meanVTEC += VTEC(idRecv, idEpoch);
          stdVTEC  += VTEC(idRecv, idEpoch)*VTEC(idRecv, idEpoch);
          countVTEC++;
        }

    Gnss::checkMaxChange(maxChangeVTEC, infoMaxChangeVTEC, printStatistics, normalEquationInfo.comm);

    // VTEC analysis
    // -------------
    if(countVTEC)
    {
      stdVTEC   = sqrt((stdVTEC-meanVTEC*meanVTEC/countVTEC)/(countVTEC-1));
      meanVTEC /= countVTEC;
      logInfo<<"  VTEC: mean = "<<meanVTEC%"%.2f +- "s<<stdVTEC%"%.2f (range "s<<minVTEC%"%.2f -- "s<<maxVTEC%"%.2f) TECU"s<<Log::endl;
    }

    return 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphere::writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    if(!fileNameVTEC.empty() && VTEC.size())
    {
      VariableList fileNameVariableList;
      addTimeVariables(fileNameVariableList);
      addVariable("station", "****", fileNameVariableList);
      evaluateTimeVariables(0, gnss().times.at(0), gnss().times.back(), fileNameVariableList);

      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write VTEC to files <"<<fileNameVTEC(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      for(auto &recv : gnss().receiver)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable())
        {
          fileNameVariableList["station"]->setValue(recv->name());

          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(recv->useable(idEpoch) && VTEC(recv->idRecv(), idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = gnss().times.at(idEpoch);
              epoch.value = VTEC(recv->idRecv(), idEpoch);
              arc.push_back(epoch);
            }
          if(arc.size())
            InstrumentFile::write(fileNameVTEC(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Vector GnssParametrizationIonosphere::slantDelay(const Gnss::ObservationEquation &eqn, Matrix &B) const
{
  try
  {
    Double STEC = eqn.dSTEC;
    if(VTEC.size())
      STEC += mapping(eqn.posRecv, eqn.elevationRecvLocal) * VTEC(eqn.receiver->idRecv(), eqn.idEpoch);

    // ionospheric correction, notation see:
    // Fritsche, M., R. Dietrich, C. Knoefel, A. Ruelke, S. Vey, M. Rothacher, and P. Steigenberger (2005),
    // Impact of higher-order ionospheric terms on GPS estimates,
    // Geophys. Res. Lett., 32, L23311, doi:10.1029/2005GL024342.
    // ----------------------------------------------------------
    // second order magentic effect
    const Vector3d piercePoint = intersection(eqn.posRecv, eqn.posTrans);
    const Rotary3d rotEarth    = Planets::celestial2TerrestrialFrame(eqn.timeRecv);
    const Vector3d b           = rotEarth.inverseRotate(magnetosphere->magenticFieldVector(eqn.timeRecv, rotEarth.rotate(piercePoint))); // magentic field vector in CRF
    const Vector3d k           = normalize(eqn.posRecv - eqn.posTrans); // line of sight
    const Double   s           = 1e16*7527.*LIGHT_VELOCITY*inner(b, k);
    // third order
    constexpr Double r = 1e16*(2437 * 0.66 * (20.-6.)/(4.55-1.38)*1e-6) * 1e16;
    // bending
    constexpr Double HF2     = 70;                                                                           // F2 layer scale height
    const     Double hmF2    = std::max(350e3, eqn.posRecv.r() - DEFAULT_R + 50e3)/1000;                     // single layer ionosphere in 350 km or 50 km above satellite
    const     Double bending =  7.5*1e31*std::exp(-2.13*eqn.elevationRecvLocal)/(HF2*std::pow(hmF2, 0.125))  // Bending, Hoque and Jakowski 2008, approx formula
                             + -Ionosphere::Ap*0.1108*1e16*std::exp(-2.1844*eqn.elevationRecvLocal)/(HF2*std::pow(hmF2, 0.3));  // + TEC difference bending effect, Hoque and Jakowski 2008

    Vector l(eqn.types.size());
    B = Vector(eqn.types.size());
    for(UInt i=0; i<eqn.types.size(); i++)
    {
      const Double f1 = eqn.types.at(i).frequency();
      const Double f2 = f1*f1;
      const Double f3 = f1*f2;
      const Double f4 = f1*f3;

      // phase observations
      if(eqn.types.at(i) == GnssType::PHASE)
      {
        if(apply1stOrder) l(i)   += -Ionosphere::Ap/f2 * STEC;     // 1. order ionosphere correction (in terms of TEC)
        if(apply1stOrder) B(i,0) += -Ionosphere::Ap/f2;
        if(apply2ndOrder) l(i)   += -s/2/f3 * STEC;                // 2. order ionosphere correction (in terms of TEC)
        if(apply2ndOrder) B(i,0) += -s/2/f3;
        if(apply3rdOrder) l(i)   += -r/3/f4 * std::pow(STEC, 2);   // 3. order ionosphere correction (in terms of TEC)
        if(apply3rdOrder) B(i,0) += -r/3/f4 * 2*STEC;
        if(applyBending)  l(i)   +=  bending/f4*std::pow(STEC, 2); // Bending, Hoque and Jakowski 2008, approx formula
        if(applyBending)  B(i,0) +=  bending/f4 * 2*STEC;
      }

      // range observations
      if(eqn.types.at(i) == GnssType::RANGE)
      {
        if(apply1stOrder) l(i)   += Ionosphere::Ap/f2 * STEC; // 1. order ionosphere correction (in terms of TEC)
        if(apply1stOrder) B(i,0) += Ionosphere::Ap/f2;
        if(apply2ndOrder) l(i)   += s/f3  * STEC;             // 2. order ionosphere correction (in terms of TEC)
        if(apply2ndOrder) B(i,0) += s/f3;
        if(apply3rdOrder) l(i)   += r/f4 * std::pow(STEC, 2); // 3. order ionosphere correction (in terms of TEC)
        if(apply3rdOrder) B(i,0) += r/f4 * 2*STEC;
      }
    }

    return l;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
