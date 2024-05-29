/***********************************************/
/**
* @file gnssParametrizationIonosphereSTEC.cpp
*
* @brief IonosphereSTEC.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereSTEC.h"

/***********************************************/

GnssParametrizationIonosphereSTEC::GnssParametrizationIonosphereSTEC(Config &config)
{
  try
  {
    readConfig(config, "name",                    name,              Config::OPTIONAL, "parameter.STEC", "used for parameter selection");
    readConfig(config, "apply2ndOrderCorrection", apply2ndOrder,     Config::DEFAULT,  "1", "apply ionospheric correction");
    readConfig(config, "apply3rdOrderCorrection", apply3rdOrder,     Config::DEFAULT,  "1", "apply ionospheric correction");
    readConfig(config, "applyBendingCorrection",  applyBending,      Config::DEFAULT,  "1", "apply ionospheric correction");
    readConfig(config, "magnetosphere",           magnetosphere,     Config::MUSTSET,  "",  "");
    readConfig(config, "nameConstraint",          nameConstraint,    Config::OPTIONAL, "constraint.STEC", "used for parameter selection");
    readConfig(config, "sigmaSTEC",               exprSigmaSTEC,     Config::DEFAULT,  "0",  "expr. for sigma [TECU] for STEC constraint, variable E (elevation) available");
    if(isCreateSchema(config)) return;

    apply1stOrder = TRUE;

    // can exprSigmaSTEC be evaluated to a constant?
    // ---------------------------------------------
    VariableList varList;
    varList.undefineVariable("E");
    exprSigmaSTEC->simplify(varList);
    try
    {
      sigmaSTEC = (exprSigmaSTEC->evaluate(varList) > 0);
    }
    catch(std::exception &/*e*/)
    {
      sigmaSTEC = -1;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereSTEC::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;
    estimateSTEC    = TRUE;
    applyConstraint = FALSE;

    if(sigmaSTEC > 0) // sigmaSTEC is constant -> fast version
    {
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
        for(auto recv : gnss->receivers)
          if(recv->isMyRank() && recv->useable(idEpoch))
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
                recv->observation(idTrans, idEpoch)->sigmaSTEC = sigmaSTEC;
    }
    else if(sigmaSTEC) // sigmaSTEC is expression
    {
      logStatus<<"compute STEC accuracy"<<Log::endl;
      VariableList varList;
      Log::Timer timer(gnss->times.size());
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
      {
        timer.loopStep(idEpoch);
        for(auto recv : gnss->receivers)
          if(recv->isMyRank() && recv->useable(idEpoch))
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                            nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);
                varList.setVariable("E", eqn.elevationRecvLocal);
                recv->observation(idTrans, idEpoch)->sigmaSTEC = exprSigmaSTEC->evaluate(varList);
              } // for(recv, trans)
      } // for(idEpoch)
      Parallel::barrier(comm);
      timer.loopEnd();
    } // if(isSigmaSTEC)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereSTEC::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    estimateSTEC    = isEnabled(normalEquationInfo, name);
    applyConstraint = FALSE;
    if(!estimateSTEC)
      return;
    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaSTEC;

    UInt count = 0;
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
            if(recv->observation(idTrans, idEpoch))
              count++;
    Parallel::reduceSum(count, 0, normalEquationInfo.comm);
    logInfo<<count%"%9i "s<<(applyConstraint ? "constrained " : "")<<"STEC parameters (preeliminated)"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereSTEC::observationCorrections(GnssObservationEquation &eqn) const
{
  try
  {
    const Double STEC = eqn.STEC;

    // ionospheric correction, notation see:
    // Fritsche, M., R. Dietrich, C. Knoefel, A. Ruelke, S. Vey, M. Rothacher, and P. Steigenberger (2005),
    // Impact of higher-order ionospheric terms on GPS estimates,
    // Geophys. Res. Lett., 32, L23311, doi:10.1029/2005GL024342.
    // ----------------------------------------------------------
    // second order magentic effect
    constexpr Double   radiusIono  = 6371e3+450e3;                          // single layer ionosphere in 450 km
    const     Double   rRecv       = std::min(eqn.posRecv.r(), radiusIono); // LEO satellites flying possibly higher
    const     Vector3d k           = normalize(eqn.posRecv-eqn.posTrans);   // direction from transmitter
    const     Double   rk          = inner(eqn.posRecv, k);
    const     Vector3d piercePoint = eqn.posRecv - (std::sqrt(rk*rk+radiusIono*radiusIono-rRecv*rRecv)+rk) * k;
    const     Rotary3d rotEarth    = Planets::celestial2TerrestrialFrame(eqn.timeRecv);
    const     Vector3d b           = rotEarth.inverseRotate(magnetosphere->magenticFieldVector(eqn.timeRecv, rotEarth.rotate(piercePoint))); // magentic field vector in CRF
    const     Double   s           = 1e16*7527.*LIGHT_VELOCITY*inner(b, k);
    // third order
    constexpr Double r = 1e16*(2437 * 0.66 * (20.-6.)/(4.55-1.38)*1e-6) * 1e16;
    // bending
    constexpr Double HF2     = 70;                                                                           // F2 layer scale height
    const     Double hmF2    = std::max(350e3, eqn.posRecv.r() - DEFAULT_R + 50e3)/1000;                     // single layer ionosphere in 350 km or 50 km above satellite
    const     Double bending =  7.5*1e31*std::exp(-2.13*eqn.elevationRecvLocal)/(HF2*std::pow(hmF2, 0.125))  // Bending, Hoque and Jakowski 2008, approx formula
                             + -Ionosphere::Ap*0.1108*1e16*std::exp(-2.1844*eqn.elevationRecvLocal)/(HF2*std::pow(hmF2, 0.3));  // + TEC difference bending effect, Hoque and Jakowski 2008

    eqn.B = Vector(eqn.types.size());
    for(UInt i=0; i<eqn.types.size(); i++)
    {
      const Double f1 = eqn.types.at(i).frequency();
      const Double f2 = f1*f1;
      const Double f3 = f1*f2;
      const Double f4 = f1*f3;

      // phase observations
      if(eqn.types.at(i) == GnssType::PHASE)
      {
        if(apply1stOrder) eqn.l(i)   -= -Ionosphere::Ap/f2 * STEC;     // 1. order ionosphere correction (in terms of TEC)
        if(apply1stOrder) eqn.B(i,0) += -Ionosphere::Ap/f2;
        if(apply2ndOrder) eqn.l(i)   -= -s/2/f3 * STEC;                // 2. order ionosphere correction (in terms of TEC)
        if(apply2ndOrder) eqn.B(i,0) += -s/2/f3;
        if(apply3rdOrder) eqn.l(i)   -= -r/3/f4 * std::pow(STEC, 2);   // 3. order ionosphere correction (in terms of TEC)
        if(apply3rdOrder) eqn.B(i,0) += -r/3/f4 * 2*STEC;
        if(applyBending)  eqn.l(i)   -=  bending/f4*std::pow(STEC, 2); // Bending, Hoque and Jakowski 2008, approx formula
        if(applyBending)  eqn.B(i,0) +=  bending/f4 * 2*STEC;
      }

      // range observations
      if(eqn.types.at(i) == GnssType::RANGE)
      {
        if(apply1stOrder) eqn.l(i)   -= Ionosphere::Ap/f2 * STEC; // 1. order ionosphere correction (in terms of TEC)
        if(apply1stOrder) eqn.B(i,0) += Ionosphere::Ap/f2;
        if(apply2ndOrder) eqn.l(i)   -= s/f3  * STEC;             // 2. order ionosphere correction (in terms of TEC)
        if(apply2ndOrder) eqn.B(i,0) += s/f3;
        if(apply3rdOrder) eqn.l(i)   -= r/f4 * std::pow(STEC, 2); // 3. order ionosphere correction (in terms of TEC)
        if(apply3rdOrder) eqn.B(i,0) += r/f4 * 2*STEC;
      }
    }
    copy(eqn.B, eqn.A.column(GnssObservationEquation::idxSTEC));
    if(!estimateSTEC)
      eqn.B = Vector();

    // add additional constraining equation
    // ------------------------------------
    if(applyConstraint && (eqn.sigmaSTEC > 0) && !std::isnan(eqn.sigmaSTEC))
    {
      // extend l, A, and B by one row
      for(Matrix &A : std::vector<std::reference_wrapper<Matrix>>{eqn.l, eqn.A, eqn.B})
      {
        const Matrix A0 = A;
        A = Matrix(A0.rows()+1, A0.columns());
        copy(A0, A.row(0, A0.rows()));
      }

      // constrain STEC;
      eqn.l(eqn.l.rows()-1)    = -eqn.dSTEC/eqn.sigmaSTEC; // constrain towards zero (0-f(x0))
      eqn.B(eqn.B.rows()-1, 0) =  1./eqn.sigmaSTEC;        // in TECU
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
