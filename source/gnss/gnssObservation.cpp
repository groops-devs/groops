/***********************************************/
/**
* @file gnssObservation.cpp
*
* @brief Code & Phase observations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-18
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrizationIonosphere.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssObservation.h"

/***********************************************/

Bool Gnss::Observation::init(const Receiver &receiver, const Transmitter &transmitter, GnssParametrizationIonospherePtr /*ionosphere*/, UInt idEpoch, Double &phaseWindupOld)
{
  try
  {
    if((!receiver.useable(idEpoch)) || (!transmitter.useable(idEpoch)))
      return FALSE;

    // position, time of transmitter & receiver
    // ----------------------------------------
    const Time     timeRecv = receiver.timeCorrected(idEpoch);
    const Vector3d posRecv  = receiver.position(idEpoch);
    Time     timeTrans;
    Vector3d posTrans;
    transmitter.transmitTime(idEpoch, timeRecv, posRecv, timeTrans, posTrans);
    const Vector3d k              = normalize(posRecv - posTrans);                    // line of sight from transmitter to receiver
    const Vector3d kRecv          = receiver.local2antennaFrame(idEpoch).transform(receiver.celestial2localFrame(idEpoch).transform(-k)); // line of sight in receiver antenna system (north, east, up)
    const Angle    azimutRecv     = kRecv.lambda();
    const Angle    elevationRecv  = kRecv.phi();
    const Vector3d kTrans         = transmitter.celestial2antennaFrame(idEpoch, timeTrans).transform(k);
    const Angle    azimutTrans    = kTrans.lambda();
    const Angle    elevationTrans = kTrans.phi();

    if(elevationRecv<receiver.elevationCutOff())
      return FALSE;

    std::vector<GnssType> types(size());
    for(UInt i=0; i<size(); i++)
      types.at(i) = at(i).type;

    // Composed signals (e.g. C2DG)
    std::vector<GnssType> typesTransmitted;
    Matrix T;
    receiver.signalComposition(idEpoch, types, typesTransmitted, T);
    if(!T.size())
      return FALSE;

    // Accuracies and antenna pattern check
    // ------------------------------------
    const Vector sigma0 = receiver.accuracy(idEpoch, azimutRecv, elevationRecv, types);
    const Vector acv    = receiver.antennaVariations(idEpoch, azimutRecv, elevationRecv, types)
                        + T * transmitter.antennaVariations(idEpoch, azimutTrans, elevationTrans, typesTransmitted);
    for(UInt i=0; i<size(); i++)
    {
      at(i).sigma0 = sigma0(i);
      at(i).sigma  = acv(i); // temporarily misuse sigma for ACV pattern nan check
    }
    auto iter = std::remove_if(obs.begin(), obs.end(), [](const auto &x)
    {
      if(((x.type == GnssType::PHASE) || (x.type == GnssType::RANGE)) && (x.sigma0 <= 0)) return TRUE; // remove Phase/Range for sigma <= 0
      return std::isnan(x.sigma0) || std::isnan(x.sigma);                                              // remove all NAN values
    });
    obs.erase(iter, obs.end());
    obs.shrink_to_fit();
    if(size()==0)
      return FALSE;

    for(UInt i=0; i<size(); i++)
      at(i).sigma = at(i).sigma0;
    for(UInt i=0; i<size(); i++)
      at(i).residuals = at(i).redundancy = 0.;

    std::sort(obs.begin(), obs.end(), [](const SingleObservation &obs1, const SingleObservation &obs2){return (obs1.type < obs2.type);});

    // phase wind-up
    // Carrier phase wind-up in GPS reflectometry, Georg Beyerle, Springer Verlag 2008
    // -------------------------------------------------------------------------------
    const Transform3d crf2arfRecv  = receiver.local2antennaFrame(idEpoch) * receiver.celestial2localFrame(idEpoch);
    const Transform3d crf2arfTrans = transmitter.celestial2antennaFrame(idEpoch, timeTrans);
    const Vector3d Tx = crf2arfRecv.transform(crossProduct(crossProduct(k, crf2arfTrans.inverseTransform(Vector3d(1,0,0))), k));
    const Vector3d Ty = crf2arfRecv.transform(crossProduct(crossProduct(k, crf2arfTrans.inverseTransform(Vector3d(0,1,0))), k));
    Double phaseWindup = atan2(Tx.y()+Ty.x(), Ty.y()-Tx.x()); // both left-handed systems
    phaseWindup -= PI/2; // to be consistent with Wu et al. (1993) and Kouba (2009) definition
    while((phaseWindupOld-phaseWindup)>PI)
      phaseWindup += 2*PI;
    while((phaseWindupOld-phaseWindup)<-PI)
      phaseWindup -= 2*PI;
    phaseWindupOld = phaseWindup;
    for(UInt i=0; i<size(); i++)
      if(at(i).type == GnssType::PHASE)
        at(i).observation -= phaseWindup/(2*PI) * (LIGHT_VELOCITY/at(i).type.frequency());

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt Gnss::Observation::index(GnssType type) const
{
  for(UInt i=0; i<size(); i++)
    if(at(i).type == type)
      return i;
  return NULLINDEX;
}

/***********************************************/

Bool Gnss::Observation::observationList(Gnss::AnalysisType analysisType, std::vector<GnssType> &types) const
{
  try
  {
    types.clear();

    // code observations
    if(analysisType & Gnss::ANALYSIS_CODE)
    {
      // need obs at two frequencies
      std::set<GnssType> typeFrequencies;
      for(UInt i=0; i<size(); i++)
        if(at(i).type == GnssType::RANGE)
          typeFrequencies.insert(at(i).type & GnssType::FREQUENCY);
      if(typeFrequencies.size() < 2)
        return FALSE;

      for(UInt i=0; i<size(); i++)
        if(at(i).sigma>0)
          if(at(i).type == GnssType::RANGE)
            types.push_back( at(i).type );
    }

    // phase observations
    if(analysisType & Gnss::ANALYSIS_PHASE)
    {
      // need obs at two frequencies
      std::set<GnssType> typeFrequencies;
      for(UInt i=0; i<size(); i++)
        if(at(i).type == GnssType::PHASE)
          typeFrequencies.insert(at(i).type & GnssType::FREQUENCY);
      if(typeFrequencies.size() < 2)
        return FALSE;

      for(UInt i=0; i<size(); i++)
        if(at(i).sigma>0)
          if(at(i).type == GnssType::PHASE)
            types.push_back( at(i).type );
    }

    // ionospheric delay as pseudo observations
    if(analysisType & Gnss::ANALYSIS_IONO)
    {
      for(UInt i=0; i<size(); i++)
        if(at(i).sigma>0)
          if(at(i).type == GnssType::IONODELAY)
            types.push_back( at(i).type );
    }

    // all available observations
//     if(analysisType & Gnss::ANALYSIS_RAW)
//     {
//       for(UInt i=0; i<size(); i++)
//         if(at(i).sigma>0)
//           if(GnnsType::index(type, at(i).type) == NULLINDEX)
//             type.push_back( at(i).type );
//     }

//     std::sort(type.begin(), type.end());
    return (types.size() != 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Observation::setDecorrelatedResiduals(const std::vector<GnssType> &types, const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy)
{
  try
  {
    for(UInt i=0; i<types.size(); i++)
    {
      const UInt idx = index(types.at(i));
      at(idx).residuals  = residuals(i,0) * at(idx).sigma;
      at(idx).redundancy = redundancy(i,0);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Observation::updateParameter(const_MatrixSliceRef x, const_MatrixSliceRef /*covariance*/)
{
  try
  {
    dSTEC += x(0,0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Gnss::ObservationEquation::compute(const Observation &observation, const Receiver &receiver_, const Transmitter &transmitter_,
                                        GnssParametrizationIonospherePtr ionosphere, UInt idEpoch_, Bool decorrelate, const std::vector<GnssType> &types_)
{
  try
  {
    const UInt obsCount = types_.size();
    types  = types_;
    l      = Vector(obsCount);
    sigma  = Vector(obsCount);
    sigma0 = Vector(obsCount);

    for(UInt i=0; i<obsCount; i++)
    {
      const UInt idx = observation.index(types.at(i));
      if(idx == NULLINDEX)
        continue;
      types.at(i) = observation.at(idx).type;
      l(i)        = observation.at(idx).observation;
      sigma(i)    = observation.at(idx).sigma;
      sigma0(i)   = observation.at(idx).sigma0;
    }

    // position, time of transmitter & receiver
    // ----------------------------------------
    receiver      = &receiver_;
    transmitter   = &transmitter_;
    idEpoch       = idEpoch_;
    timeRecv      = receiver->timeCorrected(idEpoch);
    posRecv       = receiver->position(idEpoch);
    transmitter->transmitTime(idEpoch, timeRecv, posRecv, timeTrans, posTrans);
    velocityRecv  = receiver->velocity(idEpoch);
    velocityTrans = transmitter->velocity(idEpoch, timeTrans);

    // orientation of antennas
    // -----------------------
    const Vector3d k          = normalize(posRecv - posTrans);                               // line of sight from transmitter to receiver
    const Double   rDotRecv   = inner(k, velocityRecv) /LIGHT_VELOCITY;
    const Vector3d kRecvLocal = receiver->celestial2localFrame(idEpoch).transform(-k);       // line of sight in local frame (north, east, up or vehicle)
    azimutRecvLocal           = kRecvLocal.lambda();
    elevationRecvLocal        = kRecvLocal.phi();
    const Vector3d kRecvAnt   = receiver->local2antennaFrame(idEpoch).transform(kRecvLocal); // line of sight in left-handed antenna system (useally north, east, up)
    azimutRecvAnt             = kRecvAnt.lambda();
    elevationRecvAnt          = kRecvAnt.phi();

    const Vector3d kTrans    = transmitter->celestial2antennaFrame(idEpoch, timeTrans).transform(k);
    const Double   rDotTrans = inner(k, velocityTrans)/LIGHT_VELOCITY;
    azimutTrans              = kTrans.lambda();
    elevationTrans           = kTrans.phi();

    // Corrected range
    // ---------------
    const Double r1 = posTrans.r();
    const Double r2 = posRecv.r();
    r12  = (posRecv - posTrans).r();
    r12 += 2*DEFAULT_GM/pow(LIGHT_VELOCITY,2)*log((r1+r2+r12)/(r1+r2-r12)); // curved space-time
    r12 += 2*inner(posTrans, velocityTrans)/LIGHT_VELOCITY;                 // relativistic clock correction
    r12 += receiver->troposphere(idEpoch, azimutRecvLocal, elevationRecvLocal);
    r12 -= LIGHT_VELOCITY * transmitter->clockError(idEpoch, timeTrans);
    r12 += LIGHT_VELOCITY * receiver->clockError(idEpoch);

    // approximate range
    // -----------------
    for(UInt i=0; i<obsCount; i++)
      if((types.at(i) == GnssType::RANGE) || (types.at(i) == GnssType::PHASE))
        l(i) -= r12;

    // reduce ambiguities
    // ------------------
    track = observation.track;
    if(track && track->ambiguity)
      l -= track->ambiguity->ambiguities(types);

    // Composed signals (e.g. C2DG)
    Matrix T;
    receiver->signalComposition(idEpoch, types, typesTransmitted, T);

    // antenna correction
    // ------------------
    l -= receiver->antennaVariations(idEpoch, azimutRecvAnt,  elevationRecvAnt,  types);
    l -= T * transmitter->antennaVariations(idEpoch, azimutTrans, elevationTrans, typesTransmitted);

    // reduce ionospheric effects
    // --------------------------
    dSTEC = observation.dSTEC;
    B = Matrix();
    if(ionosphere)
      l -= ionosphere->slantDelay(*this, B);

    // design matrix
    // -------------
    A = Matrix(obsCount, 9 + obsCount + T.columns());
    for(UInt i=0; i<obsCount; i++)
    {
      if((types.at(i) == GnssType::RANGE) || (types.at(i) == GnssType::PHASE))
      {
        A(i, idxPosRecv+0)  = k.x() * (1+rDotTrans);   // receiver coord x
        A(i, idxPosRecv+1)  = k.y() * (1+rDotTrans);   // receiver coord y
        A(i, idxPosRecv+2)  = k.z() * (1+rDotTrans);   // receiver coord z
        A(i, idxClockRecv)  = 1.0+rDotTrans-rDotRecv;  // receiver clock correction
        A(i, idxPosTrans+0) = k.x() * (-1+rDotTrans);  // transmitter coord x
        A(i, idxPosTrans+1) = k.y() * (-1+rDotTrans);  // transmitter coord y
        A(i, idxPosTrans+2) = k.z() * (-1+rDotTrans);  // transmitter coord z
        A(i, idxClockTrans) = -1.0;                    // transmitter clock correction
        A(i, idxRange)      = 1.0;                     // range correction (troposphere, ...)
      }

      A(i, idxUnit+i) = 1.0; // unit matrix
    }  // for(i=0..obsCount)
    copy(T, A.column(idxUnit + obsCount, T.columns()));

    // Decorrelate
    // -----------
    if(decorrelate)
      for(UInt i=0; i<obsCount; i++)
      {
        if(l.size()) l.row(i) *= 1/sigma(i);
        if(A.size()) A.row(i) *= 1/sigma(i);
        if(B.size()) B.row(i) *= 1/sigma(i);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
