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
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssObservation.h"

/***********************************************/

static void positionVelocityTime(const GnssReceiver &receiver, const GnssTransmitter &transmitter, const Rotary3d &rotCrf2Trf, UInt idEpoch,
                                 Time &timeRecv,  Vector3d &posRecv,  Vector3d &velRecv,  Angle &azimutRecv,  Angle &elevationRecv,
                                 Time &timeTrans, Vector3d &posTrans, Vector3d &velTrans, Angle &azimutTrans, Angle &elevationTrans,
                                 Vector3d &k, Vector3d &kRecv, Vector3d &kTrans)
{
  try
  {
    // receiver in celestial reference frame
    timeRecv = receiver.timeCorrected(idEpoch);
    posRecv  = rotCrf2Trf.inverseRotate(receiver.position(idEpoch));
    velRecv  = rotCrf2Trf.inverseRotate(receiver.velocity(idEpoch));
    if(receiver.isEarthFixed())
      velRecv += crossProduct(Vector3d(0., 0., 7.29211585531e-5), posRecv);

    // transmitter position and time
    posTrans = transmitter.position(idEpoch, timeRecv-seconds2time(20200e3/LIGHT_VELOCITY));
    Vector3d posOld;
    for(UInt i=0; (i<10) && ((posTrans-posOld).r() > 0.0001); i++) // iteration
    {
      timeTrans = timeRecv - seconds2time((posTrans-posRecv).r()/LIGHT_VELOCITY);
      posOld    = posTrans;
      posTrans  = transmitter.position(idEpoch, timeTrans);
    }
    velTrans = transmitter.velocity(timeTrans);

    // line of sight from transmitter to receiver
    k              = normalize(posRecv - posTrans);
    kRecv          = receiver.global2antennaFrame(idEpoch).transform(rotCrf2Trf.rotate(-k)); // line of sight in receiver antenna system (north, east, up)
    azimutRecv     = kRecv.lambda();
    elevationRecv  = kRecv.phi();
    kTrans         = transmitter.celestial2antennaFrame(idEpoch, timeTrans).transform(k);
    azimutTrans    = kTrans.lambda();
    elevationTrans = kTrans.phi();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Bool GnssObservation::init(const GnssReceiver &receiver, const GnssTransmitter &transmitter, const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                           UInt idEpoch, Angle elevationCutOff, Double &phaseWindupOld)
{
  try
  {
    if((!receiver.useable(idEpoch)) || (!transmitter.useable(idEpoch)))
      return FALSE;

    // position, time of transmitter & receiver
    // ----------------------------------------
    Time     timeRecv, timeTrans;
    Vector3d posRecv,  velRecv, posTrans, velTrans, k, kRecv, kTrans;
    Angle    azimutRecv, elevationRecv, azimutTrans, elevationTrans;
    Rotary3d rotCrf2Trf;
    if(receiver.isEarthFixed())
      rotCrf2Trf = rotationCrf2Trf(receiver.timeCorrected(idEpoch));
    positionVelocityTime(receiver, transmitter, rotCrf2Trf, idEpoch,
                         timeRecv, posRecv, velRecv, azimutRecv, elevationRecv,
                         timeTrans, posTrans, velTrans, azimutTrans, elevationTrans,
                         k, kRecv, kTrans);

    if(elevationRecv < elevationCutOff)
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
    const Vector sigma0   = receiver.accuracy(timeRecv, azimutRecv, elevationRecv, types);
    const Vector acvRecv  = receiver.antennaVariations(timeRecv, azimutRecv, elevationRecv, types);
    const Vector acvTrans = transmitter.antennaVariations(timeTrans, azimutTrans, elevationTrans, typesTransmitted);
    for(UInt i=0; i<size(); i++)
    {
      at(i).sigma0 = sigma0(i);
      at(i).sigma  = acvRecv(i); // temporarily misuse sigma for ACV pattern nan check
      for(UInt k=0; k<T.columns(); k++)
        if(T(i,k))
          at(i).sigma  += T(i,k) * acvTrans(k);
    }
    obs.erase(std::remove_if(obs.begin(), obs.end(), [](const auto &x)
    {
      if(((x.type == GnssType::PHASE) || (x.type == GnssType::RANGE)) && (x.sigma0 <= 0)) return TRUE; // remove Phase/Range for sigma <= 0
      return std::isnan(x.sigma0) || std::isnan(x.sigma);                                              // remove all NAN values
    }), obs.end());
    obs.shrink_to_fit();
    if(size()==0)
      return FALSE;

    for(UInt i=0; i<size(); i++)
      at(i).sigma = at(i).sigma0;
    for(UInt i=0; i<size(); i++)
      at(i).residuals = at(i).redundancy = 0.;

    std::sort(obs.begin(), obs.end(), [](const GnssSingleObservation &obs1, const GnssSingleObservation &obs2) {return (obs1.type < obs2.type);});

    // phase wind-up
    // Carrier phase wind-up in GPS reflectometry, Georg Beyerle, Springer Verlag 2008
    // -------------------------------------------------------------------------------
    const Transform3d crf2arfRecv  = receiver.global2antennaFrame(idEpoch) * rotCrf2Trf;
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
        at(i).observation -= phaseWindup/(2*PI) * (at(i).type.wavelength());

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssObservation::index(GnssType type) const
{
  for(UInt i=0; i<size(); i++)
    if(at(i).type == type)
      return i;
  return NULLINDEX;
}

/***********************************************/

Bool GnssObservation::observationList(GnssObservation::Group group, std::vector<GnssType> &types) const
{
  try
  {
    types.clear();

    // code observations
    if(group & RANGE)
    {
      // need obs at two frequencies
      std::set<GnssType> typeFrequencies;
      for(UInt i=0; i<size(); i++)
        if((at(i).type == GnssType::RANGE) && (at(i).sigma > 0))
          typeFrequencies.insert(at(i).type & GnssType::FREQUENCY);
      if(typeFrequencies.size() < 2)
        return FALSE;

      for(UInt i=0; i<size(); i++)
        if(at(i).sigma>0)
          if(at(i).type == GnssType::RANGE)
            types.push_back( at(i).type );
    }

    // phase observations
    if(group & PHASE)
    {
      // need obs at two frequencies
      std::set<GnssType> typeFrequencies;
      for(UInt i=0; i<size(); i++)
        if((at(i).type == GnssType::PHASE) && (at(i).sigma > 0))
          typeFrequencies.insert(at(i).type & GnssType::FREQUENCY);
      if(typeFrequencies.size() < 2)
        return FALSE;

      for(UInt i=0; i<size(); i++)
        if(at(i).sigma>0)
          if(at(i).type == GnssType::PHASE)
            types.push_back( at(i).type );
    }

    // ionospheric delay as pseudo observations
    if(group & IONO)
    {
      for(UInt i=0; i<size(); i++)
        if((at(i).type == GnssType::IONODELAY) && (at(i).sigma > 0))
          types.push_back( at(i).type );
    }

    return (types.size() != 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssObservation::setHomogenizedResiduals(const std::vector<GnssType> &types, const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy)
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

void GnssObservation::updateParameter(const_MatrixSliceRef x, const_MatrixSliceRef /*covariance*/)
{
  try
  {
    STEC  += x(0,0);
    dSTEC += x(0,0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GnssObservationEquation::compute(const GnssObservation &observation, const GnssReceiver &receiver_, const GnssTransmitter &transmitter_,
                                      const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                      UInt idEpoch_, Bool homogenize, const std::vector<GnssType> &types_)
{
  try
  {
    const UInt obsCount = types_.size();
    idEpoch     = idEpoch_;
    receiver    = &receiver_;
    transmitter = &transmitter_;
    track       = observation.track;
    STEC        = observation.STEC;
    dSTEC       = observation.dSTEC;
    sigmaSTEC   = observation.sigmaSTEC;
    types       = types_;
    rankDeficit = 0;
    l           = Vector(obsCount);
    sigma       = Vector(obsCount);
    sigma0      = Vector(obsCount);

    // position, time of transmitter & receiver
    // ----------------------------------------
    Vector3d k, kRecvAnt, kTrans;
    Rotary3d rotCrf2Trf;
    if(receiver_.isEarthFixed())
      rotCrf2Trf = rotationCrf2Trf(receiver_.timeCorrected(idEpoch));
    positionVelocityTime(receiver_, transmitter_, rotCrf2Trf, idEpoch_,
                         timeRecv,  posRecv,  velocityRecv,  azimutRecvAnt, elevationRecvAnt,
                         timeTrans, posTrans, velocityTrans, azimutTrans,   elevationTrans,
                         k, kRecvAnt, kTrans);
    const Double   rDotTrans  = inner(k, velocityTrans)/LIGHT_VELOCITY;
    const Double   rDotRecv   = inner(k, velocityRecv) /LIGHT_VELOCITY;
    const Vector3d kRecvLocal = receiver_.global2localFrame(idEpoch).transform(rotCrf2Trf.rotate(-k));
    azimutRecvLocal    = kRecvLocal.lambda();
    elevationRecvLocal = kRecvLocal.phi();

    if(!obsCount)
      return;

    // observations
    // ------------
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

    // Corrected range
    // ---------------
    const Double r1  = posTrans.r();
    const Double r2  = posRecv.r();
    Double       r12 = (posRecv - posTrans).r();
    r12 += 2*DEFAULT_GM/pow(LIGHT_VELOCITY,2)*log((r1+r2+r12)/(r1+r2-r12)); // curved space-time
    r12 += 2*inner(posTrans, velocityTrans)/LIGHT_VELOCITY;                 // relativistic clock correction
    r12 -= LIGHT_VELOCITY * transmitter->clockError(idEpoch);
    r12 += LIGHT_VELOCITY * receiver->clockError(idEpoch);
    r12 *= transmitter->scaleFactor(idEpoch);

    for(UInt i=0; i<obsCount; i++)
      if((types.at(i) == GnssType::RANGE) || (types.at(i) == GnssType::PHASE))
        l(i) -= r12;

    // Composed signals (e.g. C2DG)
    Matrix T;
    receiver->signalComposition(idEpoch, types, typesTransmitted, T);

    // design matrix
    // -------------
    A = Matrix(obsCount, 10 + obsCount + T.columns());
    B = Matrix();
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
        A(i, idxRange)      =  1.0;                    // range correction (troposphere, ...)
        A(i, idxSTEC)       = types.at(i).ionosphericFactor();
      }

      A(i, idxUnit+i) = 1.0; // unit matrix
    }  // for(i=0..obsCount)
    copy(T, A.column(idxUnit + obsCount, T.columns()));

    // antenna correction and other corrections
    // ----------------------------------------
    l -= receiver->antennaVariations(timeRecv, azimutRecvAnt,  elevationRecvAnt,  types);
    l -= T * transmitter->antennaVariations(timeTrans, azimutTrans, elevationTrans, typesTransmitted);
    if(track && track->ambiguity)
      l -= track->ambiguity->ambiguities(types);     // reduce ambiguities
    if(reduceModels)
      reduceModels(*this);

    // Homogenize
    // ----------
    if(homogenize)
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

void GnssObservationEquation::eliminateGroupParameters(Bool removeRows)
{
  try
  {
    if(!B.size())
      return;
    const Vector tau = QR_decomposition(B);
    if(!removeRows)
      rankDeficit += B.columns();
    for(Matrix &A : std::vector<std::reference_wrapper<Matrix>>{A, l})
      if(A.size())
      {
        QTransMult(B, tau, A);
        if(removeRows)
          A = A.row(B.columns(), A.rows()-B.columns());
        else
        {
          A.row(0, B.columns()).setNull();
          QMult(B, tau, A);
        }
      }
    B = Matrix();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
