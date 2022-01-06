/***********************************************/
/**
* @file integralEquation.cpp
*
* @brief Functions for the (short arc) integral equation approach.
*
* @author Torsten Mayer-Guerr
* @date 2005-02-24
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "integralEquation.h"

/***********************************************/

IntegralEquation::IntegralEquation(UInt integrationDegree, UInt interpolationDegree)
{
  init(integrationDegree, interpolationDegree);
}

/***********************************************/

void IntegralEquation::init(UInt integrationDegree, UInt interpolationDegree)
{
  try
  {
    this->interpolationDegree = interpolationDegree;
    this->integrationDegree   = integrationDegree;
    if(integrationDegree%2 == 0)
      throw(Exception("polnomial degree for integration must be odd."));

    // Es gibt grad Gleichungssysteme
    // grad/2 am Anfang des Bahnbogens (es fehlen Punkte zum integrieren)
    // grad/2 am Ende des Bahnbogens (es fehlen Punkte zum integrieren)
    // 1 fuer den Rest des Bogens
    W.resize(integrationDegree);

    UInt half = (integrationDegree-1)/2;
    for(UInt k=0; k<half; k++)
    {
      W.at(k) = Matrix(integrationDegree+1, integrationDegree+1);
      for(UInt i=0; i<W.at(k).columns(); i++)
        for(UInt n=0; n<W.at(k).rows(); n++)
          W.at(k)(n,i) = ((n==0) ? 1.0 : pow(static_cast<Double>(i)-k, n));
      inverse(W.at(k));
    }

    for(UInt k=0; k<=half; k++)
    {
      W.at(k+half) = Matrix(integrationDegree+1, integrationDegree+1);
      for(UInt i=0; i<W.at(k+half).columns(); i++)
        for(UInt n=0; n<W.at(k+half).rows(); n++)
          W.at(k+half)(n,i) = ((n==0) ? 1.0 : pow(static_cast<Double>(i)-half-k, n));
      inverse(W.at(k+half));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
/***********************************************/

IntegralEquation::Arc  IntegralEquation::integrateArc(const OrbitArc &orbit, const std::vector<Rotary3d> &rotEarth, GravityfieldPtr gradientfield,
                                                      const const_MatrixSlice &g, const const_MatrixSlice &accCalibration, Bool computeVelocity, Bool computeAcceleration) const
{
  try
  {
    UInt   posCount = orbit.size();
    Double T        = (orbit.at(posCount-1).time - orbit.at(0).time).seconds();

    // Gravitationsgradient = Tensor
    // (im raumfesten System)
    // -----------------------------
    std::vector<Matrix> tensor(posCount);
    for(UInt i=0; i<posCount; i++)
    {
      Time     time     = orbit.at(i).time;
      Vector3d posEarth = rotEarth.at(i).rotate(orbit.at(i).position);
      tensor.at(i)      = rotEarth.at(i).inverseRotate(gradientfield->gravityGradient(time, posEarth)).matrix();
    }

    Matrix A(3*posCount, 3*posCount + 6 + g.columns());
    MatrixSlice VPos        (A.column(0,            3*posCount));
    MatrixSlice VPosBoundary(A.column(3*posCount,   6));
    MatrixSlice vPos        (A.column(3*posCount+6, g.columns()));

    // Integral mit dem Kern K(tau, tau')
    // und Drehung in das erdfeste System
    // ----------------------------------
    if(IntegrationPos.rows()!=3*posCount)
      IntegrationPos = integrationMatrixPosition(T, posCount);
    for(UInt i=0; i<posCount; i++)
      matMult(1., IntegrationPos.slice(3, 3*i, 3*posCount-3, 3), rotEarth.at(i).matrix().trans(), VPos.slice(3, 3*i, 3*posCount-3, 3));

    // Referenzpositionen und Randwerte
    // --------------------------------
    Vector   vBound(3*posCount);
    Vector   vGpsPos(3*posCount);
    Vector3d posStart = orbit.at(0).position;
    Vector3d posEnd   = orbit.at(posCount-1).position;
    for(UInt i=0; i<posCount; i++)
    {
      Double tau = static_cast<Double>(i)/(posCount-1.0);

      vGpsPos(3*i+0) = orbit.at(i).position.x();
      vGpsPos(3*i+1) = orbit.at(i).position.y();
      vGpsPos(3*i+2) = orbit.at(i).position.z();

      vBound(3*i+0) = (1-tau) * posStart.x() + tau * posEnd.x();
      vBound(3*i+1) = (1-tau) * posStart.y() + tau * posEnd.y();
      vBound(3*i+2) = (1-tau) * posStart.z() + tau * posEnd.z();
    }

    // Referenzbeschleunigungen integrieren
    // Reference = (Gravity*g-(vPos-vBound))
    // -------------------------------------
    matMult(1., VPos, g, vPos);
    vBound -= vGpsPos;
    for(UInt j=0; j<g.columns(); j++)
      vPos.column(j) += vBound;

    // partielle Ableitungen nach den Randwerten berechnen (Boundary)
    // --------------------------------------------------------------
    for(UInt i=0; i<posCount; i++)
    {
      Double tau = i/(posCount-1.0);
      VPosBoundary(3*i+0, 0) = VPosBoundary(3*i+1, 1) = VPosBoundary(3*i+2, 2) = (1-tau); // r_anf
      VPosBoundary(3*i+0, 3) = VPosBoundary(3*i+1, 4) = VPosBoundary(3*i+2, 5) = tau;     // r_end
    }

    // Indirekten Effekt beruecksichtigen
    // (Aenderung der Position aendert die Feldstaerke entlang der Bahn)
    // (I-K*tensor)^-1
    // ----------------------------
    // Einheitsmatrix
    Matrix Inv(3*posCount, 3*posCount);
    for(UInt i=0; i<Inv.rows(); i++)
      Inv(i,i) = 1.0;
    // I +- K*Tensor
    for(UInt i=0; i<posCount; i++)
      matMult(-1., IntegrationPos.slice(3,3*i,3*posCount-3,3), tensor.at(i), Inv.slice(3,3*i,3*posCount-3,3));

    // (I-K*tensor)^-1(Gravity*g+vBound-vPos)+pos
    // wobei A = (Gravity, Boundary, Reference)
    // und Reference = (IntegrationPos*g+Randwerte-pos)
    // ------------------------------------------------
    solveInPlace(Inv, A);

    // Referenzpositionen
    // Reference = (I-K*tensor)^-1(Gravity*g+vBound-vPos)+pos
    // --------------------
    Matrix vDeltaPos = vPos;
    for(UInt j=0; j<vPos.columns(); j++)
      vPos.column(j) += vGpsPos;

    Arc arc;
    arc.vPos         = vPos;
    arc.VPos         = VPos;
    arc.VPosBoundary = VPosBoundary;

    // =============================================================================

    // Beobachtungsgleichungen für Vektorielle Beschleunigungen aus Beschleunigungen
    // -----------------------------------------------------------------------------
    if(computeVelocity || computeAcceleration)
    {
      arc.VAcc = Matrix(3*posCount, 3*posCount);
      // gedrehte Einheitsmatrix
      for(UInt i=0; i<posCount; i++)
        copy(rotEarth.at(i).matrix().trans(), arc.VAcc.slice(3*i,3*i,3,3));

      // Referenzbeschleuigungen
      // -----------------------
      arc.vAcc = arc.VAcc * g;

      // indirekter Effekt
      for(UInt i=0; i<posCount; i++)
        matMult(1., tensor.at(i), vDeltaPos.row(3*i,3), arc.vAcc.row(3*i,3));

      // indirekter Effekt
      // -----------------
      for(UInt i=0; i<posCount; i++)
        matMult(1., tensor.at(i), arc.VPos.row(3*i,3), arc.VAcc.row(3*i,3));

      // Randwerte (nur indirekter Effekt)
      // ---------------------------------
      arc.VAccBoundary = Matrix(3*posCount, 6);
      for(UInt i=0; i<posCount; i++)
        matMult(1., tensor.at(i), arc.VPosBoundary.row(3*i,3), arc.VAccBoundary.row(3*i,3));
    }

    // =============================================================================

    // Beobachtungsgleichungen für Vektorielle Geschwindigkeiten aus Beschleunigungen
    // ------------------------------------------------------------------------------
    if(computeVelocity)
    {
      // Integral mit dem Kern (d/dt K(tau, tau'))
      if(IntegrationVel.rows()!=3*posCount)
        IntegrationVel = integrationMatrixVelocity(T, posCount);
      arc.VVel = IntegrationVel * arc.VAcc;

      // Randwerte
      // ---------
      arc.VVelBoundary = Matrix(3*posCount, 6);
      for(UInt i=0; i<posCount; i++)
      {
        arc.VVelBoundary(3*i+0, 0) = arc.VVelBoundary(3*i+1, 1) = arc.VVelBoundary(3*i+2, 2) = -1/T; // r_anf
        arc.VVelBoundary(3*i+0, 3) = arc.VVelBoundary(3*i+1, 4) = arc.VVelBoundary(3*i+2, 5) =  1/T; // r_end
      }

      // indirekter Effekt
      matMult(1., IntegrationVel, arc.VAccBoundary, arc.VVelBoundary);

      // Referenzgeschwindigkeiten
      // -------------------------
      arc.vVel = IntegrationVel * arc.vAcc;

        // Randwerte - pos
      for(UInt i=0; i<posCount; i++)
        for(UInt j=0; j<arc.vVel.columns(); j++)
        {
          arc.vVel(3*i+0,j) += (arc.vPos(3*posCount-3,j) - arc.vPos(0,j))/T; //(1/T) * pos1End.x() - (1/T) * pos1Start.x();
          arc.vVel(3*i+1,j) += (arc.vPos(3*posCount-2,j) - arc.vPos(1,j))/T; //(1/T) * pos1End.y() - (1/T) * pos1Start.y();
          arc.vVel(3*i+2,j) += (arc.vPos(3*posCount-1,j) - arc.vPos(2,j))/T; //(1/T) * pos1End.z() - (1/T) * pos1Start.z();
        }
    }

    // =============================================================================

    // Fit integrated orbit to given orbit
    // -----------------------------------
    Matrix B(3*posCount, 6+accCalibration.columns());
    copy(VPosBoundary, B.column(0,6));
    if(accCalibration.size())
      matMult(1., VPos, accCalibration, B.column(6,accCalibration.columns()));
    Vector x = leastSquares(B, vDeltaPos);

    if(arc.vPos.size()) matMult(-1, arc.VPosBoundary, x.row(0,6), arc.vPos);
    if(arc.vVel.size()) matMult(-1, arc.VVelBoundary, x.row(0,6), arc.vVel);
    if(arc.vAcc.size()) matMult(-1, arc.VAccBoundary, x.row(0,6), arc.vAcc);

    if(accCalibration.size())
    {
      Vector y = accCalibration*x.row(6,accCalibration.columns());
      if(arc.vPos.size()) matMult(-1, arc.VPos, y, arc.vPos);
      if(arc.vVel.size()) matMult(-1, arc.VVel, y, arc.vVel);
      if(arc.vAcc.size()) matMult(-1, arc.VAcc, y, arc.vAcc);
    }

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// interpolate to observations times
// ---------------------------------
void IntegralEquation::interpolateArc(const std::vector<OrbitArc> &pod, const OrbitArc &orbit, const Arc &arc, Matrix &l, Matrix &VPos, Matrix &VPosBoundary) const
{
  try
  {
    // observations
    l = Matrix(3*pod.at(0).size(), arc.vPos.columns());
    for(UInt j=0; j<l.columns(); j++)
      for(UInt k=0; k<pod.at(j).size(); k++)
      {
        l(3*k+0,j) = pod.at(j).at(k).position.x();
        l(3*k+1,j) = pod.at(j).at(k).position.y();
        l(3*k+2,j) = pod.at(j).at(k).position.z();
      }

    Polynomial polynomial(orbit.times(), interpolationDegree);
    auto timesPod = pod.at(0).times();
    l            -= polynomial.interpolate(timesPod, arc.vPos,         3);
    VPos          = polynomial.interpolate(timesPod, arc.VPos,         3);
    VPosBoundary  = polynomial.interpolate(timesPod, arc.VPosBoundary, 3);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Matrix IntegralEquation::integrationMatrixPosition(Double T, UInt posCount) const
{
  try
  {
    Matrix D(3*posCount, 3*posCount);
    UInt   unt = posCount-1;    // Anzahl Integrationsintervalle
    UInt   integrationDegree = W.size();

    for(UInt i=1; i<posCount-1; i++)
    {
      Double tau = i/static_cast<Double>(unt);
      Vector factor(posCount);

      UInt nr        = 0;
      UInt nullpunkt = 0;

      for(UInt k=0; k<unt; k++)    // Integrationsintervall
      {
        for(UInt j=0; j<W.at(nr).columns(); j++) // Position im Interpolationspolynom
          for(UInt n=0; n<W.at(nr).rows(); n++)  // Grad im Interpolationspolynoms
            if(i>k)
              factor(k+j-nullpunkt) += -T*T/(unt*unt) * (1-tau) * (1./(n+2.) + k /(n+1.))      * W.at(nr)(j,n);
            else
              factor(k+j-nullpunkt) +=  T*T/(unt*unt) *   tau *   (1./(n+2.) - (unt-k)/(n+1.)) * W.at(nr)(j,n);

        // Polynommatrix auswaehlen
        if(k<integrationDegree/2)
          nr++, nullpunkt++;
        if(k>=unt-integrationDegree/2-1)
          nr++, nullpunkt++;
      }

      for(UInt k=0; k<posCount; k++)
        D(3*i+0,3*k+0) = D(3*i+1,3*k+1) = D(3*i+2,3*k+2) = factor(k);
    }
    return D;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix IntegralEquation::integrationMatrixVelocity(Double T, UInt posCount) const
{
  try
  {
    Matrix D(3*posCount, 3*posCount);
    UInt   unt = posCount-1;    // Anzahl Integrationsintervalle
    UInt   integrationDegree = W.size();

    for(UInt i=0; i<posCount; i++)
    {
      Vector factor(posCount);

      UInt nr        = 0;
      UInt nullpunkt = 0;

      for(UInt k=0; k<unt; k++)    // Integrationsintervall
      {
        for(UInt j=0; j<W.at(nr).columns(); j++) // Position im Interpolationspolynom
          for(UInt n=0; n<W.at(nr).rows(); n++)  // Grad im Interpolationspolynoms
            if(i>k)
              factor(k+j-nullpunkt) += T/(unt*unt) * (1./(n+2.) + k /(n+1.))      * W.at(nr)(j,n);
            else
              factor(k+j-nullpunkt) += T/(unt*unt) * (1./(n+2.) - (unt-k)/(n+1.)) * W.at(nr)(j,n);

        // Polynommatrix auswaehlen
        if(k<integrationDegree/2)
          nr++, nullpunkt++;
        if(k>=unt-integrationDegree/2-1)
          nr++, nullpunkt++;
      }

      for(UInt k=0; k<posCount; k++)
        D(3*i+0,3*k+0) = D(3*i+1,3*k+1) = D(3*i+2,3*k+2) = factor(k);
    }
    return D;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
