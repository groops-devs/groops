/***********************************************/
/**
* @file tidesPole.h
*
* @brief Centrifugal effect of polar motion.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2014-05-23
*
*/
/***********************************************/

#ifndef __GROOPS_TIDESPOLE__
#define __GROOPS_TIDESPOLE__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesPole = R"(
\subsection{PoleTide}\label{tidesType:poleTide}
The potential coefficients of the solid Earth pole tide according to the
IERS2003 conventions are given by
\begin{equation}
\begin{split}
\Delta c_{21} &= s\cdot(m_1 + o\cdot m_2), \\
\Delta s_{21} &= s\cdot(m_2 - o\cdot m_1),
\end{split}
\end{equation}
with $s$ is the \config{scale}, $o$ is the \config{outPhase} and
$(m_1,m_2)$ are the wobble variables in seconds of arc.
They are related to the polar motion variables $(x_p,y_p)$ according to
\begin{equation}
\begin{split}
m_1 &=  (x_p - \bar{x}_p), \\
m_2 &= -(y_p - \bar{y}_p),
\end{split}
\end{equation}
The mean pole $(\bar{x}_p, \bar{y}_p)$ is approximated by a polynomial
read from \configFile{inputfileMeanPole}{meanPolarMotion}.

The displacment is calculated with
\begin{equation}
\begin{split}
S_r          &= -v\sin2\vartheta(m_1\cos\lambda+m_2\sin\lambda),\\
S_\vartheta &= -h\cos2\vartheta(m_1\cos\lambda+m_2\sin\lambda),\\
S_\lambda   &=  h\cos\vartheta(m_1\sin\lambda-m_2\cos\lambda),
\end{split}
\end{equation}
where $h$ is the \config{horizontalDisplacement}
and $v$ is the \config{verticalDisplacement}.

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "files/fileMeanPolarMotion.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Centrifugal effect of polar motion.
* @ingroup tidesGroup
* @see Tides */
class TidesPole : public TidesBase
{
  Double          scale,     outPhase;
  Double          horizontalDisplacement, verticalDisplacement;
  Double          factor;
  MeanPolarMotion meanPole;

  void pole(const Time &time, EarthRotationPtr rotation, Double &m1, Double &m2) const;

public:
  TidesPole(Config &config);

  SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
  Vector3d deformation(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                       Double gravity, const Vector &hn, const Vector &ln) const override;
  void     deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                       EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                       std::vector<std::vector<Vector3d>> &disp) const override;
};

/***********************************************/

#endif /* __GROOPS_TIDES__ */
