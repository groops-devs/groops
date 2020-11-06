/***********************************************/
/**
* @file digitalFilterIntegral.h
*
* @brief Numerical integration using polynomial approximation.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERINTEGRAL__
#define __GROOPS_DIGITALFILTERINTEGRAL__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterIntegral = R"(
\subsection{Integral}
Numerical integration using polynomial approximation.
The input time series is approximated by a moving polynomial of degree \config{polynomialDegree}
by solving
\begin{equation}
  \begin{bmatrix} x(t_k+\tau_0) \\ \vdots \\ x(t_k+\tau_M) \end{bmatrix}
  =
  \begin{bmatrix}
  1      & \tau_0 & \tau_0^2 & \cdots & \tau_0^M \\
  \vdots & \vdots & \vdots   &        & \vdots   \\
  1      & \tau_M & \tau_M^2 & \cdots & \tau_M^M \\
  \end{bmatrix}%^{-1}
  \begin{bmatrix}
  a_0 \\ \vdots \\ a_M
  \end{bmatrix}
  \qquad\text{with}\quad
  \tau_j =  (j-M/2)\cdot \Delta t,
\end{equation}
for each time step $t_k$ ($\Delta t$ is the \config{sampling} of the time series).
The numerical integral for each time step $t_k$ is approximated by the center interval of the estimated polynomial.

\fig{!hb}{0.7}{DigitalFilter_integral}{fig:DigitalFilterIntegral}{Numerical integration by polynomial approximation.}

\config{polynomialDegree} should be even to avoid a phase shift.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Numerical integration using polynomial approximation.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterIntegral : public DigitalFilterARMA
{
public:
  DigitalFilterIntegral(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterIntegral::DigitalFilterIntegral(Config &config)
{
  try
  {
    UInt   degree;
    Double dt;

    readConfig(config, "polynomialDegree", degree,  Config::MUSTSET,  "7",   "degree of approximation polynomial");
    readConfig(config, "sampling",         dt,      Config::DEFAULT,  "1.0", "assumed time step between points");
    readConfig(config, "padType",          padType, Config::MUSTSET,  "",    "");
    if(isCreateSchema(config)) return;

    Matrix F(degree+1, degree+1);
    Matrix D(degree+1, degree+1);
    Vector v(degree+1);
    for(UInt j=0; j<F.rows(); j++)
    {
      D(j, j) = 1.0/std::pow(dt, j);
      v(j) = dt/(j+1.);
      F(j,0) = 1.;
      for(UInt n=1; n<F.columns(); n++)
        F(j,n) = std::pow(0.5*degree-j, n);
    }
    inverse(F);

    bnStartIndex = bn.rows()/2;
    bn = Vector(degree+1);
    matMult(1.0, v.trans(), D*F, bn.trans());
    an = Vector(2);
    an(0) = 1.0; an(1) = -1.0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
