/***********************************************/
/**
* @file digitalFilterDerivative.h
*
* @brief Numerical differentiation using polynomial approximation.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERDERIVATIVE__
#define __GROOPS_DIGITALFILTERDERIVATIVE__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterDerivative = R"(
\subsection{Derivative}
Symmetric MA filter for numerical differentiation using polynomial approximation. The input time series is approximated by a moving polynomial of degree \config{polynomialDegree}, by solving
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
The filter coefficients for the $k$-th derivative are obtained by taking the appropriate row of the inverse coefficient matrix $\mathbf{W}$:
\begin{equation}
  b_n = \prod_{i=0}^{k-1} (k-i) \mathbf{w}_{2,:}.
\end{equation}
The \config{polynomialDegree} shoud be even if no phase shift should be introduced.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Numerical differentiation using polynomial approximation.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterDerivative : public DigitalFilterARMA
{
public:
  DigitalFilterDerivative(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterDerivative::DigitalFilterDerivative(Config &config)
{
  try
  {
    UInt   degree;
    UInt   derivative;
    Double dt = 1.0;

    readConfig(config, "polynomialDegree", degree,     Config::MUSTSET,  "8",   "degree of approximation polynomial");
    readConfig(config, "derivative",       derivative, Config::DEFAULT,  "1",   "take kth derivative");
    readConfig(config, "sampling",         dt,         Config::DEFAULT,  "1.0", "assumed time step between points");
    readConfig(config, "padType",          padType,    Config::MUSTSET,  "",    "");
    if(isCreateSchema(config)) return;

    Matrix F(degree+1, degree+1);
    for(UInt j=0; j<F.rows(); j++)
    {
      F(j,0) = 1.;
      for(UInt n=0; n<F.columns(); n++)
       F(j,n) = std::pow(0.5*degree-j, n);
    }
    inverse(F);

    Double factor = 1.0;
    for(UInt n=1; n<=derivative; n++)
      factor *= n/dt;

    an = Vector(1, 1.);
    bn = Vector(F.row(derivative).trans()*factor);
    bnStartIndex = bn.size()/2; // non-causal filter
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
