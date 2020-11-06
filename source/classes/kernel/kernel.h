/***********************************************/
/**
* @file kernel.h
*
* @brief Isotropic harmonic integral kernels.
* E.g. the Abel-Poisson kernel for the continuation of the potential.
* Can also be used as Radial Basis Functions for splines and wavelets.
*
* @author Torsten Mayer-Guerr
* @date 2003-09-03
*
*/
/***********************************************/

#ifndef __GROOPS_KERNEL__
#define __GROOPS_KERNEL__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernel = R"(
\section{Kernel}\label{kernelType}
Kernel defines harmonic isotropic integral kernels $K$.
\begin{equation}
T(P) = \frac{1}{4\pi}\int_\Omega K(P,Q)\cdot f(Q)\,d\Omega(Q),
\end{equation}
where $T$ is the (disturbance)potential and $f$ is a functional on the spherical surface~$\Omega$.
The Kernel can be exapanded into a series of (fully normalized) legendre polynomials
\begin{equation}\label{eq.kernel}
K(\cos\psi,r,R) = \sum_n \left(\frac{R}{r}\right)^{n+1}
k_n\sqrt{2n+1}\bar{P}_n(\cos\psi).
\end{equation}
On the one hand the kernel defines the type of the functionals~$f$ that are measured
or have to be computed, e.g. gravity anomalies given by the Stokes-kernel.
On the other hand the kernel functions can be used as basis functions to represent
the gravity field, e.g. as spline functions or wavelets.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup kernelGroup Kernel
* @brief Isotropic harmonic integral kernels.
* @ingroup classesGroup
* The interface is given by @ref Kernel.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Kernel;
typedef std::shared_ptr<Kernel> KernelPtr;
class GravityfieldBase;

/***** CLASS ***********************************/

/** @brief Isotropic harmonic integral kernels.
* E.g. the Abel-Poisson kernel for the continuation of the potential.
* Can also be used as Radial Basis Functions for splines and wavelets.
* An Instance of this class can be created by @ref readConfig. */
class Kernel
{
public:
  /// Destructor.
  virtual ~Kernel() {}

  /** @brief Function value of the Kernel.
  * @param p Computational point.
  * @param q Source point. */
  virtual Double kernel(const Vector3d &p, const Vector3d &q) const;

  /** @brief Radial Derivative.
  * @f$ \frac{\partial K}{\partial r} @f$.
  * Applied to the computational point @a p.
  * @param p Computational point.
  * @param q Source point. */
  virtual Double radialDerivative(const Vector3d &p, const Vector3d &q) const;

  /** @brief Gradient of the Kernel.
  * @f$ \nabla K @f$.
  * Applied to the computational point @a p.
  * @param p Computational point.
  * @param q Source point. */
  virtual Vector3d gradient(const Vector3d &p, const Vector3d &q) const;

  /** @brief Gradient of the gradient (Tensor).
  * @f$ \nabla \nabla K @f$.
  * Applied to the computational point @a p.
  * @param p Computational point.
  * @param q Source point. */
  virtual Tensor3d gradientGradient(const Vector3d &p, const Vector3d &q) const;

  /** @brief Apply inverse kernel to a given Kernel.
  * Convolution of Kernel of this class with the passed @a kernel.
  * Example: Is this class the Stokes kernel, the gravity anomalies of the other kernel is returned:
  * @f[ -\frac{\partial K}{\partial r}-\frac{2K}{r} @f]
  * @param kernel This will be convoluted.
  * @param p Computational point.
  * @param q Source point. */
  virtual Double inverseKernel(const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;

  /** @brief Apply inverse kernel to a gravity field.
  * Example: Is this class the stokes kernel,
  * gravity anomalies are computed from the given gravity field:
  * @f[ -\frac{\partial T}{\partial r}-\frac{2T}{r} @f]
  * @param time For time variable gravity fields.
  * @param p Computational point.
  * @param field Apply kernel to this field. */
  virtual Double inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const;

  /** @brief Legendre coefficients.
  * The Kernel can be represented by a series of LegendrePolynomial:
  * @f[ K(\cos\psi,r,R) = \sum_n (R/r)^{n+1} k_n \sqrt{2n+1}P_n(\cos\psi)  @f]
  * This function returns the inverse coefficients @f$ k_n @f$.
  * Example: The Stokes kernel gives @f$ r/(n-1) @f$.
  * @param p Computational point.
  * @param degree maximum degree.
  * @return vector with size @a degree+1. */
  virtual Vector coefficients(const Vector3d &p, UInt degree) const=0;

  /** @brief Legendre coefficients of the inverse Kernel.
  * The inverse Kernel can be represented by a series of LegendrePolynomial:
  * @f[ K^(\cos\psi,r,R) = \sum_n (R/r)^{n+1} 1/k_n  \sqrt{2n+1}P_n(\cos\psi)   @f]
  * This function returns the inverse coefficients @f$ 1/k_n @f$.
  * Example: The Stokes kernel gives @f$ (n-1)/r @f$.
  * @param p Computational point.
  * @param degree maximum degree.
  * @param interior functional is below masses.
  * @return vector with size @a degree+1. */
  virtual Vector inverseCoefficients(const Vector3d &p, UInt degree, Bool interior=FALSE) const=0;

  /** @brief Minimum for bandlimited kernels.
  * Otherwise 0 is returned. */
  virtual UInt minDegree() const {return  0;}

  /** @brief Maximum degree for bandlimited kernels.
  * Otherwise INFINITYDEGREE is returned. */
  virtual UInt maxDegree() const {return INFINITYDEGREE;}

  /** @brief creates an derived instance of this class. */
  static KernelPtr create(Config &config, const std::string &name);

protected:
  Double   kernel          (Vector3d const &p, Vector3d const &q, const Vector &kn) const;
  Double   radialDerivative(Vector3d const &p, Vector3d const &q, const Vector &kn) const;
  Vector3d gradient        (Vector3d const &p, Vector3d const &q, const Vector &kn) const;
  Tensor3d gradientGradient(Vector3d const &p, Vector3d const &q, const Vector &kn) const;

  Vector computeFactors                   (Double r, Double R, const Vector &kn) const;
  Vector computeFactorsRadialDerivative   (Double r, Double R, const Vector &kn) const;
  Vector computeFactorsRadialDerivative2nd(Double r, Double R, const Vector &kn) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Kernel.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a kernel is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] kernel Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Kernel */
template<> Bool readConfig(Config &config, const std::string &name, KernelPtr &kernel, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
