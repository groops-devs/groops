/***********************************************/
/**
* @file fourier.cpp
*
* @brief FFT-functions.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2004-10-25
*/
/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/fourier.h"

/***********************************************/

static inline void butterfly2(std::complex<Double> *f, const std::vector<std::complex<Double>> &twiddles, UInt step, UInt m)
{
  for(UInt i=0; i<m; i++)
  {
    const auto t = f[i+m] * twiddles[i*step];
    f[i+m]  = f[i]-t;
    f[i]   += t;
  }
}

/***********************************************/

static inline void butterfly3(std::complex<Double> *f, const std::vector<std::complex<Double>> &twiddles, UInt step, UInt m)
{
  const Double epi3 = twiddles[step*m].imag();
  for(UInt i=0; i<m; i++)
  {
    const auto t1 = f[i+  m] * twiddles[  i*step];
    const auto t2 = f[i+2*m] * twiddles[2*i*step];
    const auto t3 = t1 + t2;
    const auto t0 = std::complex<Double>(-epi3*(t1.imag()-t2.imag()), epi3*(t1.real()-t2.real()));
    f[i+m]    = f[i] - 0.5*t3;
    f[i+2*m]  = f[i+m] - t0;
    f[i]     += t3;
    f[i+m]   += t0;
  }
}

/***********************************************/

static inline void butterfly4(std::complex<Double> *f, const std::vector<std::complex<Double>> &twiddles, UInt step, UInt m, Bool inverse)
{
  for(UInt i=0; i<m; i++)
  {
    const auto t0 = f[i+  m] * twiddles[  i*step];
    const auto t1 = f[i+2*m] * twiddles[2*i*step];
    const auto t2 = f[i+3*m] * twiddles[3*i*step];
    const auto t3 = t0   + t2;
    const auto t4 = f[i] - t1;
    const auto t5 = std::complex<Double>(+t0.imag()-t2.imag(), -t0.real()+t2.real()) * ((inverse) ? -1. : 1.);
    f[i]    += t1;
    f[i+2*m] = f[i] - t3;
    f[i]    += t3;
    f[i+m]   = t4 + t5;
    f[i+3*m] = t4 - t5;
  }
}

/***********************************************/

static inline void butterfly5(std::complex<Double> *f, const std::vector<std::complex<Double>> &twiddles, UInt step, UInt m)
{
  const std::complex<Double> ya = twiddles[step*m];
  const std::complex<Double> yb = twiddles[step*2*m];
  for(UInt i=0; i<m; i++)
  {
    const auto t0  = f[i+1*m] * twiddles[  i*step];
    const auto t1  = f[i+2*m] * twiddles[2*i*step];
    const auto t2  = f[i+3*m] * twiddles[3*i*step];
    const auto t3  = f[i+4*m] * twiddles[4*i*step];
    const auto t4  = t0 + t3;
    const auto t5  = t1 + t2;
    const auto t6  = t1 - t2;
    const auto t7  = t0 - t3;
    const auto t8  = f[i] + t4 * ya.real() + t5 * yb.real();
    const auto t9  = f[i] + t4 * yb.real() + t5 * ya.real();
    const auto t10 = std::complex<Double>(+t7.imag(), -t7.real()) * ya.imag() + std::complex<Double>(t6.imag(), -t6.real()) * yb.imag();
    const auto t11 = std::complex<Double>(-t7.imag(), +t7.real()) * yb.imag() + std::complex<Double>(t6.imag(), -t6.real()) * ya.imag();
    f[i]     += t4 + t5;
    f[i+1*m]  = t8 - t10;
    f[i+2*m]  = t9 + t11;
    f[i+3*m]  = t9 - t11;
    f[i+4*m]  = t8 + t10;
  }
}

/***********************************************/

// perform the butterfly for one stage of a mixed radix FFT.
static inline void butterfly(std::complex<Double> *f, const std::vector<std::complex<Double>> &twiddles, UInt step, UInt m, UInt p)
{
  std::vector<std::complex<Double>> t(p);
  for(UInt u=0; u<m; u++)
  {
    for(UInt i=0; i<p; i++)
      t[i] = f[i*m+u];
    for(UInt i=0; i<p; i++)
    {
      f[i*m+u] = t[0];
      for(UInt k=1; k<p; k++)
        f[i*m+u] += t[k] * twiddles[(k*step*(i*m+u))%twiddles.size()];
    }
  }
}

/***********************************************/

static inline std::vector<UInt> computeRadix(UInt n)
{
  std::vector<UInt> factors;
  UInt p = 4; // factor out powers of 4, powers of 2, then any remaining primes
  do
  {
    while(n % p)
    {
      if(p == 4)      p  = 2;
      else if(p == 2) p  = 3;
      else            p += 2;
    }
    n /= p;
    factors.push_back(p);
    factors.push_back(n);
  }
  while(n > 1);

  return factors;
}

/***********************************************/

// recursive call: DFT of size m*p performed by doing p instances of smaller DFTs of size m, each one takes a decimated version of the input
static inline void recursiveFft(Bool inverse, std::complex<Double> *f, const std::complex<Double> *input, const UInt *factors, const std::vector<std::complex<Double>> &twiddles, UInt step)
{
  const UInt p = *(factors++); // the radix
  const UInt m = *(factors++); // stage's fft length/p

  if(m>1)
  {
    for(UInt i=0; i<p; i++)
      recursiveFft(inverse, f+(i*m), input+(i*step), factors, twiddles, p*step);
  }
  else
    for(UInt i=0; i<p; i++)
      f[i] = input[i*step];

  // recombine the p smaller DFTs
  switch(p)
  {
    case 2:  butterfly2(f, twiddles, step, m);          break;
    case 3:  butterfly3(f, twiddles, step, m);          break;
    case 4:  butterfly4(f, twiddles, step, m, inverse); break;
    case 5:  butterfly5(f, twiddles, step, m);          break;
    default: butterfly (f, twiddles, step, m, p);       break;
  }
}

/***********************************************/
/***********************************************/

std::vector<std::complex<Double>> Fourier::fft(const Vector &data)
{
  try
  {
    const UInt count = data.rows();

    std::vector<std::complex<Double>> input(count);
    for(UInt i=0; i<count; i++)
      input[i] = data(i);

    // compute twiddle factors
    std::vector<std::complex<Double>> twiddles(count);
    for(UInt i=0; i<twiddles.size(); i++)
      twiddles[i] = std::exp(std::complex<Double>(0, -2*PI*i/count));

    std::vector<std::complex<Double>> F(twiddles.size());
    std::vector<UInt> factors = computeRadix(twiddles.size());
    recursiveFft(FALSE/*inverse*/, F.data(), input.data(), factors.data(), twiddles, 1);
    F.resize((count+2)/2);
    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Fourier::synthesis(const std::vector<std::complex<Double>> &F, Bool countEven)
{
  try
  {
    const UInt count = 2*F.size() - (countEven ? 2 : 1);

    // extent input symmetric
    std::vector<std::complex<Double>> F2(count);
    F2[0] = F[0];
    for(UInt i=1; i<F.size(); i++)
    {
      F2[i]       = F[i];
      F2[count-i] = std::conj(F[i]);
    }

    // compute twiddle factors
    std::vector<std::complex<Double>> twiddles(count);
    for(UInt i=0; i<twiddles.size(); i++)
      twiddles[i] = std::exp(std::complex<Double>(0, 2*PI*i/count));

    std::vector<std::complex<Double>> F3(twiddles.size());
    std::vector<UInt> factors = computeRadix(twiddles.size());
    recursiveFft(TRUE/*inverse*/, F3.data(), F2.data(), factors.data(), twiddles, 1);

    Vector data(count);
    for(UInt i=0; i<F3.size(); i++)
      data(i) = (1./count)*F3[i].real();
    return data;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Fourier::frequencies(UInt count, Double dt)
{
  Vector frequencies((count+2)/2);
  for(UInt k=0; k<(count+2)/2; k++)
    frequencies(k) = static_cast<Double>(k)/(count*dt);
  return frequencies;
}

/***********************************************/

void Fourier::complex2AmplitudePhase(const std::vector<std::complex<Double>> &F, Vector &amplitude, Vector &phase)
{
  amplitude = Vector(F.size());
  phase     = Vector(F.size());
  for(UInt k = 0; k<F.size(); k++)
  {
    phase(k)     = std::atan2(F.at(k).imag(), F.at(k).real());
    amplitude(k) = std::abs(F.at(k));
  }
}

/***********************************************/
/***********************************************/

static void cosTransformation(Vector &x)
{
  const UInt n = x.rows();
  if(n<2)
    return;

  Double c = x(0)-x(n-1);
  x(0) += x(n-1);
  for(UInt i=1; i<n/2; i++)
  {
    c += 2*std::cos(PI*i/(n-1)) * (x(i)-x(n-1-i));
    const Double t1 = (x(i)+x(n-1-i));
    const Double t2 = (x(i)-x(n-1-i)) * 2 * std::sin(PI*i/(n-1));
    x(i)     = t1 - t2;
    x(n-1-i) = t1 + t2;
  }
  if(n%2)
    x(n/2) *= 2;

  const auto F = Fourier::fft(x.row(0, n-1));

  x(0) = F[0].real();
  x(1) = c;
  for(UInt i=1; i<n/2; i++)
  {
    x(2*i+0) = F[i].real();
    x(2*i+1) = x(2*i-1) - F[i].imag();
  }
  if(n%2)
    x(n-1) = F.back().real();
}

/***********************************************/

Vector Fourier::covariance2psd(const Vector &cov, Double dt)
{
  try
  {
    Vector psd = 2*dt*cov; // one sided PSD
    cosTransformation(psd);
    return psd;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Fourier::psd2covariance(const Vector &psd, Double dt)
{
  try
  {
    Vector cov = 0.25/(dt*(psd.rows()-1))*psd; // from  one sided PSD
    cosTransformation(cov);
    return cov;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
