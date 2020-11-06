/***********************************************/
/**
* @file digitalFilter2FrequencyResponse.cpp
*
* @brief Amplitude and phase response of a filter cascade
*
* @author Andreas Kvas
* @date 2017-02-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute amplitude-, phase-, group delay and frequency response of a \configClass{digitalFilter}{digitalFilterType} cascade.
The \configFile{outputfileResponse}{matrix} is a matrix with following columns:
freq $[Hz]$, ampl, phase $[rad]$, group delay $[-]$, real, imag.

When \config{unwrapPhase} is set to true, $2\pi$ jumps of the phase response are removed before writing the output to file.

The response of the filter cascade is given by the product of each individual frequency response:
\begin{equation}
  H(f) = \prod_f H_j(f).
\end{equation}
Amplitude and phase response are computed from the frequency response via
\begin{equation}
  A(f) = |H(f)| \hspace{5pt}\text{and}\hspace{5pt} \Phi(f) = \arctan \frac{\mathcal{I}(H(f))}{\mathcal{R}(H(f))}.
\end{equation}
The group delay is computed by numerically differentiating the phase response
\begin{equation}
  \tau_g(f_k) = \frac{1}{2} \left[\frac{\Phi(f_k) - \Phi(f_{k-1})}{2\pi(f_k-f_{k-1})} + \frac{\Phi(f_{k+1}) - \Phi(f_{k})}{2\pi(f_{k+1}-f_{k})}\right] \approx \frac{d\Phi}{df}\frac{df}{d\omega}.
\end{equation}
The frequency vector for a \config{length} $N$ and a \config{sampling} $\Delta t$ is given by
\begin{equation}
  f_k = \frac{k}{N \Delta t}, \hspace{15pt} k \in \{0, \dots, \left\lfloor\frac{N+2}{2}\right\rfloor-1\}.
\end{equation}

See also \program{DigitalFilter2ImpulseResponse}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Amplitude and phase response of a filter cascade.
* @ingroup programsGroup */
class DigitalFilter2FrequencyResponse
{
private:
  Vector groupDelay(const Vector &f, const Vector &phase);
  Vector unwrap(const Vector &phase);
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(DigitalFilter2FrequencyResponse, SINGLEPROCESS, "amplitude and phase response of a filter cascade", Misc)

/***********************************************/

void DigitalFilter2FrequencyResponse::run(Config& config)
{
  try
  {
    FileName         fileNameOut;
    UInt             length;
    Double           sampling;
    DigitalFilterPtr filter;
    Bool             skipZero = FALSE;
    Bool             unwrapPhase = FALSE;

    readConfig(config, "outputfileResponse", fileNameOut, Config::MUSTSET,  "",    "columns: freq [Hz], ampl, phase [rad], group delay [-], real, imag");
    readConfig(config, "digitalFilter",      filter,      Config::MUSTSET,  "",    "");
    readConfig(config, "length",             length,      Config::DEFAULT,  "512", "length of the data series in time domain");
    readConfig(config, "sampling",           sampling,    Config::DEFAULT,  "1.0", "sampling to determine frequency [seconds]");
    readConfig(config, "skipZeroFrequency",  skipZero,    Config::DEFAULT,  "0",   "omit zero frequency when writing to file");
    readConfig(config, "unwrapPhase",        unwrapPhase, Config::DEFAULT,  "0",   "unwrap phase response");
    if(isCreateSchema(config)) return;

    logStatus<<"compute frequency response"<<Log::endl;
    Matrix A((length+2)/2, 6);
    copy(Fourier::frequencies(length, sampling), A.column(0));

    auto F = filter->frequencyResponse(length);
    Vector amplitude, phase;
    Fourier::complex2AmplitudePhase(F, amplitude, phase);
    Vector unwrappedPhase = unwrap(phase);
    Vector grpDel = groupDelay(A.column(0), unwrappedPhase);
    if(unwrapPhase)
      phase = unwrappedPhase;

    copy(amplitude,  A.column(1));
    copy(phase,      A.column(2));
    copy(grpDel,     A.column(3));
    for(UInt i=0; i<A.rows(); i++)
    {
      A(i, 4) = F.at(i).real();
      A(i, 5) = F.at(i).imag();
    }

    logStatus<<"write response to <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector DigitalFilter2FrequencyResponse::unwrap(const Vector &phase)
{
  Vector uPhase = phase;

  for(UInt n = 1; n<uPhase.rows(); n++)
  {
    while( (uPhase(n) - uPhase(n-1)) <= PI)
      uPhase(n) += 2*PI;

    while( (uPhase(n) - uPhase(n-1)) >= PI)
      uPhase(n) -= 2*PI;
  }

  return uPhase;
}

/***********************************************/

Vector DigitalFilter2FrequencyResponse::groupDelay(const Vector &f, const Vector &uPhase)
{
  Vector grpDel(uPhase.rows());

  grpDel(0) = -(uPhase(1)-uPhase(0))/(f(1)-f(0));
  for(UInt n = 1; n<uPhase.rows()-1; n++)
    grpDel(n) = -( (uPhase(n)-uPhase(n-1))/(f(n)-f(n-1)) + (uPhase(n+1)-uPhase(n))/(f(n+1)-f(n)))*0.5;
  grpDel(grpDel.rows()-1) = -(uPhase(grpDel.rows()-1)-uPhase(grpDel.rows()-2))/(f(grpDel.rows()-1)-f(grpDel.rows()-2));

  return grpDel/(2*PI);
}

/***********************************************/
