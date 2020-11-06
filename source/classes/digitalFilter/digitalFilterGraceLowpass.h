/***********************************************/
/**
* @file digitalFilterGraceLowpass.h
*
* @brief Low pass filter as used in GRACE processing.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERGRACELOWPASS__
#define __GROOPS_DIGITALFILTERGRACELOWPASS__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterGraceLowpass = R"(
\subsection{GraceLowpass}
Low pass and differentation filter as used for GRACE KBR and ACC data in the Level1A processing.

\fig{!hb}{0.8}{DigitalFilter_graceLowpass}{fig:DigitalFilterGraceLowpass}{Amplitude response of the low pass filter used in the L1A processing.}
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Low pass filter as used in GRACE processing.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterGraceLowpass : public DigitalFilterARMA
{
public:
  DigitalFilterGraceLowpass(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterGraceLowpass::DigitalFilterGraceLowpass(Config &config)
{
  try
  {
    Double      fs, B, Tf;
    UInt        Nc;
    Double      f0;
    Bool        reduceFit;
    UInt        type = 0;
    std::string choice;

    readConfig(config, "rawDataRate",         fs,        Config::MUSTSET,  "10.0",    "sampling frequency in Hz (fs).");
    readConfig(config, "convolutionNumber",   Nc,        Config::DEFAULT,  "7",       "number of self convolutions of the filter kernel");
    readConfig(config, "fitInterval",         Tf,        Config::DEFAULT,  "70.7",    "length of the filter kernel [seconds]");
    readConfig(config, "lowPassBandwith",     B,         Config::DEFAULT,  "0.1",     "target low pass bandwidth");
    readConfig(config, "normFrequency",       f0,        Config::DEFAULT,  "0.37e-3", "norm filter at this frequency [Hz] (default: GRACE dominant (J2) signal frequency)");
    readConfig(config, "reduceQuadraticFit",  reduceFit, Config::DEFAULT,  "1",       "remove->filter->restore quadratic fit");
    if(readConfigChoice(config, "derivative", choice,    Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "derivative1st", choice, "range rate"))         type = 1;
      if(readConfigChoiceElement(config, "derivative2nd", choice, "range acceleration")) type = 2;
      endChoice(config);
    }
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    const Int NB = static_cast<Int>(round(B*Tf));   // number of frequency bins in the passband
    const Int Nf = static_cast<Int>(round(fs*Tf));  // number of raw data points in the fit intervall
    const Int Nh = (Nf-1)/2;

    // self convoluted (Nc times) boxcar in frequency domain -> discrete sinc
    Vector Fk(2*Nh+1);
    for(Int k=-Nh; k<=Nh; k++)
      for(Int m=-NB; m<=NB; m++)
        Fk(k+Nh) += (k!=m) ? (pow(sin(PI*(k-m)/Nc)/sin(PI*(k-m)/Nf), Nc)) : pow(Nf/Nc, Nc);

    // back transformation into time domain: discrete cos transformation
    Vector fn(2*Nh+1);
    for(Int n=-Nh; n<=Nh; n++)
      for(Int k=-Nh; k<=Nh; k++)
        if(type == 0)
          fn(n+Nh) += Fk(Nh-k) * cos(2*PI*k*n/Nf);
        else if(type == 1)
          fn(n+Nh) += Fk(Nh-k) * sin(2*PI*k*n/Nf) * (-2*PI*k/Tf); // 1st derivative
        else if(type == 2)
          fn(n+Nh) += Fk(Nh-k) * cos(2*PI*k*n/Nf) * (-2*PI*k/Tf)*(2*PI*k/Tf); // 2nd derivative

    // norm filter (preserve dominant (J2) signal frequency)
    Vector f0n(2*Nh+1);
    for(Int n=-Nh; n<=Nh; n++)
      for(Int k=-Nh; k<=Nh; k++)
        f0n(n+Nh) += Fk(Nh-k) * cos(2*PI*k*n/Nf);
    Double Fnorm = 0.0;
    for(Int n=-Nh; n<=Nh; n++)
      Fnorm += cos(2*PI*f0*n/fs) * f0n(Nh-n);
    fn *= 1/Fnorm;

    // Reduce Quadratic Fit
    // --------------------
    if(reduceFit)
    {
      Matrix Qt(3, fn.rows()); // = Q^T
      for(UInt z=0; z<Qt.columns(); z++)
      {
        const Double t = ((fn.rows()-1.)/2.-z)/fs;
        Qt(0,z) = 1.0;
        Qt(1,z) = t;
        Qt(2,z) = t*t;
      }

      // (Q^T*Q)^-1 * Q^T
      Matrix QtQ(Qt.rows(), Matrix::SYMMETRIC);
      rankKUpdate(1., Qt.trans(), QtQ);
      Matrix QtQQt = solve(QtQ,Qt);

      matMult(-1, QtQQt.trans(), (Qt*fn), fn); // remove
      if(type == 0)
        fn += QtQQt.trans().column(0);   // restore
      else if(type == 1)
        fn += QtQQt.trans().column(1);   // restore 1st derivative
      else if(type == 2)
        fn += 2*QtQQt.trans().column(2); // restore 2nd derivative
    }

    bn = fn;
    bnStartIndex = bn.rows()/2;
    an = Vector(1); an(0) = 1.0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
