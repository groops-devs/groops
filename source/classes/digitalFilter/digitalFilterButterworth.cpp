/***********************************************/
/**
* @file digitalFilterButterworth.cpp
*
* @brief Butterworth filter.
*
* @author Andreas Kvas
* @date 2017-02-07
*
*/
/***********************************************/

#include "base/import.h"
#include "base/fourier.h"
#include "classes/digitalFilter/digitalFilterButterworth.h"

/***********************************************/

DigitalFilterButterworth::DigitalFilterButterworth(Config& config)
{
  try
  {
    UInt   order;       // filter order
    UInt   nparams = 1; // number of parameters in adjustment
    Double Wn, Vn;      // normalized cutoff frequency
    std::function<std::vector<std::complex<Double>>(const std::vector<std::complex<Double>>&)> frequencyMapping;
    std::string choice;

    readConfig(config, "order", order, Config::MUSTSET, "", "filter order");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", "filter type"))
    {
      if(readConfigChoiceElement(config, "lowpass",  choice))
      {
        readConfig(config, "Wn", Wn, Config::MUSTSET, "", "normalized cutoff frequency (f_c / f_nyq)");
        frequencyMapping = std::bind(lowpassFrequencyScaling, std::placeholders::_1, Wn);
        nparams = order;
      }
      if(readConfigChoiceElement(config, "highpass", choice))
      {
        readConfig(config, "Wn", Wn, Config::MUSTSET, "", "normalized cutoff frequency (f_c / f_nyq)");
        frequencyMapping = std::bind(highpassFrequencyScaling, std::placeholders::_1, Wn);
        nparams = order;
      }
      if(readConfigChoiceElement(config, "bandpass", choice))
      {
        readConfig(config, "Wn1", Wn, Config::MUSTSET, "", "lower normalized cutoff frequency (f_c / f_nyq)");
        readConfig(config, "Wn2", Vn, Config::MUSTSET, "", "upper normalized cutoff frequency (f_c / f_nyq)");
        frequencyMapping = std::bind(bandpassFrequencyScaling, std::placeholders::_1, Wn, Vn);
        nparams = 2*order;
      }
      if(readConfigChoiceElement(config, "bandstop", choice))
      {
        readConfig(config, "Wn1", Wn, Config::MUSTSET, "", "lower normalized cutoff frequency (f_c / f_nyq)");
        readConfig(config, "Wn2", Vn, Config::MUSTSET, "", "upper normalized cutoff frequency (f_c / f_nyq)");
        frequencyMapping = std::bind(bandstopFrequencyScaling, std::placeholders::_1, Wn, Vn);
        nparams = 2*order;
      }
      endChoice(config);
    }
    readConfig(config, "backwardDirection", backward,          Config::DEFAULT,  "0", "apply filter in backward direction");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    Vector f = Fourier::frequencies(1022, 1.);
    std::vector<std::complex<Double>> z(f.rows()); // evaluation frequencies for transfer function

    for(UInt k = 0; k<f.rows(); k++)
      z.at(k) = std::exp(std::complex<Double>(0.,1.) * f(k) * 2.0*PI);

    // compute impulse response from transfer functions
    std::vector< std::complex<Double> > H =  transferFunction(z, frequencyMapping, order);
    H.back() = std::complex<Double>(H.back().real(), 0); // IFFT of real symmetric signal
    Vector h = Fourier::synthesis(H, TRUE/*even*/);

    // set up least squares adjustment
    Matrix H0(h.rows(), nparams);
    for(UInt k=0; k<nparams; k++) // columnwise circulant shift
    {
      copy(h.row(0, h.rows()-k-1), H0.slice(k+1, k, h.rows()-k-1, 1)); // diagonal and below
      copy(h.row(h.rows()-k-1, k+1), H0.slice(0, k, k+1, 1));          // above diagonal
    }

    // solve least squares via QR decomposition
    Vector l = h.row(nparams+1, h.rows()-nparams-1);
    Matrix A = H0.row(nparams+1, H0.rows()-nparams-1);
    Vector a_hat = leastSquares(A, l);

    // set numerator and denominator polynomials
    bnStartIndex = 0; // causal filter
    bn =  h.row(0, nparams+1);
    matMult(-1.0, H0.row(0, nparams+1), a_hat, bn);

    an = Vector(a_hat.rows()+1);
    an(0) = 1.0;
    axpy(-1.0, a_hat, an.row(1, a_hat.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::transferFunction(const std::vector<std::complex<Double> > &z, std::function<std::vector<std::complex<Double>>(const std::vector<std::complex<Double>>&)> &mapping, UInt order, Double g0)
{
  std::vector< std::complex<Double> > s  = mapping(z); // map z to s-plane
  std::vector< std::complex<Double> > Bn = butterworthPolynomial(s, order);

  std::vector< std::complex<Double> > H(z.size());
  for(UInt l = 0; l<H.size(); l++)
    H.at(l) = g0/Bn.at(l);

  return H;
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::butterworthPolynomial(const std::vector< std::complex<Double> > &s, UInt order)
{
  std::vector< std::complex<Double> > Bn(s.size(), std::complex<Double>(1.0, 0.0));
  if((order%2) != 0) // odd orders
    for(UInt l = 0; l<Bn.size(); l++)
      Bn.at(l) = s.at(l) + 1.0;

  for(UInt k = 1; k<order/2+1; k++)
    for(UInt l = 0; l<Bn.size(); l++)
      Bn.at(l) *= (s.at(l)*s.at(l) - 2.0*s.at(l)*std::cos( (2.0*k+order-1.0)/(2.0*order)*PI) + 1.0);

  return Bn;
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::lowpassFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn)
{
  std::vector< std::complex<Double> > s(z.size());
  for(UInt l = 0; l<s.size(); l++)
    s.at(l) = 1.0/std::tan(0.5*PI*Wn)*(z.at(l)-1.0)/(z.at(l)+1.0);

  return s;
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::highpassFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn)
{
  std::vector< std::complex<Double> > s(z.size());
  for(UInt l = 0; l<s.size(); l++)
    s.at(l) = std::tan(0.5*PI*Wn)*(z.at(l)+1.0)/(z.at(l)-1.0);

  return s;
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::bandpassFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn, Double Vn)
{
  Double omega = std::sqrt(std::tan(0.5*PI*Wn)*std::tan(0.5*PI*Vn));
  Double Q = omega/(std::tan(0.5*PI*Vn) - std::tan(0.5*PI*Wn));

  std::vector< std::complex<Double> > s(z.size());
  for(UInt l = 0; l<s.size(); l++)
    s.at(l) = Q*((z.at(l)-1.0)/(z.at(l)+1.0)/omega + omega*(z.at(l)+1.0)/(z.at(l)-1.0) );

  return s;
}

/***********************************************/

std::vector< std::complex<Double> > DigitalFilterButterworth::bandstopFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn, Double Vn)
{
  Double omega = std::sqrt(std::tan(0.5*PI*Wn)*std::tan(0.5*PI*Vn));
  Double Q = omega/(std::tan(0.5*PI*Vn) - std::tan(0.5*PI*Wn));

  std::vector< std::complex<Double> > s(z.size());
  for(UInt l = 0; l<s.size(); l++)
    s.at(l) = 1.0/(Q*((z.at(l)-1.0)/(z.at(l)+1.0)/omega + omega*(z.at(l)+1.0)/(z.at(l)-1.0)));

  return s;
}

/***********************************************/

