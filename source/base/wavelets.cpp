/***********************************************/
/**
* @file wavelets.cpp
*
* @brief Basic functions for wavelet decomposition
*
* @author Andreas Kvas
* @date 2017-06-28
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/wavelets.h"

/***********************************************/

namespace Wavelets
{

/***********************************************/

Vector lowpass(const Vector &w)
{
  Vector b(w.rows());
  for(UInt k=0; k<b.rows(); k++)
    b(k) = w(w.rows()-1-k);
  return b;
}

/***********************************************/

Vector highpass(const Vector &w)
{
  Vector b(w.rows());
  for(UInt k=0; k<b.rows(); k++)
    b(k) = w(k) * ((k%2) ? 1. : -1.);
  return b;
}

/***********************************************/

void halfbandFilter(const Matrix &input, const Vector &wl, Matrix &detailCoefficients, Matrix &approxCoefficients)
{
  // high/lowpass filters
  const Vector hn = lowpass(wl);
  const Vector gn = highpass(wl);

  if(hn.size() != gn.size())
    throw(Exception("Length of high- and low-pass filters do not match."));

  // filter order
  const UInt order = hn.rows();

  // symmetric padding
  Matrix padded(2*order+input.rows(), input.columns());
  copy(input, padded.row(order, input.rows()));
  for(UInt k = 0; k < order; k++)
  {
    copy(input.row(k), padded.row(order - 1 - k));
    copy(input.row(input.rows() - 1 - k), padded.row(input.rows() + order + k));
  }

  // initialize state
  Matrix stateLowpass(order, padded.columns());
  Matrix stateHighpass(order, padded.columns());

  // filter input signal
  Matrix outputLowpass(padded.rows(), padded.columns());
  Matrix outputHighpass(padded.rows(), padded.columns());

  for(UInt k = 0; k<padded.rows(); k++)
  {
    // yn = b0 * xn + s0
    axpy(hn(0), padded.row(k), outputLowpass.row(k));
    axpy(1.0, stateLowpass.row(0), outputLowpass.row(k));

    axpy(gn(0), padded.row(k), outputHighpass.row(k));
    axpy(1.0, stateHighpass.row(0), outputHighpass.row(k));

    // s_i = b_i+1 * xn + s_i+1 - a_i+1*yn
    for(UInt i = 0; i<order-1; i++)
    {
      copy(stateLowpass.row(i+1), stateLowpass.row(i));
      if(i+1<hn.rows())
        axpy(hn(i+1), padded.row(k), stateLowpass.row(i));

      copy(stateHighpass.row(i+1), stateHighpass.row(i));
      if(i+1<gn.rows())
        axpy(gn(i+1), padded.row(k), stateHighpass.row(i));
    }
  }

  // compute output
  const UInt outputSize = (padded.rows()-order-1)/2;
  detailCoefficients = Matrix(outputSize, padded.columns());
  approxCoefficients = Matrix(outputSize, padded.columns());

  for(UInt k = 0; k<outputSize ; k++)
  {
    copy(outputLowpass.row(order+1+2*k), approxCoefficients.row(k));
    copy(outputHighpass.row(order+1+2*k), detailCoefficients.row(k));
  }
}

/***********************************************/

std::vector<Matrix> waveletTransform(const Matrix &input, const Vector &wl, UInt maxLevel)
{
  maxLevel = std::min(maxLevel, static_cast<UInt>(std::log2(static_cast<Double>(input.rows())/static_cast<Double>(wl.rows()-1))));

  std::vector<Matrix> levels;

  Matrix detail, approx;
  approx = input;
  for(UInt k = 0; k<maxLevel; k++)
  {
    Matrix tmp(approx);
    halfbandFilter(tmp, wl, detail, approx);
    levels.push_back(detail);
  }
  levels.push_back(approx);

  return levels;
}

/**************************************************************************/

void discreteWaveletTransform(std::vector<Double> &signal, const Vector &wl, UInt J, std::vector<Double> &coefficients, std::vector<UInt> &flag, std::vector<UInt> &length)
{
  try
  {
    UInt temp_len = signal.size();
    if((temp_len % 2) != 0)
    {
      Double temp =signal.at(temp_len - 1);
      signal.push_back(temp);
      flag.push_back(1);
      temp_len++;
    }
    else
      flag.push_back(0);
    length.push_back(temp_len);
    flag.push_back(J);
    // flag contains symmetric extension length

    std::vector<Double> originalSignal, approxs, details;
    originalSignal = signal;

    // Storing Filter Values
    std::vector<Double> lpDecom,hpDecom,lpRecon,hpRecon;
    filtCoef(wl,lpDecom,hpDecom,lpRecon,hpRecon);

    for(UInt iter = 0; iter < J; iter++)
    {
      UInt factor = 2; // Downsampling Factor is 2
      UInt lf = lpDecom.size();
      symmExtention(signal,lf-1);

      std::vector<Double> cA;
      //Low Pass Branch Computation
      convfft(signal,lpDecom,cA);
      cA.erase(cA.begin(),cA.begin()+lf);
      cA.erase(cA.end()-lf+1,cA.end());
      downSampling(cA, factor, approxs);

      //High Pass Branch Computation
      std::vector<Double> cD;
      convfft(signal,hpDecom,cD);
      cD.erase(cD.begin(),cD.begin()+lf);
      cD.erase(cD.end()-lf+1,cD.end());
      downSampling(cD, factor, details);

      coefficients.insert(coefficients.begin(),details.begin(),details.end());
      length.insert(length.begin(),details.size());

      if(iter == J-1 )
      {
        coefficients.insert(coefficients.begin(),approxs.begin(),approxs.end());
        length.insert(length.begin(),approxs.size());
      }

      signal.clear();
      signal = approxs;
      approxs.clear();
      details.clear();
     }
     signal = originalSignal;
  }

  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/********************************************************************/

void inverseDiscreteWaveletTransform(std::vector<Double> &dwtop, const Vector &wl, const std::vector<UInt> &flag, std::vector<Double> &recoveredSignal, std::vector<UInt> &length)
{
  try
  {
    UInt J = flag.at(1);
    UInt lf;

    std::vector<Double> approxs;
    std::vector<Double> details;
    UInt appLen = length.at(0);
    UInt detLen = length.at(1);

    std::vector<Double>::iterator dwt;
    dwt = dwtop.begin();
    approxs.assign(dwt,dwtop.begin()+appLen);
    details.assign(dwtop.begin()+appLen, dwtop.begin()+ 2* appLen);

    for(UInt i = 0; i < J; i++)
    {
      UInt factor = 2; // Upsampling Factor
      // Storing Filter Values
      std::vector<Double> lpDecom,hpDecom, lpRecon, hpRecon;
      filtCoef(wl, lpDecom,hpDecom, lpRecon, hpRecon);
      lf = lpRecon.size();

      // Operations in the Low Frequency branch of the Synthesis Filter Bank
      std::vector<Double> X_lp;
      std::vector<Double> cA_up;
      upSampling(approxs, factor,cA_up);
      cA_up.pop_back();
      convfft(cA_up, lpRecon, X_lp);

      // Operations in the High Frequency branch of the Synthesis Filter Bank
      std::vector<Double> X_hp;
      std::vector<Double> cD_up;
      upSampling(details, factor, cD_up);
      cD_up.pop_back();
      convfft(cD_up, hpRecon, X_hp);

      appLen += detLen;
      recoveredSignal.resize(X_lp.size());
      transform (X_lp.begin(), X_lp.end(), X_hp.begin(), X_hp.begin(), [](Double i, Double j ){return i+j;});
      recoveredSignal = X_hp;
      recoveredSignal.erase(recoveredSignal.begin(),recoveredSignal.begin()+lf-2);
      recoveredSignal.erase(recoveredSignal.end()-(lf - 2),recoveredSignal.end());

      approxs.clear();
      details.clear();

      if(i < J - 1)
      {
        detLen = length.at(i+2);
        for (UInt l = 0; l < detLen;l++)
        {
          Double temp = dwtop.at(appLen + l);
          details.push_back(temp);
        }
      }
      approxs = recoveredSignal;
      for(UInt iter1 = 0; iter1 < (approxs.size() - detLen); iter1++)
        approxs.pop_back();
    }
    // Remove Padding
    UInt paddLength = flag.at(0);
    recoveredSignal.erase(recoveredSignal.end()- paddLength,recoveredSignal.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/********************************************/

void convfft(std::vector<Double> &a, std::vector<Double> &b, std::vector<Double> &c)
{
  try
  {
     UInt N = a.size() + b.size() - 1;
     std::vector<std::complex<Double> > inp, filt;
     for(UInt i=0; i < a.size(); i++)
     {
       Double temp = a.at(i);
       inp.push_back(std::complex<Double>(temp,0));
     }
     for(UInt i =0; i < b.size(); i++)
     {
       Double temp = b.at(i);
       filt.push_back(std::complex<Double>(temp,0));
     }
     //apply Fast Fourier Transform
     dfft(inp,1,N);
     dfft(filt,1,N);
     std::vector<std::complex<Double> > temp;
     UInt K=inp.size();
     for(UInt i =0; i < K; i++)
     {
       std::complex<Double> mult = inp.at(i)*filt.at(i);
       temp.push_back(mult);
     }
     dfft(temp, -1 , K);
     for(UInt i =0; i < N ; i++)
     {
       Double temp1 =real(temp.at(i));
       c.push_back(temp1);
     }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***************************************************************/

void dfft(std::vector<std::complex<Double> > &data, int sign,UInt N)
{
  try
  {
    Double pi = -PI;
    if(sign == 1 || sign == -1)
      pi = sign * pi;
    UInt len = data.size();
    std::vector<std::complex<Double> >::iterator it;
    it = data.end();
    if(len != N)
    {
      UInt al = N - len;
      data.insert(it, al, std::complex<Double> (0,0));
    }

    UInt K = (UInt) pow(2.0,ceil(log10(static_cast<Double>(N))/log10(2.0)));
    std::vector<std::complex<Double> >::iterator it1;
    it1 = data.end();
    if(N < K)
    {
      UInt al = K - N;
      data.insert(it1,al,std::complex<Double>(0,0));
      N = K;
    }

    bitreverse(data);

    for(UInt iter = 1; iter < N; iter <<=1)
    {
      const UInt step = iter << 1;
      const Double theta =  pi / Double(iter);
      Double wtemp = sin(theta * .5);

      //Multipliers
      Double wreal = -2 * wtemp * wtemp;
      Double wimag = sin(theta);

      //Factors
      Double wr = 1.0;
      Double wi = 0.0;

      //Iteration through two loops
      for(UInt m = 0; m < iter; m++)
      {
        //Iteration within m
        for(UInt i = m; i < N; i += step)
        {
          //jth position
          const UInt j = i + iter;
          Double tempr = wr * real(data.at(j)) - wi * imag(data.at(j));
          Double tempi = wr * imag(data.at(j)) + wi * real(data.at(j));
          std::complex<Double> temp(tempr,tempi);
          data.at(j)  = data.at(i) - temp;
          data.at(i) += temp;
        }

        wtemp = wr;
        wr += wr * wreal - wi * wimag;
        wi += wi * wreal + wtemp * wimag ;
      }
    }

    if(sign == -1)
    {
      Double scale = 1.0/Double(N);
      for (UInt i = 0; i < N; i++)
        data.at(i) *= scale;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/**************************************************************/

void bitreverse(std::vector<std::complex<Double> > &signal)
{
  try
  {
    UInt len = signal.size();
    UInt N   = (UInt) pow(2.0,ceil(log10(static_cast<Double>(len))/log10(2.0)));
    UInt rev = 0;
    //Processing Input Data
    for(UInt iter = 0; iter < N; ++iter)
    {
      if (rev > iter)
      {
        //Replacing current values with reversed values
        Double tempr = real(signal.at(rev));
        Double tempi = imag(signal.at(rev));
        std::complex<Double> temp(tempr,tempi);
        signal.at(rev)  = signal.at(iter);
        signal.at(iter) = temp;
      }
      //Using filter "filt" such that the value of reverse changes with each iteration
      UInt filt = N;
      while (rev & (filt >>= 1))
        rev &= ~filt;
      rev |= filt;
    }
  }

  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/**************************************************************/

void filtCoef(const Vector &wl, std::vector<Double> &lpDecom, std::vector<Double> &hpDecom, std::vector<Double> &lpRecon, std::vector<Double> &hpRecon)
{
  try
  {
    Vector gn =  highpass(wl);
    Vector hn =  lowpass (wl);
    UInt length = gn.size();
    for(UInt i = 0; i < length; i++)
    {
      lpDecom.push_back(hn(i));
      hpDecom.push_back(gn(i));
      lpRecon.push_back(hn(length-i-1));
      hpRecon.push_back(gn(length-i-1));
    }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/

void symmExtention(std::vector<Double> &sig, UInt paddLength)
{
  try
  {
    UInt len = sig.size();
    for(UInt i =0; i < paddLength; i++)
    {
      double temp1= sig.at(i * 2);
      double temp2= sig.at(len - 1);
      sig.insert(sig.begin(),temp1);
      sig.insert(sig.end(),temp2);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/****************************************************************/

void downSampling(std::vector<Double> &signal, UInt factor, std::vector<Double> &signalNew)
{
  try
  {
    UInt len = signal.size();
    UInt lenNew = (UInt) ceil( (double) len / (double) factor);
    for(UInt i = 0; i < lenNew; i++)
    {
      double temp = signal.at(i*factor);
      signalNew.push_back(temp);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void upSampling(std::vector<Double> &signal, UInt factor, std::vector<Double> &signalNew)
{
  try
  {
    UInt len     = signal.size();
    UInt lenNew  = (UInt) ceil( (Double) len * (Double) factor);
    for(UInt i = 0; i < lenNew; i++)
      if(i % factor == 0)
      {
        Double temp = signal.at(i/factor);
        signalNew.push_back(temp);
      }
    else
      signalNew.push_back(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

} // end namespace Wavelets

/***********************************************/
