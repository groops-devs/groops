/***********************************************/
/**
* @file DigitalFilterNotch.cpp
*
* Implemented after Sophocles J. Orfanidis, "Introduction To Signal Processing", 2009.
*
* @brief Notch filter.
*
* @author Andreas Kvas
* @date 2017-02-07
*
*/
/***********************************************/

#include "base/import.h"
#include "base/fourier.h"
#include "classes/digitalFilter/digitalFilterNotch.h"
#include "files/fileMatrix.h"

/***********************************************/

DigitalFilterNotch::DigitalFilterNotch(Config& config)
{
  try
  {
    Double w0;   // normalized notch frequency
    Double bw;   // bandwidth at -3db

    readConfig(config, "notchFrequency",    w0,                Config::MUSTSET,  "",  "normalized notch frequency w_n = (f_n/f_nyq)");
    readConfig(config, "bandWidth",         bw,                Config::MUSTSET,  "",  "bandwidth at -3db. Quality factor of filter Q = w_n/bw");
    readConfig(config, "backwardDirection", backward,          Config::DEFAULT,  "0", "apply filter in backward direction");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    w0 *= PI;
    bw *= PI;

    Double gb = 1/std::sqrt(2); // |H(w)|^2 = 1/2, p. 575 11.3.1

    Double beta = (std::sqrt(1.0-gb*gb)/gb)*std::tan(bw/2.0); // p. 575 11.3.4
    Double gain = 1.0/(1.0+beta); // p. 575 11.3.6

    bn = Vector(3, gain);
    bn(1) = -2.0*std::cos(w0)*gain; // bn = gain*[1.0, -2*cos(w0), 1.0], p. 575 11.3.7
    bnStartIndex = 0;

    an = Vector(3, 1.0);
    an(1) = -2.0*gain*std::cos(w0);
    an(2) = 2.0*gain-1.0; // an = [1.0, -2*gain*cos(w0), 2*gain-1.0], p. 575 11.3.7
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
