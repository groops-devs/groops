/***********************************************/
/**
* @file doodsonHarmonicsCalculateAdmittance.cpp
*
* @brief Computes the admittance function from TGP.
*
* @author Torsten Mayer-Guerr
* @date 2009-04-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes the admittance function to interpolate minor tides from
tides given in \configFile{inputfileDoodsonHarmonics}{doodsonHarmonic}
using \configFile{inputfileTideGeneratingPotential}{tideGeneratingPotential}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "base/tideGeneratingPotential.h"
#include "files/fileAdmittance.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileTideGeneratingPotential.h"

/***** CLASS ***********************************/

/** @brief Computes the admittance function from TGP.
* @ingroup programsGroup */
class DoodsonHarmonicsCalculateAdmittance
{
  void admitMatrix(UInt degreeInterpolation, UInt degreeExtrapolation, Int idBand,
                   const TideGeneratingPotential &tgpMajor, const std::vector<UInt> &majorIdx,
                   const TideGeneratingPotential &tgp, MatrixSliceRef A) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonicsCalculateAdmittance, SINGLEPROCESS, "computes the admittance function from TGP", DoodsonHarmonics)

/***********************************************/

void DoodsonHarmonicsCalculateAdmittance::run(Config &config)
{
  try
  {
    FileName fileNameAdmit, fileNameDoodson, fileNameTGP;
    Double   threshold;
    UInt     degreeInterpolation, degreeExtrapolation = MAX_UINT;
    std::vector<Doodson> doodsonExclude;

    readConfig(config, "outputfileAdmittance",             fileNameAdmit,       Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileDoodsonHarmonics",        fileNameDoodson,     Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP,         Config::OPTIONAL, "{groopsDataDir}/tides/generatingTide_HW95.txt", "TGP");
    readConfig(config, "threshold",                        threshold,           Config::DEFAULT,  "1e-4", "[m^2/s^2] only interpolate tides with TGP greater than threshold");
    readConfig(config, "degreeInterpolation",              degreeInterpolation, Config::DEFAULT,  "1",    "polynomial degree for interpolation");
    readConfig(config, "degreeExtrapolation",              degreeExtrapolation, Config::OPTIONAL, "1",    "polynomial degree for extrapolation");
    readConfig(config, "excludeDoodsonForInterpolation",   doodsonExclude,      Config::OPTIONAL, R"(["S1", "S2"])", "major tides not used for interpolation");
    if(isCreateSchema(config)) return;

    // =====================================================

    // read ocean tide file
    // --------------------
    logStatus<<"read doodson harmonics file <"<<fileNameDoodson<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameDoodson, d);

    logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
    TideGeneratingPotential tgp;
    readFileTideGeneratingPotential(fileNameTGP, tgp);

    // =====================================================

    // eliminate small amplitudes (except all major tides)
    tgp.erase(std::remove_if(tgp.begin(), tgp.end(), [&](const TideGeneratingConstituent &x)
                             {
                               if(std::find(d.doodson.begin(), d.doodson.end(), x) != d.doodson.end())
                                 return FALSE; // major tides must be included
                               return (x.amplitude() < threshold);
                             }), tgp.end());

    // add missing tides with zero admittance
    for(auto &dood: d.doodson)
    {
      auto iter = std::lower_bound(tgp.begin(), tgp.end(), dood);
      if((iter == tgp.end()) || (*iter != dood))
        tgp.insert(iter, TideGeneratingConstituent(dood, 0., 0.));
    }

    // set tgp to zero for exluding doodson
    for(auto &dood: doodsonExclude)
    {
      auto iter = std::lower_bound(tgp.begin(), tgp.end(), dood);
      if((iter != tgp.end()) && (*iter == dood))
        iter->c = iter->s = 0.;
    }

    TideGeneratingPotential tgpMajor;
    std::copy_if(tgp.begin(), tgp.end(), std::back_inserter(tgpMajor), [&](const auto &x) {return std::find(d.doodson.begin(), d.doodson.end(), x) != d.doodson.end();});

    // =====================================================

    // interpolate band by band (long, diurnal, semi diurnal)
    // ------------------------------------------------------
    Matrix A(d.doodson.size(), tgp.size());
    logStatus<<"Interpolation  from "<<A.rows()<<" major tides to "<<A.columns()<<" tidal frequencies"<<Log::endl;
    for(Int idBand=0; idBand<3; idBand++)
    {
      TideGeneratingPotential tgpMajorBand;
      std::vector<UInt>       idxMajorBand;
      for(UInt i=0; i<tgpMajor.size(); i++)
        if((tgpMajor.at(i).d[0] == idBand) && tgpMajor.at(i).admit())
        {
          tgpMajorBand.push_back(tgpMajor.at(i));
          idxMajorBand.push_back(i);
        }

      if(tgpMajorBand.size())
        admitMatrix(degreeInterpolation, degreeExtrapolation, idBand, tgpMajorBand, idxMajorBand, tgp, A);
    }

    // use major tides directly
    // ------------------------
    for(UInt i=0; i<d.doodson.size(); i++)
    {
      const UInt k = std::distance(tgp.begin(), std::lower_bound(tgp.begin(), tgp.end(), d.doodson.at(i)));
      A.column(k).setNull();
      A(i, k) = 1.0;
    }

    // =====================================================

    // write results
    // -------------
    logStatus<<"writing admittance file <"<<fileNameAdmit<<">"<<Log::endl;
    Admittance admit;
    admit.doodsonMajor.insert(admit.doodsonMajor.begin(), tgpMajor.begin(), tgpMajor.end());
    admit.doodsonMinor.insert(admit.doodsonMinor.begin(), tgp.begin(), tgp.end());
    admit.admittance = A;
    writeFileAdmittance(fileNameAdmit, admit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

void DoodsonHarmonicsCalculateAdmittance::admitMatrix(UInt degreeInterpolation, UInt degreeExtrapolation, Int idBand,
                                                      const TideGeneratingPotential &tgpMajor, const std::vector<UInt> &majorIdx,
                                                      const TideGeneratingPotential &tgp, MatrixSliceRef A) const
{
  try
  {
    // extrapolation polynomial
    Matrix B;
    if(degreeExtrapolation != MAX_UINT)
    {
      B = Matrix(tgpMajor.size(), degreeExtrapolation+1);
      for(UInt i=0; i<tgpMajor.size(); i++)
        for(UInt n=0; n<=degreeExtrapolation; n++)
          B(i, n) = std::pow(tgpMajor.at(i).frequency(), n);
      Matrix L(tgpMajor.size(), tgpMajor.size());
      for(UInt i=0; i<L.columns(); i++)
        L(i, i) = 1./tgpMajor.at(i).admit();
      B = leastSquares(B, L);
    }

    // interpolation polynomial
    std::vector<Matrix> P(tgpMajor.size()-degreeInterpolation);
    UInt idx = 0;
    for(UInt i=0; i<P.size(); i++)
    {
      P.at(i) = Matrix(degreeInterpolation+1, degreeInterpolation+1);
      for(UInt k=0; k<=degreeInterpolation; k++)
        for(UInt n=0; n<=degreeInterpolation; n++)
          P.at(i)(k, n) =  std::pow(tgpMajor.at(k+idx).frequency(), n);
      inverse(P.at(i));
      for(UInt n=0; n<=degreeInterpolation; n++)
        for(UInt k=0; k<=degreeInterpolation; k++)
          P.at(i)(n, k) /= tgpMajor.at(k+idx).admit();
      idx++;
    }

    idx = 0;
    for(UInt i=0; i<tgp.size(); i++)
      if((tgp.at(i).d[0] == idBand) && tgp.at(i).admit())
      {
        const Double freq = tgp.at(i).frequency();
        if((idx < P.size()-1) && (freq > tgpMajor.at(idx+degreeInterpolation).frequency()))
          idx++;

        if(((tgpMajor.at(0).frequency() <= freq) && (freq <= tgpMajor.back().frequency())) || (degreeExtrapolation == MAX_UINT))
        {
          // interpolation
          for(UInt n=0; n<=degreeInterpolation; n++)
            for(UInt k=0; k<=degreeInterpolation; k++)
              A(majorIdx.at(k+idx), i) += tgp.at(i).admit() * std::pow(freq, n) * P.at(idx)(n, k);
        }
        else
        {
          // extrapolation
          for(UInt n=0; n<=degreeExtrapolation; n++)
            for(UInt k=0; k<tgpMajor.size(); k++)
              A(majorIdx.at(k), i) += tgp.at(i).admit() * std::pow(freq, n) * B(n, k);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
