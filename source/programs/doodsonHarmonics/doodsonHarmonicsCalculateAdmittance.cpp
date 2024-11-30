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
  void admitMatrix(UInt degreeInterpolation, UInt degreeExtrapolation, UInt degree, Int order,
                   const TideGeneratingPotential &tgpMajor, const std::vector<UInt> &majorIdx,
                   const TideGeneratingPotential &tgp, MatrixSliceRef A) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonicsCalculateAdmittance, SINGLEPROCESS, "computes the admittance function from TGP", DoodsonHarmonics)

/***********************************************/

void DoodsonHarmonicsCalculateAdmittance::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
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
    readConfig(config, "excludeDoodsonForInterpolation",   doodsonExclude,      Config::OPTIONAL, R"(["164.554", "164.555", "164.556", "164.566"])", "major tides not used for interpolation");
    if(isCreateSchema(config)) return;

    // =====================================================

    // read ocean tide file
    // --------------------
    logStatus<<"read doodson harmonics file <"<<fileNameDoodson<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameDoodson, d);
    if(!std::is_sorted(d.doodson.begin(), d.doodson.end()))
      throw(Exception("must be sorted"));

    logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
    TideGeneratingPotential tgp;
    readFileTideGeneratingPotential(fileNameTGP, tgp);

    // =====================================================

    auto isInList = [](const std::vector<Doodson> &list, const Doodson &doodson)
    {
      return std::find(list.begin(), list.end(), doodson) != list.end();
    };

    // eliminate small amplitudes (except all major tides)
    tgp.erase(std::remove_if(tgp.begin(), tgp.end(), [&](auto &t)
                             {return (t.amplitude() < threshold) && !isInList(d.doodson, t);}), tgp.end());

    // add missing tides with zero admittance
    for(auto &dood: d.doodson)
    {
      auto iter = std::lower_bound(tgp.begin(), tgp.end(), dood, [](const Doodson &tgp, const Doodson &dood){return tgp < dood;});
      if((iter == tgp.end()) || (dynamic_cast<const Doodson&>(*iter) != dood))
        tgp.insert(iter, TideGeneratingConstituent(dood, 0, 0., 0.));
    }

    // set tgp to zero for exluding doodson
    for(auto &t : tgp)
      if(isInList(doodsonExclude, t))
        t.c = t.s = 0.;

    TideGeneratingPotential tgpMajor;
    std::vector<UInt>       idxMajor;
    for(UInt i=0; i<tgp.size(); i++)
      if(isInList(d.doodson, tgp.at(i)))
      {
        tgpMajor.push_back(tgp.at(i));
        idxMajor.push_back(i);
      }
    if(tgpMajor.size() != d.doodson.size())
      throw(Exception("major tides not unique"));

    // =====================================================

    // interpolate each dgeree and order
    // ---------------------------------
    Matrix A(d.doodson.size(), tgp.size());
    logStatus<<"Interpolation  from "<<A.rows()<<" major tides to "<<A.columns()<<" tidal frequencies"<<Log::endl;
    for(UInt degree=2; degree<=3; degree++)
      for(UInt order=0; order<=degree; order++)
      {
        TideGeneratingPotential tgpMajorBand;
        std::vector<UInt>       idxMajorBand;
        for(UInt i=0; i<tgpMajor.size(); i++)
          if((tgpMajor.at(i).degree == degree) && (tgpMajor.at(i).d[0] == Int(order)) && tgpMajor.at(i).admit())
          {
            tgpMajorBand.push_back(tgpMajor.at(i));
            idxMajorBand.push_back(i);
          }

        if(tgpMajorBand.size())
          admitMatrix(degreeInterpolation, degreeExtrapolation, degree, order, tgpMajorBand, idxMajorBand, tgp, A);
      }

    // use major tides directly
    // ------------------------
    for(UInt k=0; k<d.doodson.size(); k++)
    {
      A.column(idxMajor.at(k)).setNull();
      A(k, idxMajor.at(k)) = 1.0;
    }

    // delete unused minor tides
    // -------------------------
    UInt countUnused = 0;
    for(UInt i=0; i<A.columns(); i++)
      if(isStrictlyZero(A.column(i)))
       countUnused++;
    Matrix Afull(A);
    A = Matrix(d.doodson.size(), Afull.columns()-countUnused);
    UInt idx = 0;
    for(UInt i=0; i<Afull.columns(); i++)
      if(isStrictlyZero(Afull.column(i)))
        tgp.erase(tgp.begin()+idx);
      else
        copy(Afull.column(i), A.column(idx++));

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

void DoodsonHarmonicsCalculateAdmittance::admitMatrix(UInt degreeInterpolation, UInt degreeExtrapolation, UInt degree, Int order,
                                                      const TideGeneratingPotential &tgpMajor, const std::vector<UInt> &majorIdx,
                                                      const TideGeneratingPotential &tgp, MatrixSliceRef A) const
{
  try
  {
    // extrapolation polynomial
    Matrix B;
    if(degreeExtrapolation != MAX_UINT)
    {
      degreeExtrapolation = std::min(degreeExtrapolation, tgpMajor.size()-1);
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
    degreeInterpolation = std::min(degreeInterpolation, tgpMajor.size()-1);
    std::vector<Matrix> P(tgpMajor.size()-degreeInterpolation);
    for(UInt i=0; i<P.size(); i++)
    {
      P.at(i) = Matrix(degreeInterpolation+1, degreeInterpolation+1);
      for(UInt k=0; k<=degreeInterpolation; k++)
        for(UInt n=0; n<=degreeInterpolation; n++)
          P.at(i)(k, n) = std::pow(tgpMajor.at(k+i).frequency(), n);
      inverse(P.at(i));
      for(UInt k=0; k<=degreeInterpolation; k++)
        P.at(i).column(k) *= 1./tgpMajor.at(k+i).admit();
    }

    UInt idx = 0;
    for(UInt i=0; i<tgp.size(); i++)
      if((tgp.at(i).degree == degree) && (tgp.at(i).d[0] == order) && tgp.at(i).admit())
      {
        const Double freq = tgp.at(i).frequency();
        if((idx < P.size()-1) && (freq > tgpMajor.at(idx+(degreeInterpolation+1)/2).frequency()))
          idx++;

        const Double dfreq = PI/27.3216; // [rad/day] half tidal group
        if(((tgpMajor.at(0).frequency()-dfreq <= freq) && (freq <= tgpMajor.back().frequency()+dfreq)) || (degreeExtrapolation == MAX_UINT))
        {
          // interpolation
          for(UInt n=0; n<=degreeInterpolation; n++)
            for(UInt k=0; k<=degreeInterpolation; k++)
              A(majorIdx.at(k+idx), i) += tgp.at(i).admit() * std::pow(freq, n) * P.at(idx)(n, k);
        }
        else
        {
          // extrapolation
          for(UInt n=0; n<B.rows(); n++)
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
