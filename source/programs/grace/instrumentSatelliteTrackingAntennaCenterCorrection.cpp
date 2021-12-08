/***********************************************/
/**
* @file instrumentSatelliteTrackingAntennaCenterCorrection.cpp
*
* @brief Compute antenna center correction from orbit configuration.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the correction due to offset of the antenna center relative the center of mass.
The offsets $\M c_A$ and $\M c_B$ in \configFile{inputfileAntennaCenters}{matrix} are given in the satellite
reference frame. These offsets are rotated into the the inertial frame with $\M D_A$ and $\M D_B$ from
\configFile{inputfileStarCamera}{instrument} and projected onto the line of sight (LOS)
\begin{equation}
  \rho_{AOC} = \M e_{AB}\cdot(\M D_A\,\M c_A - \M D_B\,\M c_B),
\end{equation}
with the unit vector in line of sight direction
\begin{equation}
  \M e_{AB} = \frac{\M r_B - \M r_A}{\left\lVert{\M r_B - \M r_A}\right\rVert}.
\end{equation}

This program is also able to compute antenna offset correction covariance matrix,
which is the propagated variance of each satellite orientation observation to the AOC:
\begin{equation}
  \pmb{\Sigma}_{\Delta\rho_{AOC}}^{Sat} = \M F_{\pmb{\alpha}}^{Sat}\pmb{\Sigma}_{\pmb{\alpha\alpha}}^{Sat}(\M F_{\pmb{\alpha}}^{Sat})^T,
\end{equation}
where $\M F_{\pmb{\alpha}}^{Sat}$ is the partial derivative of the AOC w.r.t. the small angle rotations $\pmb{\alpha}$,
\begin{equation}
  \M F_{\pmb{\alpha}}^{Sat} =
  \begin{bmatrix}
  \frac{\partial \Delta\rho_{AOC}^{Sat}(t_1)}{\partial\pmb{\alpha}(t_1)} &        &  \\
          & \ddots &  \\
          &        & \frac{\partial \Delta\rho_{AOC}^{Sat}(t_N)}{\partial\pmb{\alpha}(t_N)}
  \end{bmatrix}.
\end{equation}
The contribution of indivudal satellite to AOC covariance matrix can be reported by
\configFile{outputfileSatelliteTrackingCovariance1/2}{matrix}. The complete AOC covariance matrix
\configFile{outputfileSatelliteTrackingCovariance}{matrix} considers the total influence of both satellites:
\begin{equation}
  \pmb{\Sigma}_{\Delta\rho_{AOC}} = \pmb{\Sigma}_{\Delta\rho_{AOC}}^{Sat1} + \pmb{\Sigma}_{\Delta\rho_{AOC}}^{Sat2}.
\end{equation}
The corrections for the range-rates and range-acceleration are computed by differentiating
an interpolation polynomial of degree \config{interpolationDegree}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "misc/grace/graceKBandGeometry.h"

/***** CLASS ***********************************/

/** @brief Compute antenna center correction from orbit configuration.
* @ingroup programsGroup */
class InstrumentSatelliteTrackingAntennaCenterCorrection
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentSatelliteTrackingAntennaCenterCorrection, PARALLEL, "compute antenna center correction from orbit configuration", Grace, Instrument)

/***********************************************/

void InstrumentSatelliteTrackingAntennaCenterCorrection::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName    outSSTName;
    FileName    outSSTCovarianceName, outSSTCovarianceName1, outSSTCovarianceName2;
    FileName    inCovarianceSca1, inCovarianceSca2;
    FileName    inCovarianceVarianceFactors;
    UInt        indexVarianceFactor1, indexVarianceFactor2;
    FileName    orbit1Name,  starCamera1Name;
    FileName    orbit2Name,  starCamera2Name;
    Vector3d    center1, center2;
    UInt        degree;
    UInt        sstType = 99;
    std::string choice;

    readConfig(config, "outputfileSatelliteTracking", outSSTName, Config::OPTIONAL, "", "corrections for range, range-rate, and range-accelerations");
    if(readConfigSequence(config, "variancePropagationStarCamera", Config::OPTIONAL, "", ""))
    {
      if(readConfigChoice(config, "outputFile", choice, Config::MUSTSET, "", ""))
      {
        if(readConfigChoiceElement(config, "combined",   choice, ""))
          readConfig(config, "outputfileSatelliteTrackingCovariance",  outSSTCovarianceName,  Config::MUSTSET, "grace_sigmaAoc.dat", "Covariance for chosen SST type");
        if(readConfigChoiceElement(config, "individual", choice, ""))
        {
          readConfig(config, "outputfileSatelliteTrackingCovariance1", outSSTCovarianceName1, Config::MUSTSET, "grace1_sigmaAoc.dat", "Covariance for chosen SST type");
          readConfig(config, "outputfileSatelliteTrackingCovariance2", outSSTCovarianceName2, Config::MUSTSET, "grace2_sigmaAoc.dat", "Covariance for chosen SST type");
        }
        endChoice(config);
      }
      readConfig(config, "inputFileCovarianceSca1", inCovarianceSca1, Config::MUSTSET, "", "Full covariance matrix from sensor fusion");
      readConfig(config, "inputFileCovarianceSca2", inCovarianceSca2, Config::MUSTSET, "", "Full covariance matrix from sensor fusion");
      if(readConfigSequence(config, "varianceFactors", Config::OPTIONAL, "", ""))
      {
        readConfig(config, "sigmasCovarianceMatrixArc", inCovarianceVarianceFactors, Config::MUSTSET, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
        readConfig(config, "indexVarianceFactor1",      indexVarianceFactor1,        Config::MUSTSET, "",  "in sigmasCovarianceMatrixArc");
        readConfig(config, "indexVarianceFactor2",      indexVarianceFactor2,        Config::MUSTSET, "",  "in sigmasCovarianceMatrixArc");
        endSequence(config);
      }
      if(readConfigChoice(config, "sstType", choice, Config::MUSTSET, "", ""))
      {
        if(readConfigChoiceElement(config, "range",             choice, "")) sstType = 0;
        if(readConfigChoiceElement(config, "rangeRate",         choice, "")) sstType = 1;
        if(readConfigChoiceElement(config, "rangeAcceleration", choice, "")) sstType = 2;
        endChoice(config);
      }
      endSequence(config);
    }

    readConfig(config, "inputfileOrbit1",             orbit1Name,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit2",             orbit2Name,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera1",        starCamera1Name, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera2",        starCamera2Name, Config::MUSTSET,  "", "");
    if(readConfigChoice(config, "antennaCenters",     choice,          Config::MUSTSET, "", "KBR antenna phase center"))
    {
      if(readConfigChoiceElement(config, "value", choice, ""))
      {
        readConfig(config, "center1X", center1.x(), Config::DEFAULT,   "1.4451172588", "x-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Y", center1.y(), Config::DEFAULT,  "-0.0004233040", "y-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Z", center1.z(), Config::DEFAULT,   "0.0022786600", "z-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center2X", center2.x(), Config::DEFAULT,   "1.4443870350", "x-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Y", center2.y(), Config::DEFAULT,   "0.0005761203", "y-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Z", center2.z(), Config::DEFAULT,   "0.0033040887", "z-coordinate of antenna position in SRF [m] for GRACEB");
      }
      if(readConfigChoiceElement(config, "file",  choice, ""))
      {
        FileName fileName;
        readConfig(config, "inputAntennaCenters", fileName, Config::MUSTSET, "", "");
        if(!isCreateSchema(config))
        {
          Matrix x;
          readFileMatrix(fileName, x);
          center1 = Vector3d(x(0,0), x(1,0), x(2,0));
          center2 = Vector3d(x(3,0), x(4,0), x(5,0));
        }
      }
      endChoice(config);
    }
    readConfig(config, "interpolationDegree", degree, Config::DEFAULT,  "2", "differentiation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    Double varianceFactor1 = 1.0;
    Double varianceFactor2 = 1.0;

    if(!inCovarianceVarianceFactors.empty())
    {
      logStatus<<"read variance factors"<<Log::endl;
      Vector sigmas;
      readFileMatrix(inCovarianceVarianceFactors, sigmas);
      varianceFactor1 = sigmas(indexVarianceFactor1);
      varianceFactor2 = sigmas(indexVarianceFactor2);
    }

    logStatus<<"read orbit and star camera data and generate antenna offset corrections"<<Log::endl;
    InstrumentFile  orbit1File(orbit1Name);
    InstrumentFile  orbit2File(orbit2Name);
    InstrumentFile  starCamera1File(starCamera1Name);
    InstrumentFile  starCamera2File(starCamera2Name);
    InstrumentFile::checkArcCount({orbit1File, orbit2File, starCamera1File, starCamera2File});

    std::vector<Arc> arcList(orbit1File.arcCount(), Epoch::SATELLITETRACKING);
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      // Variance propagation, if desired
      if(!outSSTCovarianceName.empty() || !outSSTCovarianceName1.empty())
      {
        Matrix CovarianceSca1, CovarianceSca2;
        readFileMatrix(inCovarianceSca1.appendBaseName(".arc"+arcNo%"%03i"s), CovarianceSca1);
        readFileMatrix(inCovarianceSca2.appendBaseName(".arc"+arcNo%"%03i"s), CovarianceSca2);

        fillSymmetric(CovarianceSca1); CovarianceSca1.setType(Matrix::GENERAL);
        fillSymmetric(CovarianceSca2); CovarianceSca2.setType(Matrix::GENERAL);

        Matrix partial1, partial2;
        GraceKBandGeometry::partialOfAntennaCenterCorrectionWrtRollPitchYaw(orbit1File.readArc(arcNo), orbit2File.readArc(arcNo),
                                                                            starCamera1File.readArc(arcNo), starCamera2File.readArc(arcNo),
                                                                            center1, center2,
                                                                            partial1, partial2);

        const Double dt = medianSampling(orbit1File.readArc(arcNo).times()).seconds();
        Matrix Covariance1 = std::pow(varianceFactor1,2) * GraceKBandGeometry::variancePropagationStarCamera2SatelliteRanging(partial1, CovarianceSca1, dt, degree, sstType);
        Matrix Covariance2 = std::pow(varianceFactor2,2) * GraceKBandGeometry::variancePropagationStarCamera2SatelliteRanging(partial2, CovarianceSca2, dt, degree, sstType);

        if(!outSSTCovarianceName.empty())
          writeFileMatrix(outSSTCovarianceName.appendBaseName(".arc"+arcNo%"%03i"s),  Covariance1 + Covariance2);
        if(!outSSTCovarianceName1.empty())
          writeFileMatrix(outSSTCovarianceName1.appendBaseName(".arc"+arcNo%"%03i"s), Covariance1);
        if(!outSSTCovarianceName2.empty())
          writeFileMatrix(outSSTCovarianceName2.appendBaseName(".arc"+arcNo%"%03i"s), Covariance2);
      }

      // antenna center correction
      return GraceKBandGeometry::antennaCenterCorrection(orbit1File.readArc(arcNo), orbit2File.readArc(arcNo),
                                                         starCamera1File.readArc(arcNo), starCamera2File.readArc(arcNo),
                                                         center1, center2, degree);
    }, comm); // forEach

    if(Parallel::isMaster(comm) && !outSSTName.empty())
    {
      logStatus<<"write tracking data to file <"<<outSSTName<<">"<<Log::endl;
      InstrumentFile::write(outSSTName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
