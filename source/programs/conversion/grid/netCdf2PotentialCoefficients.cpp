/***********************************************/
/**
* @file netCdf2PotentialCoefficients.cpp
*
* @brief Convert a NetCDF file to a sequence of PotentialCoefficients files
*
* @author Andreas Kvas
* @date 2019-10-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts a COARDS compliant NetCDF file into potential coefficients by least squares.
If multiple \config{variableNameData} are given the grids values are accumulated before the adjustment.

See also \program{NetCdfInfo}, \program{NetCdf2GridRectangular}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedData.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Convert COARDS compliant grids to PotentialCoefficients
* @ingroup programsConversionGroup */
class NetCdf2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NetCdf2PotentialCoefficients, PARALLEL, "Convert a NetCDF file to a sequence of PotentialCoefficients files", Conversion, Grid)

/***********************************************/

void NetCdf2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName    outName, inName;
    std::string lonName, latName, timeName, loopVar;
    std::vector<std::string> dataName;
    Double      a, f, GM, R;
    KernelPtr   kernel;
    UInt        minDegree, maxDegree;
    Double      noData = NAN_EXPR;

    readConfig(config, "outputfilePotentialCoefficients", outName,   Config::MUSTSET,  "", "One file for each epoch in the NetCDF file. Use loopTimeVariable as template.");
    readConfig(config, "loopTimeVariable",                loopVar,   Config::MUSTSET,  "loopTime", "");
    readConfig(config, "inputfileNetCdf",                 inName,    Config::MUSTSET,  "", "");
    readConfig(config, "variableNameLongitude",           lonName,   Config::MUSTSET,  "lon",  "name of NetCDF variable");
    readConfig(config, "variableNameLatitude",            latName,   Config::MUSTSET,  "lat",  "name of NetCDF variable");
    readConfig(config, "variableNameTime",                timeName,  Config::OPTIONAL, "time", "name of NetCDF variable (leave blank for static grids)");
    readConfig(config, "variableNameData",                dataName,  Config::OPTIONAL, "", "name of NetCDF variable");
    readConfig(config, "noDataValue",                     noData,    Config::OPTIONAL, "", "no data value for data variables");
    readConfig(config, "semimajorAxis",                   a,         Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference semimajor axis for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",               f,         Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    readConfig(config, "kernel",                          kernel,    Config::MUSTSET,  "",  "kernel in which the grid values are given");
    readConfig(config, "minDegree",                       minDegree, Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",                       maxDegree, Config::MUSTSET,  "",  "");
    readConfig(config, "GM",                              GM,        Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,         Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius for potential coefficients");
    if(isCreateSchema(config)) return;

#ifdef NOLIB_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else
    // set up grid
    // -----------
    GriddedDataRectangular grid;
    std::vector<Time> epochs(1);
    Matrix l;
    if(Parallel::isMaster(comm))
    {
      // open netCDF file
      // ----------------
      logStatus<<"read file <"<<inName<<">"<<Log::endl;
      NetCdf::InFile file(inName);
      NetCdf::Variable  lon    = file.variable(lonName);
      NetCdf::Variable  lat    = file.variable(latName);
      NetCdf::Dimension dimLon = lon.dimensions().at(0);
      NetCdf::Dimension dimLat = lat.dimensions().at(0);

      // set up grid
      // -----------
      grid.ellipsoid  = Ellipsoid(a, f);
      grid.longitudes = NetCdf::convertAngles(lon.values());
      grid.latitudes  = NetCdf::convertAngles(lat.values());
      grid.heights.resize(grid.latitudes.size(), 0.0);

      // set up time axis
      // ----------------
      NetCdf::Dimension dimTime;
      if(!timeName.empty())
      {
        auto var = file.variable(timeName);
        dimTime  = var.dimensions().at(0);
        epochs   = NetCdf::convertTimes(var.values(), var.attribute("units").value());
      }

      l = Matrix(grid.longitudes.size()*grid.latitudes.size(), epochs.size());

      // data variables
      // --------------
      Single::forEach(epochs.size(), [&](UInt idEpoch)
      {
        for(const std::string &name : dataName)
        {
          auto var  = file.variable(name);
          auto dims = var.dimensions();

          std::vector<UInt> start(dims.size(), 0);
          std::vector<UInt> count;
          for(auto &dim : dims)
            count.push_back(dim.length());

          if((dims.size() != 2) && (timeName.empty() || (dims.size() != 3)))
            throw(Exception("variable <"+name+"> has wrong dimensions"));

          if(!timeName.empty() && (dims.size() > 2))
          {
            if(dims.at(0) != dimTime)
              throw(Exception("variable <"+name+"> must have time as first dimension"));
            start.at(0) = idEpoch;
            count.at(0) = 1;
          }

          Vector values = var.values(start, count);
          for(UInt k=0; k<values.rows(); k++)
            if(std::isnan(values(k)) || (values(k) == noData))
              values(k) = 0.0;

          if((dims.at(dims.size()-1) == dimLat) && (dims.at(dims.size()-2) == dimLon))
            axpy(1., flatten(reshape(values, count.at(dims.size()-1), count.at(dims.size()-2)).trans()), l.column(idEpoch));
          else if((dims.at(dims.size()-1) == dimLon) && (dims.at(dims.size()-2) == dimLat))
            axpy(1., values, l.column(idEpoch));
          else
            throw(Exception("variable <"+name+"> must have ("+latName+", "+lonName+") dimensions"));
        }
      });
    }
    logStatus<<"distribute grid values"<<Log::endl;
    Parallel::broadCast(grid, 0, comm);
    Parallel::broadCast(l, 0, comm);

    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius, dLambda, dPhi;
    grid.geocentric(lambda, phi, radius, dLambda, dPhi);

    // set up normal equations
    // -----------------------
    logStatus<<"compute potential coefficients order-by-order"<<Log::endl;
    Matrix cosml(lambda.size(), maxDegree+1);
    Matrix sinml(lambda.size(), maxDegree+1);
    for(UInt m=0; m <maxDegree+1; m++)
      for(UInt j=0; j<lambda.size(); j++)
      {
        cosml(j, m) = std::cos(m*static_cast<Double>(lambda.at(j)));
        sinml(j, m) = std::sin(m*static_cast<Double>(lambda.at(j)));
      }

    std::vector<Matrix> N, n;
    N.push_back(Matrix(maxDegree+1, Matrix::SYMMETRIC));
    n.push_back(Matrix(maxDegree+1, l.columns()));
    for(UInt m=1; m<maxDegree+1; m++)
    {
      N.push_back(Matrix(2*(maxDegree+1-m), Matrix::SYMMETRIC));
      n.push_back(Matrix(2*(maxDegree+1-m), l.columns()));
    }

    logStatus<<"accumulate normal equations"<<Log::endl;
    Parallel::forEach(phi.size(), [&](UInt i)
    {
      Matrix Pnm = SphericalHarmonics::Pnm(Angle(PI*0.5 - phi.at(i)), radius.at(i)/R, maxDegree);
      Vector kn  = kernel->inverseCoefficients(polar(Angle(0.0), phi.at(i), radius.at(i)), maxDegree);
      for(UInt n=0; n<maxDegree+1; n++)
        Pnm.slice(n,0,1,n+1) *= GM/R * kn(n);

      for(UInt m=0; m<maxDegree+1; m++)
      {
        Matrix A;
        if(m == 0)
        {
          A = cosml.column(0)*Pnm.column(0).trans();
        }
        else
        {
          A = Matrix(lambda.size(), 2*(maxDegree+1-m));
          matMult(1.0, cosml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(0, maxDegree+1-m));
          matMult(1.0, sinml.column(m), Pnm.slice(m, m, maxDegree+1-m, 1).trans(), A.column(maxDegree+1-m, maxDegree+1-m));
        }
        const Double weight = dPhi.at(i)*std::cos(phi.at(i))/2.;
        rankKUpdate(weight, A, N.at(m));
        matMult(weight, A.trans(), l.row(i*lambda.size(), lambda.size()), n.at(m));
      }
    }, comm);

    logStatus<<"solve normal equations"<<Log::endl;
    std::vector<Matrix> x;
    for(UInt m=0; m<maxDegree+1; m++)
    {
      Parallel::reduceSum(N.at(m), 0, comm);
      Parallel::reduceSum(n.at(m), 0, comm);
      if(Parallel::isMaster(comm))
      {
        for(UInt k=0; k<N.at(m).rows(); k++)
          if(N.at(m)(k, k) == 0.0)
            N.at(m)(k, k) = 1.0;
        x.push_back(solve(N.at(m), n.at(m)));
      }
    }

    logStatus<<"write results to <"<<outName<<">"<<Log::endl;
    if(Parallel::isMaster(comm))
    {
      VariableList fileNameVariableList;
      addVariable(loopVar, fileNameVariableList);

      Matrix Cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix Snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      for(UInt k=0; k<epochs.size(); k++)
      {
        copy(x.at(0).column(k), Cnm.column(0));
        for(UInt m=1; m<maxDegree+1; m++)
        {
          copy(x.at(m).slice(0,             k, maxDegree+1-m, 1), Cnm.slice(m, m, maxDegree+1-m, 1));
          copy(x.at(m).slice(maxDegree+1-m, k, maxDegree+1-m, 1), Snm.slice(m, m, maxDegree+1-m, 1));
        }

        fileNameVariableList[loopVar]->setValue(epochs.at(k).mjd());
        writeFileSphericalHarmonics(outName(fileNameVariableList), SphericalHarmonics(GM, R, Cnm, Snm).get(maxDegree, minDegree));
      }
    }
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
