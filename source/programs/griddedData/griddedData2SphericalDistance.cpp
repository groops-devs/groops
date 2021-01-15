/***********************************************/
/**
* @file griddedData2SphericalDistance.cpp
*
* @brief Spherical distance between all point pairs of two grids.
*
* @author Andreas KVas
* @date 2018-01-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute the spherical distance on the unit sphere in radians between all point pairs of two grids.
The spherical distance is computed by
\begin{equation}
  \psi_{12} = \arccos(\M n_1 \cdot \M n_2),
\end{equation}
where $\M n_i$ is the (normalized) position. This implies that all points are projected onto the unit sphere.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Spherical distance between all point pairs of two grids.
* @ingroup programsGroup */
class GriddedData2SphericalDistance
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2SphericalDistance, SINGLEPROCESS, "spherical distance between all point pairs of two grids.", Grid, Matrix)

/***********************************************/

void GriddedData2SphericalDistance::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameDistance;
    GridPtr  gridPtr1, gridPtr2;

    readConfig(config, "outputfileMatrix", fileNameDistance, Config::MUSTSET, "", "matrix containing the spherical distance between all point pairs [rad]");
    readConfig(config, "grid1",             gridPtr1,         Config::MUSTSET,  "", "");
    readConfig(config, "grid2",             gridPtr2,         Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    std::vector<Vector3d> points1 = gridPtr1->points();
    std::vector<Vector3d> points2 = gridPtr2->points();
    std::for_each(points1.begin(), points1.end(), [](Vector3d &p) {p.normalize();});
    std::for_each(points2.begin(), points2.end(), [](Vector3d &p) {p.normalize();});

    logStatus<<"compute distance"<<Log::endl;
    Matrix S(points1.size(), points2.size());
    for(UInt i=0; i<points1.size(); i++)
      for(UInt k=0; k<points2.size(); k++)
        S(i,k) = std::acos(inner(points1.at(i), points2.at(k)));

    logStatus<<"writing matrix to file <"<<fileNameDistance<<">"<<Log::endl;
    writeFileMatrix(fileNameDistance, S);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
