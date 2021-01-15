/***********************************************/
/**
* @file normalsRegularizationBorders.cpp
*
* @brief Two regularization matrices for inside und outside of Border.
* The diagonal is saved as Vector.
*
* @author Annette Eicker
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program sets up two regularization matrices for two different regional areas.
For a given set of points defined by \configClass{grid}{gridType} it is evaluated, whether each point
(corresponding to an unknown parameter of a respective parameterization by space localizing basis functions)
is inside or outside a certain area given by \configClass{border}{borderType}.
Each regularization matrix is a diagonal matrix, one of them features a one if the
point is inside, and a zero if the point lies outside the area. The other matrix features
a zero if the point is inside, and a one if the point lies outside the area
This results in two regularization matrices with
\begin{equation}
\M R_1+\M R_2=\M I.
\end{equation}
The two matrices are provided as vectors of the diagonal
in the output files \configFile{outputfileOutside}{matrix} and \configFile{outputfileInside}{matrix}.
The regularization matrices are then used by \configClass{normalEquation:regularization}{normalEquationType:regularization}.
As an example, the two different areas could be oceanic regions on the one hand and continental areas on the other hand.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/grid/grid.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Two regularization matrices for inside und outside of Border.
* The diagonal is saved as Vector.
* @ingroup programsGroup */
class NormalsRegularizationBorders
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsRegularizationBorders, SINGLEPROCESS, "Two regularization matrices for inside und outside of Border", NormalEquation)

/***********************************************/

void NormalsRegularizationBorders::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName  insideName, outsideName;
    GridPtr   grid;
    BorderPtr border;

    readConfig(config, "outputfileInside",  insideName,  Config::MUSTSET,  "", "");
    readConfig(config, "outputfileOutside", outsideName, Config::MUSTSET,  "", "");
    readConfig(config, "grid",              grid,        Config::MUSTSET, "", "nodal point distribution of parameters, e.g harmonics splines");
    readConfig(config, "border",            border,      Config::MUSTSET, "", "regularization areas, e.g land and ocean");
    if(isCreateSchema(config)) return;


    std::vector<Vector3d> points = grid->points();
    logInfo<<"  nodal points: "<<points.size()<<Log::endl;

    Vector  diag1(points.size()); // fill with 1 or 0
    Vector  diag2(points.size()); // fill with 1 or 0

    for(UInt k=0; k<points.size(); k++)
    {
      if(border->isInnerPoint(points.at(k)))
        diag1(k)=1;
      else
        diag2(k)=1;
    }

    logStatus<<"write diagonal matrices to file"<<Log::endl;
    if(!insideName.empty())
      writeFileMatrix(insideName, diag1);
    if(!outsideName.empty())
      writeFileMatrix(outsideName, diag2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
