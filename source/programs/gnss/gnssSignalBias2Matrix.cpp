/***********************************************/
/**
* @file gnssSignalBias2Matrix.cpp
*
* @brief Computes signal biases for a given type list.
*
* @author Torsten Mayer-Guerr
* @date 2019-03-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes signal biases for a given list of \configClass{types}{gnssType}.
If the type list is empty, all types contained in \configFile{inputfileSignalBias}{gnssSignalBias} are used.
The resulting \configFile{outputfileMatrix}{matrix} contains a vector with an entry for each type.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Computes signal biases for a given type list.
* @ingroup programsGroup */
class GnssSignalBias2Matrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssSignalBias2Matrix, SINGLEPROCESS, "Computes signal biases for a given type list", Gnss)


/***********************************************/

void GnssSignalBias2Matrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameSignalBias, fileNameMatrix, fileNameTypes;
    std::vector<GnssType> types;

    readConfig(config, "outputfileMatrix",    fileNameMatrix,     Config::MUSTSET,  "", "");
    readConfig(config, "outputfileTypes",     fileNameTypes,      Config::OPTIONAL, "", "ASCII list of types");
    readConfig(config, "inputfileSignalBias", fileNameSignalBias, Config::MUSTSET,  "", "");
    readConfig(config, "types",               types,              Config::OPTIONAL, "", "if not set, all types in the file are used");
    if(isCreateSchema(config)) return;

    logStatus<<"read signal bias file <"<<fileNameSignalBias<<">"<<Log::endl;
    GnssSignalBias signalBias;
    readFileGnssSignalBias(fileNameSignalBias, signalBias);

    if(types.size() == 0)
      types = signalBias.types;
    Vector bias = signalBias.compute(types);

    if(!fileNameTypes.empty())
    {
      logStatus<<"write types to file <"<<fileNameTypes<<">"<<Log::endl;
      std::vector<std::string> typeNames;
      for(const auto &type : types)
        typeNames.push_back(type.str());
      writeFileStringList(fileNameTypes, typeNames);
    }

    logStatus<<"write matrix to file <"<<fileNameMatrix<<">"<<Log::endl;
    writeFileMatrix(fileNameMatrix, bias);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
