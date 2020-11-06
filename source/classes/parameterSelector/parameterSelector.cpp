/***********************************************/
/**
* @file parameterSelector.cpp
*
* @brief Index vector from selected parameters.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#define DOCSTRING_ParameterSelector

#include "base/import.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/parameterSelector/parameterSelectorComplement.h"
#include "classes/parameterSelector/parameterSelectorMatrix.h"
#include "classes/parameterSelector/parameterSelectorNames.h"
#include "classes/parameterSelector/parameterSelectorRange.h"
#include "classes/parameterSelector/parameterSelectorWildcard.h"
#include "classes/parameterSelector/parameterSelectorZeros.h"
#include "classes/parameterSelector/parameterSelector.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParameterSelector, "parameterSelectorType",
                      ParameterSelectorWildcard,
                      ParameterSelectorNames,
                      ParameterSelectorRange,
                      ParameterSelectorMatrix,
                      ParameterSelectorZeros,
                      ParameterSelectorComplement)

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParameterSelector, "parameterSelectorType")

/***********************************************/

ParameterSelector::ParameterSelector(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "selected parameters"))
    {
      renameDeprecatedChoice(config, type, "name", "wildcard", date2time(2020, 5, 30));

      if(readConfigChoiceElement(config, "wildcard",   type, "parameter name matching"))
        parameters.push_back(new ParameterSelectorWildcard(config));
      if(readConfigChoiceElement(config, "names",      type, "manual list of parameter names"))
        parameters.push_back(new ParameterSelectorNames(config));
      if(readConfigChoiceElement(config, "range",      type, "range of parameters"))
        parameters.push_back(new ParameterSelectorRange(config));
      if(readConfigChoiceElement(config, "matrix",     type, "matrix containing parameter indexes"))
        parameters.push_back(new ParameterSelectorMatrix(config));
      if(readConfigChoiceElement(config, "zeros",      type, "additional zero parameters"))
        parameters.push_back(new ParameterSelectorZeros(config));
      if(readConfigChoiceElement(config, "complement", type, "all parameters except those selected"))
        parameters.push_back(new ParameterSelectorComplement(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }

    varList = config.getVarList();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParameterSelector::~ParameterSelector()
{
  for(UInt i=0; i<parameters.size(); i++)
    delete parameters.at(i);
}

/***********************************************/

std::vector<UInt> ParameterSelector::indexVector(const std::vector<ParameterName> &parameterNames)
{
  try
  {
    std::vector<UInt> vector;
    for(const auto &param : parameters)
    {
      std::vector<UInt> vec = param->indexVector(parameterNames, varList);
      vector.insert(vector.end(), vec.begin(), vec.end());
    }
    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<UInt> ParameterSelector::indexVectorComplement(std::vector<UInt> vector, UInt referenceLength)
{
  try
  {
    std::vector<UInt> complement;
    std::vector<UInt> all(referenceLength);
    std::iota(all.begin(), all.end(), 0);

    if(!vector.size())
      return all;

    std::sort(vector.begin(), vector.end());
    std::set_difference(all.begin(), all.end(), vector.begin(), vector.end(), std::inserter(complement, complement.end()));

    return complement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
