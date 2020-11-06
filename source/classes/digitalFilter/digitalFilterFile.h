/***********************************************/
/**
* @file digitalFilterFile.h
*
* @brief Filter from file.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERFILE__
#define __GROOPS_DIGITALFILTERFILE__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterFile = R"(
\subsection{File}
Read filter coefficients of \eqref{digitalFilterType:arma} from a coefficient file.
One column might define the index $n$
of the coefficients $a_n$ and $b_n$ in the other columns.
)";
#endif

/***********************************************/

#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Filter from file.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterFile : public DigitalFilterARMA
{
public:
  DigitalFilterFile(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterFile::DigitalFilterFile(Config &config)
{
  try
  {
    FileName fileName;
    ExpressionVariablePtr exprId, exprMA, exprAR;

    readConfig(config, "inputfileMatrix",   fileName,          Config::MUSTSET,   "",      "matrix with filter coefficients");
    readConfig(config, "index",             exprId,            Config::OPTIONAL, "data0", "index of coefficients (input columns are named data0, data1, ...)");
    readConfig(config, "bn",                exprMA,            Config::OPTIONAL, "data1", "MA coefficients (moving average) (input columns are named data0, data1, ...)");
    readConfig(config, "an",                exprAR,            Config::OPTIONAL, "data2", "AR coefficients (autoregressive) (input columns are named data0, data1, ...)");
    readConfig(config, "backwardDirection", backward,          Config::DEFAULT,   "0",     "apply filter in backward direction");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,   "0",     "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,   "",      "");
    if(isCreateSchema(config)) return;

    Matrix A;
    readFileMatrix(fileName, A);

    // create data variables
    // ---------------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    if(exprId) exprId->usedVariables(varList, usedVariables);
    if(exprMA) exprMA->usedVariables(varList, usedVariables);
    if(exprAR) exprAR->usedVariables(varList, usedVariables);
    addDataVariables(A, varList, usedVariables);

    std::vector< std::pair<Int, Double> > bk;
    std::vector< std::pair<Int, Double> > ak;
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      Int idx = exprId ? static_cast<Int>(exprId->evaluate(varList)) : static_cast<Int>(i);
      if(exprMA && exprMA->evaluate(varList) != 0.0) bk.push_back({idx, exprMA->evaluate(varList)});
      if(exprAR && exprAR->evaluate(varList) != 0.0) ak.push_back({idx, exprAR->evaluate(varList)});
    }
    if(std::any_of(ak.begin(), ak.end(), [](const std::pair<Int, Double> &a) { return a.first<0; }))
      throw(Exception("negative indicies are not allowed for AR coefficients (non-causal filter)"));

    // allocate filter coefficient vectors
    // -----------------------------------
    auto indexComparable = [](const std::pair<Int, Double> &b1, const std::pair<Int, Double> &b2) { return b1.first < b2.first; };
    if(bk.size()>0)
    {
      auto minMaxIndex = std::minmax_element(bk.begin(), bk.end(), indexComparable);
      bnStartIndex = std::max(-(*minMaxIndex.first).first, 0);
      bn = Vector((*minMaxIndex.second).first + bnStartIndex + 1);
      for(auto &b : bk) bn(b.first+bnStartIndex) = b.second;
    }
    if(ak.size()>0)
    {
      auto maxIndex = std::max_element(ak.begin(), ak.end(), indexComparable);
      an = Vector((*maxIndex).first + 1);
      for(auto &a : ak) an(a.first) = a.second;
    }
    else
      an = Vector(1, 1.0);

    if(an(0) != 1.0)
    {
      logWarning<<"a0 coefficient set to one"<<Log::endl;
      an(0) = 1.0;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
