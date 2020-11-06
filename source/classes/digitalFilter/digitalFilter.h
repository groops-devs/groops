/***********************************************/
/**
* @file digitalFilter.h
*
* @brief Digital filter implementation.
*
* @author Matthias Ellmer
* @author Andreas Kvas
* @date 2015-10-29
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTER__
#define __GROOPS_DIGITALFILTER__


// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilter = R"(
\section{DigitalFilter}\label{digitalFilterType}
Digital filter implementation for the filtering of equally spaced time series. This class implements the filter equations as
\begin{equation}\label{digitalFilterType:arma}
  \sum_{l=0}^Q a_l y_{n-l} = \sum_{k=-p_0}^{P-p_0-1} b_k x_{n-k}, \hspace{25pt} a_0 = 1,
\end{equation}
where $Q$ is the autoregressive (AR) order and $P$ is the moving average (MA) order. Note that the MA part can also be non-causal.
The characteristics of a filter cascade can be computed by the programs \program{DigitalFilter2FrequencyResponse} and \program{DigitalFilter2ImpulseResponse}.
To apply a filter cascade to a time series (or an instrument file ) use \program{InstrumentFilter}.
Each filter can be applyed in forward and backward direction by setting \config{backwardDirection}.
If the same filter is applied in both directions, the combined filter has zero phase and the squared magnitude response.
Setting \config{inFrequencyDomain} to true applies the transfer function of the filter to the DFT of the input and synthesizes the result, i.e.:
\begin{equation}
  y_n = \mathcal{F}^{-1}\{H\cdot\mathcal{F}\{x_n\}\}.
\end{equation}
This is equivalent to setting \config{padType} to \config{periodic}.

To reduce warmup effects, the input time series can be padded by chosing a \config{padType}:
\begin{itemize}
\item \config{none}: no padding is applied
\item \config{zero}: zeros are appended at the beginning and end of the input time series
\item \config{constant}: the beginning of the input time series is padded with the first value, the end is padded with the last value
\item \config{periodic}: periodic continuation of the input time series (i.,e. the beginning is padded with the last epochs and the end is padded with the first epochs)
\item \config{symmetric}: beginning and end are reflected around the first and last epoch respectively
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup digitalFilterGroup DigitalFilter
* @brief Digital filter implementation.
* @ingroup classesGroup
* The interface is given by @ref DigitalFilter.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class DigitalFilter;
class DigitalFilterBase;
typedef std::shared_ptr<DigitalFilter> DigitalFilterPtr;

/***** CLASS ***********************************/

/** @brief Digital filter implementation.
* An Instance can be created by @ref readConfig. */
class DigitalFilter
{
private:
  std::vector<DigitalFilterBase*> filters;

public:
  /** @brief Constructor from config. */
  DigitalFilter(Config &config, const std::string &name);

  /// Destructor.
 ~DigitalFilter();

  /** @brief Filter a signal along the first dimension.
  * The signal is filtered along the first axis. Transpose the input to filter
  * along the second. Take caution of warmup effects.
  * @param input   The time series to be filtered.
  * @returns The filtered time series  */
  Matrix filter(const_MatrixSliceRef input) const;

  /** @brief Frequency response evaluated at equally spaced spectral lines.
  * @param length of data in time domain
  * @return complex representation of the frequency response  */
  std::vector<std::complex<Double>> frequencyResponse(UInt length) const;

  /** @brief creates an derived instance of this class. */
  static DigitalFilterPtr create(Config &config, const std::string &name) {return DigitalFilterPtr(new DigitalFilter(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class DigitalFilter.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a class without times is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] digitalFilter Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates DigitalFilter */
template<> Bool readConfig(Config &config, const std::string &name, DigitalFilterPtr &digitalFilter, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***** CLASS ***********************************/

// Internal class
class DigitalFilterBase
{
public:
  enum class PadType {NONE, ZERO, CONSTANT, PERIODIC, SYMMETRIC};

  static Matrix pad (const_MatrixSliceRef input, UInt length, PadType padType);
  static Matrix trim(const_MatrixSliceRef input, UInt length, PadType padType);

public:
  virtual ~DigitalFilterBase() {}
  virtual Matrix filter(const_MatrixSliceRef input) const = 0;
  virtual std::vector<std::complex<Double>> frequencyResponse(UInt length) const = 0;
};

/***** CLASS ***********************************/

// Internal class: default implementation with ARMA coefficients
class DigitalFilterARMA : public DigitalFilterBase
{
protected:
  Vector  an, bn; // AR, MA coefficients
  UInt    bnStartIndex;
  Bool    inFrequencyDomain;
  Bool    backward;
  PadType padType;

  virtual UInt warmup() const;

public:
  DigitalFilterARMA() : inFrequencyDomain(FALSE), backward(FALSE), padType(PadType::NONE) {}
  virtual ~DigitalFilterARMA() {}
  virtual Matrix filter(const_MatrixSliceRef input) const;
  virtual std::vector<std::complex<Double>> frequencyResponse(UInt length) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class PadType.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a class identity filter is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] padType Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates DigitalFilter */
template<> Bool readConfig(Config &config, const std::string &name, DigitalFilterBase::PadType &padType, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

// inline std::string DigitalFilterPadding_typeName() {return "paddingType";}
// void DigitalFilterPadding_registerConfigSchema(Config &config);

/// @}

/***********************************************/

#endif /* __GROOPS_DIGITALFILTER__ */
