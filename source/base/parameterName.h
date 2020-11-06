/***********************************************/
/**
* @file parameterName.h
*
* @brief Parameter name representation.
*
* @author Sebastian Strasser
* @date 2017-05-23
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAME__
#define __GROOPS_PARAMETERNAME__

#include "base/importStd.h"
#include "base/time.h"

/***** CLASS ***********************************/

/**
* @brief Parameter name representation.
* @ingroup base */
class ParameterName
{
public:
  std::string object;    //!< Object this parameter refers to, e.g. graceA, G023, earth, ...
  std::string type;      //!< Type of this parameter, e.g. accBias, position.x, ...
  std::string temporal;  //!< Temporal representation of this parameter, e.g. trend, polynomial.degree1, ...
  std::string interval;  //!< Interval/epoch this parameter represents, e.g. 2017-01-01_00-00-00_2017-01-02_00-00-00, 2018-01-01_00-00-00

  /** @brief Constructor. All parameters are optional.
  * @param object    Object this parameter refers to, e.g. graceA, G023, earth, ...
  * @param type      Type of this parameter, e.g. accBias, position.x, ...
  * @param temporal  Temporal representation of this parameter, e.g. trend, polynomial.degree1, ...
  * @param timeStart Start of interval or point in time which this parameter represents
  * @param timeEnd   End of interval which this parameter represents  */
  ParameterName(const std::string &object="", const std::string &type="", const std::string &temporal="", const Time &timeStart=Time(), const Time &timeEnd=Time());

  /** @brief Constructor.
  * @param object   Object this parameter refers to, e.g. graceA, G023, earth, ...
  * @param type     Type of this parameter, e.g. accBias, position.x, ...
  * @param temporal Temporal representation of this parameter, e.g. trend, polynomial.degree1, ...
  * @param interval Interval/epoch this parameter represents, e.g. 2017-01-01_00-00-00_2017-01-02_00-00-00, 2018-01-01_00-00-00 */
  ParameterName(const std::string &object, const std::string &type, const std::string &temporal, const std::string &interval);

  /** @brief Returns the separator between parts of the parameter name. */
  static constexpr Char sep = ':';

  /** @brief Returns ParameterName object built from a single string @p str separated by a separator @p sep. */
  static ParameterName fromStr(const std::string &str);

  /** @brief Returns string representation of the ParameterName with parts separated by separator @p sep. */
  std::string str() const { return object + sep + type + sep + temporal + sep + interval; }

  /** @brief Fill empty parts of this parameter with parts from @p name. Returns FALSE if parameter names cannot be combined, TRUE otherwise. */
  Bool combine(const ParameterName &name);

  /** @brief Returns TRUE if the parameter name is empty, FALSE otherwise. */
  Bool empty() const { return (object.empty() && type.empty() && interval.empty() && temporal.empty()); }

  /** @brief Returns TRUE if parameter names match. Empty parts always match with non-empty parts. */
  Bool fuzzyMatch(const ParameterName &other) const;

  /** @brief Returns TRUE if parameter names match exactly. */
  Bool operator==(const ParameterName &other) const;

  /** @brief Returns TRUE if parameter names do not match exactly. */
  Bool operator!=(const ParameterName &other) const { return !(*this == other); }

  /** @brief Sort operator. */
  Bool operator<(const ParameterName &other) const;
};

/***********************************************/

#endif

