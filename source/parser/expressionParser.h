/***********************************************/
/**
* @file expressionParser.h
*
* @brief Mathemtical expression parser
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @author Andreas Kvas
* @date 2009-05-17
*
*/
/***********************************************/

#ifndef __GROOPS_EXPRESSIONPARSER__
#define __GROOPS_EXPRESSIONPARSER__

// Latex documentation
#ifdef DOCSTRING_Parser
static const char *docstringParserExpression = R"(
\subsection{Mathematical expression parser}\label{general.parser:expression}
In all input fields that accept numbers (int, uint, double, angle, time) numerical
expressions are also allowed. Values from the global section can be used as variables. The following
operations and functions are defined:
\begin{itemize}
\item Constants:    \verb|pi()|, \verb|rho()=180/pi()|, \verb|nan()|, \verb|c()|: light velocity,
                    \verb|G()|: gravitational constant, \verb|GM()|: gravitational constant of the Earth, \verb|R()|: reference radius of the Earth
\item Mathematical: \verb|+|, \verb|-|, \verb|*|, \verb|/|, \verb|^|
\item Comparison:   \verb|==|, \verb|!=|, \verb|<|, \verb|<=|, \verb|>|, \verb|>=|, result is 1 or 0
\item Logical:      not \verb|!|, and \verb|&&|, \verb'||', or \verb|isnan(x)|, result is 1 or 0
\item Functions:    \verb|sqrt(x)|, \verb|exp(x)|,
                    \verb|sin(x)|,  \verb|cos(x)|, \verb|tan(x)|,
                    \verb|asin(x)|,  \verb|acos(x)|,  \verb|atan(x)|,
                    \verb|abs(x)|,  \verb|round(x)|,  \verb|ceil(x)|,  \verb|floor(x)|,
                    \verb|deg2rad(x)|, \verb|rad2deg(x)|
\item Functions with 2 arguments: \verb|atan2(y,x)|, \verb|min(x,y)|, \verb|max(x,y)|, \verb|mod(x,y)|
\item Time functions: \verb|now()|: local time in MJD, \verb|date2mjd(year, month, day)|, \verb|gps2utc(mjd)|, \verb|utc2gps(mjd)|, \verb|dayofyear(mjd)|, \verb|decimalyear(mjd)|
\item Condition: \verb|if(c,x,y)|: If the first argument is true (not 0), the second argument is evalutated, otherwise the third.
\end{itemize}
)";
#endif

/***********************************************/

#include <map>
#include "base/importStd.h"

/***** TYPES ***********************************/

class Expression;
class ExpressionVariable;
class VariableList;
typedef std::shared_ptr<Expression> ExpressionPtr;
typedef std::shared_ptr<ExpressionVariable> ExpressionVariablePtr;

/***** CLASS ***********************************/

/** @brief Mathematcial Expression
* @ingroup parserGroup
* Represented as a tree.
* Can be created by the standard mathemtical operators.
* Example:
* @code ExpressionPtr expr = exprValue(5)+(sin(exprVar(x))^2); @endcode
* expr->string() results in "5+sin(x)^2" */
class Expression
{
public:
  enum Priority : UInt {NONE, LOGICAL_OR, LOGICAL_AND, EQUALITY, RELATION, ADDITIVE, MULTIPLICATIVE, UNARY, EXPONENTIAL, FUNCTION, VALUE};

  /// Destructor
  virtual ~Expression() {}

  Expression() = default;
  Expression(const Expression &x) = delete;
  Expression &operator=(const Expression &x) = delete;

  /** @brief Parse an expression given as string.
  * Example: "3+5*sin(1.0)+x^2".
  * @return expression represented as a tree. */
  static ExpressionPtr parse(const std::string &text);

  /** @brief Add all used variable names to @a usedName. */
  virtual void usedVariables(const VariableList &varList, std::set<std::string> &usedName) const = 0;

  /** @brief Simplifies expression.
  * @param varList values of the variables contained in the expression.
  * @param[out] resolved expression can be evaluated directly (is a constant value).
  * @return Simplified expression. */
  virtual ExpressionPtr simplify(const VariableList &varList, Bool &resolved) const = 0;

  /** @brief Calculate thep result of an expression.
  * Example: "5*x" results in 10 if variable x=5 is given.
  * @param varList values of the variables contained in the expression. */
  virtual Double evaluate(const VariableList &varList) const = 0;

  /** @brief Derivative of an expression.
  * @param var derivative with respect to this variable. */
  virtual ExpressionPtr derivative(const std::string &var) const = 0;

  /** @brief Deep copy of an expression. */
  virtual ExpressionPtr clone() const = 0;

  /** @brief String representation of an expression. */
  virtual std::string string() const = 0;

  /// Internal.
  std::string string(UInt aprio) const;

  /// Internal.
  virtual UInt priority() const = 0;
};

/***** CLASS ***********************************/

/** @brief Variable for expressions
* @ingroup parserGroup
* Set a value to variable @a name to evaluate expressions.
* @see Expression */
class ExpressionVariable
{
  enum Status {UNDEFINED, TEXT, EXPRESSION, VALUE, CIRCULAR};

  std::string   _name;
  Status        status;
  std::string   text;
  Double        value;
  ExpressionPtr expr;

public:
  /** @brief Constructor: variable with undefined value. */
  ExpressionVariable(const std::string &name="_undefined_");

  /** @brief Constructor: variable with value. */
  ExpressionVariable(const std::string &name, Double value);

  /** @brief Constructor: variable with parseable expression text. */
  ExpressionVariable(const std::string &name, const std::string &text);

  /** @brief Constructor: variable with expression. */
  ExpressionVariable(const std::string &name, ExpressionPtr expr);

  /// Copy constructor.
  ExpressionVariable(const ExpressionVariable &x);

  /// Assignement.
  ExpressionVariable &operator=(const ExpressionVariable &x);

  /** @brief Name of the variable. */
  const std::string &name() const {return _name;}

  /** @brief The unparsed content of the variable. */
  std::string getText() const;

  /** @brief set value of this variable. */
  void setValue(Double value_) {status = VALUE; value = value_;}

  /** @brief set parseable expression text. */
  void setValue(const std::string &text_) {status = TEXT; text = text_;}

  /** @brief set status of this variable undefined. */
  void setUndefined() {status = UNDEFINED;}

  /** @brief Set the name of the variable from the expression text.
  * The text of the variable must be given in the form 'name [= expr]'.
  * The variable @a name is set to 'name' and removed from the text.
  * The 'expr' is not evaluated. If 'expr' is not given the value is set to zero.
  * The @a StringParser with the @a varList is called before. */
  void parseVariableName(const VariableList &varList);

  /** @brief Add all used variable names to @a usedName. */
  void usedVariables(const VariableList &varList, std::set<std::string> &usedName);

  /** @brief Evaluate as much as possible.
  * The @a StringParser with the @a varList is called before and
  * if not all variables {names} can be resolved an exception is thrown.
  * Should be called before evaluation of the expression for long data lists
  * to accelerate the computation. Variables which changes the value after
  * this call (e.g. data0) must be undefined
  * or not in the @a varList when this function is called. */
  void simplify(const VariableList &varList);

  /** @brief Derivative of an expression.
  * The @a StringParser with the @a varList is called before and
  * if not all variables {names} can be resolved an exception is thrown. */
  ExpressionVariablePtr derivative(const std::string &varName, const VariableList &varList) const;

  /** @brief Calculate thep result of an expression.
  * Example: "5*x" results in 10 if variable x=5 is given.
  * The @a StringParser with the @a varList is called before and
  * if not all variables {names} can be resolved an exception is thrown.
  * @param varList values of the variables contained in the expression. */
  Double evaluate(const VariableList &varList) const;

  /** @brief Returns the result of the @a StringParser. */
  std::string getParsedText(const VariableList &varList, Bool &resolved) const;
};


/***** CLASS ***********************************/

/** @brief List of variables for expressions
* @ingroup parserGroup
* @see ExpressionVariable */
class VariableList
{
  std::map<std::string, ExpressionVariablePtr> map;

public:
  VariableList() {}
 ~VariableList() {}

  /// Copy constructor.
  VariableList(const VariableList &x);

  /// Assignement.
  VariableList &operator=(const VariableList &x);

  /// Concatenate.
  VariableList &operator+=(const VariableList &x);

  ExpressionVariablePtr operator[](const std::string &name);
  ExpressionVariablePtr find(const std::string &name) const;
  ExpressionVariablePtr addVariable(ExpressionVariablePtr var);
  void                  clear() { map.clear(); }
  UInt                  erase(const std::string &name) { return map.erase(name); }
  auto                  begin() const { return map.begin(); } // BEWARE: VariableList key and ExpressionVariable name must be identical
  auto                  end()   const { return map.end(); }
};

/***** FUNCTIONS *******************************/

/** @brief Add a variable with undefined value to the variable list.
* @relates VariableList */
void addVariable(const std::string &name, VariableList &varList);

/** @brief Add a variable with value to the variable list.
* @relates VariableList */
void addVariable(const std::string &name, Double value, VariableList &varList);

/** @brief Add a variable with a parseable text to the variable list.
* @relates VariableList */
void addVariable(const std::string &name, const std::string &text, VariableList &varList);

/** @brief Add a variable to the variable list.
* @relates VariableList */
void addVariable(ExpressionVariablePtr var, VariableList &varList);

/***********************************************/

#endif /* __GROOPS__ */
