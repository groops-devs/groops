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
expressions are also allowed. Declared variables can be accessed via their name. The following
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
\item Condition: \verb|if(c,x,y)|: If the first argument is true (not 0), the second argument is evaluated, otherwise the third.
\end{itemize}
)";
#endif

/***********************************************/

#include <map>
#include "base/importStd.h"

/***** TYPES ***********************************/

class Expression;
class ExpressionVariable;
typedef std::shared_ptr<Expression> ExpressionPtr;
typedef std::shared_ptr<ExpressionVariable> ExpressionVariablePtr;

/***** CLASS ***********************************/

/** @brief List of variables for expressions
* With smart memmory management. The variables are only copied if needed (copy on write, late copy).
* This is not threat save.
* @ingroup parserGroup
* @see ExpressionVariable */
class VariableList
{
  mutable std::map<std::string, ExpressionVariablePtr> map;

  ExpressionVariablePtr getVariable(const std::string &name);

public:
  VariableList() {}                                         //!< Constructor.
  VariableList(const VariableList &) = default;             //!< Copy constructor.
 ~VariableList() {}                                         //!< Destructor
  VariableList &operator=(const VariableList &) = default;  //!< Assignment.

  /// Concatenate.
  VariableList &operator+=(const VariableList &x);

  void addVariable(ExpressionVariablePtr var);

  /** @brief Set @a value of a variable with @a name. The variable is created when needed.*/
  void setVariable(const std::string &name, Double value);

  /** @brief Set @a text of a variable with @a name. The variable is created when needed.*/
  void setVariable(const std::string &name, const std::string &text);

  /** @brief Set status of a variable with @a name undefined. The variable is created when needed.*/
  void undefineVariable(const std::string &name);

  /** @brief Remove the variable with @a name.*/
  void eraseVariable(const std::string &name) {map.erase(name);}

  /** @brief Find variable with @a name.*/
  std::shared_ptr<const ExpressionVariable> find(const std::string &name) const;

  friend class ExpressionVariable;
};

/***** CLASS ***********************************/

/** @brief Variable for expressions
* @ingroup parserGroup
* Set a value to variable @a name to evaluate expressions.
* @see Expression */
class ExpressionVariable
{
public:
  class Func
  {
  public:
    virtual ~Func() {}
    virtual Double operator()() const = 0;
  };

private:
  enum Status {UNDEFINED, TEXT, EXPRESSION, VALUE, FUNC, CIRCULAR};

  std::string    _name;
  mutable Status status;
  std::string    text;
  mutable Double value;
  ExpressionPtr  expr;
  VariableList   varList;
  mutable std::shared_ptr<Func> func;
  Bool           isSimplified;

  /** @brief The unparsed content of the variable. */
  std::string getText() const;

public:
  /** @brief Constructor: variable with undefined value. */
  ExpressionVariable(const std::string &name);

  /** @brief Constructor: variable with value. */
  ExpressionVariable(const std::string &name, Double value);

  /** @brief Constructor: variable with parseable expression text. */
  ExpressionVariable(const std::string &name, const std::string &text);

  /** @brief Constructor: variable with parseable expression text. */
  ExpressionVariable(const std::string &name, const std::string &text, const VariableList &varList);

  /** @brief Constructor: variable with expression. */
  ExpressionVariable(const std::string &name, ExpressionPtr expr, const VariableList &varList);

  /** @brief Constructor: variable value computed with @a func.
  * The @a func is evaluated only if needed and only once. For computational expensive calculations. */
  ExpressionVariable(const std::string &name, const std::shared_ptr<Func> &func);

  /// Copy constructor.
  ExpressionVariable(const ExpressionVariable &x);

  ExpressionVariable &operator=(const ExpressionVariable &) = delete;

  /** @brief Parse an expression given as string.
   * Example: "3+5*sin(1.0)+x^2". */
  static Double parse(const std::string &text, const VariableList &varList);

  /** @brief Name of the variable. */
  const std::string &name() const {return _name;}

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
  void parseVariableName();

  /** @brief Evaluate as much as possible.
  * The @a StringParser with the @a varList is called before and
  * if not all variables {names} can be resolved an exception is thrown.
  * Should be called before evaluation of the expression for long data lists
  * to accelerate the computation. Variables which changes the value after
  * this call (e.g. data0) must be undefined in the @a varList. */
  void simplify(VariableList &varList);

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

  /** @brief Returns the result of the @a StringParser.
  * @a resolved is set to FALSE if not all variables defined. Untouched by success. */
  std::string getParsedText(const VariableList &varList, Bool &resolved) const;
};

/***********************************************/

#endif /* __GROOPS__ */
