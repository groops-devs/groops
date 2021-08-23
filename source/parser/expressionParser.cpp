/***********************************************/
/**
* @file expressionParser.cpp
*
* @brief mathemtical expression parser
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @author Andreas Kvas
* @date 2009-05-17
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/time.h"
#include "parser/stringParser.h"
#include "parser/expressionParser.h"

/***********************************************/

ExpressionPtr sqrt(const ExpressionPtr &ob);
ExpressionPtr exp (const ExpressionPtr &ob);
ExpressionPtr log (const ExpressionPtr &ob);
ExpressionPtr sin (const ExpressionPtr &ob);
ExpressionPtr cos (const ExpressionPtr &ob);
ExpressionPtr tan (const ExpressionPtr &ob);
ExpressionPtr asin(const ExpressionPtr &ob);
ExpressionPtr acos(const ExpressionPtr &ob);
ExpressionPtr atan(const ExpressionPtr &ob);
ExpressionPtr operator+(const ExpressionPtr &l, const ExpressionPtr &r);
ExpressionPtr operator-(const ExpressionPtr &l, const ExpressionPtr &r);
ExpressionPtr operator*(const ExpressionPtr &l, const ExpressionPtr &r);
ExpressionPtr operator/(const ExpressionPtr &l, const ExpressionPtr &r);
ExpressionPtr operator^(const ExpressionPtr &l, const ExpressionPtr &r);
ExpressionPtr operator-(const ExpressionPtr &ob);

/***********************************************/

// Generate named list of constants and functions with one, two or three parameters
class ExpressionFunction0;
class ExpressionFunction1;
class ExpressionFunction2;
class ExpressionFunction3;

class FunctionList
{
public:
  static std::vector<std::shared_ptr<ExpressionFunction0>> func0List;
  static std::vector<std::shared_ptr<ExpressionFunction1>> func1List;
  static std::vector<std::shared_ptr<ExpressionFunction2>> func2List;
  static std::vector<std::shared_ptr<ExpressionFunction3>> func3List;
};

std::vector<std::shared_ptr<ExpressionFunction0>> FunctionList::func0List;
std::vector<std::shared_ptr<ExpressionFunction1>> FunctionList::func1List;
std::vector<std::shared_ptr<ExpressionFunction2>> FunctionList::func2List;
std::vector<std::shared_ptr<ExpressionFunction3>> FunctionList::func3List;

#define FUNCLIST0(_name) class FuncList##_name : public FunctionList\
{public: FuncList##_name() {func0List.push_back(std::make_shared<_name>());}};\
static FuncList##_name funcList##_name;

#define FUNCLIST1(_name) class FuncList##_name : public FunctionList\
{public: FuncList##_name() {func1List.push_back(std::make_shared<_name>(ExpressionPtr()));}};\
static FuncList##_name funcList##_name;

#define FUNCLIST2(_name) class FuncList##_name : public FunctionList\
{public: FuncList##_name() {func2List.push_back(std::make_shared<_name>(ExpressionPtr(), ExpressionPtr()));}};\
static FuncList##_name funcList##_name;

#define FUNCLIST3(_name) class FuncList##_name : public FunctionList\
{public: FuncList##_name() {func3List.push_back(std::make_shared<_name>(ExpressionPtr(), ExpressionPtr(), ExpressionPtr()));}};\
static FuncList##_name funcList##_name;


/***********************************************/

class ExpressionValue : public Expression
{
  Double value;

public:
  explicit ExpressionValue(Double v) : value(v) {}
  std::string   string() const override;
  void          usedVariables(const VariableList &/*varList*/, std::set<std::string> &/*usedName*/) const override {}
  ExpressionPtr simplify(const VariableList &/*varList*/, Bool &resolved) const override {resolved = TRUE; return clone();}
  Double        evaluate(const VariableList &/*varList*/) const override {return value;}
  ExpressionPtr derivative(const std::string &/*var*/) const override {return std::make_shared<ExpressionValue>(0);}
  ExpressionPtr clone() const override {return std::make_shared<ExpressionValue>(value);}
  UInt          priority() const override {return Expression::Priority::VALUE;}
};

inline ExpressionPtr exprValue(Double v) {return std::make_shared<ExpressionValue>(v);}

/***********************************************/

class ExpressionVar : public Expression
{
  std::string name;

public:
  explicit ExpressionVar(const std::string &_name) : name(_name) {}
  std::string   string() const override {return name;}
  void          usedVariables(const VariableList &varList, std::set<std::string> &usedName) const override;
  ExpressionPtr simplify(const VariableList &varList, Bool &resolved) const override;
  Double        evaluate(const VariableList &varList) const override;
  ExpressionPtr derivative(const std::string &var) const override  {return exprValue((var == this->name) ? 1 : 0);}
  ExpressionPtr clone() const override {return std::make_shared<ExpressionVar>(name);}
  UInt          priority() const override {return Expression::Priority::VALUE;}
};

inline ExpressionPtr exprVar(const std::string &name) {return std::make_shared<ExpressionVar>(name);}

/***********************************************/
/***********************************************/

// Template for constants
class ExpressionFunction0 : public Expression
{
protected:
  Double value;
  ExpressionFunction0(Double v=0) : value(v) {}

public:
  virtual ~ExpressionFunction0() {}
  virtual std::string name() const = 0;
  std::string   string() const override {return name()+"()";}
  virtual ExpressionPtr create() const = 0;
  void          usedVariables(const VariableList &/*varList*/, std::set<std::string> &/*usedName*/) const override {}
  ExpressionPtr simplify(const VariableList &/*varList*/, Bool &resolved) const override {resolved = TRUE; return clone();}
  Double        evaluate(const VariableList &/*varList*/) const override {return value;}
  ExpressionPtr derivative(const std::string &/*var*/) const override {return exprValue(0);}
  ExpressionPtr clone()  const override {return create();}
  UInt          priority() const override {return Expression::Priority::FUNCTION;}
};

/***********************************************/

// Template for functions with one argument
class ExpressionFunction1 : public Expression
{
protected:
  ExpressionPtr operand;
  ExpressionFunction1(const ExpressionPtr &op) : operand(op) {}

public:
  virtual ~ExpressionFunction1() {}
  virtual std::string   name() const {return "operator";}
  virtual std::string   string() const override {return name()+"("+operand->string()+")";}
  virtual ExpressionPtr create(const ExpressionPtr &ob) const = 0;
  virtual void          usedVariables(const VariableList &varList, std::set<std::string> &usedName) const override {operand->usedVariables(varList, usedName);}
  virtual ExpressionPtr simplify(const VariableList &varList, Bool &resolved) const override;
  virtual ExpressionPtr derivative(const std::string &/*var*/) const override {throw(Exception("Derivative not defined for \""+name()+"\"."));}
  virtual ExpressionPtr clone()  const override {return create(operand->clone());}
  virtual UInt          priority() const override {return Expression::Priority::FUNCTION;}
};

/***********************************************/

// Template for functions with two arguments
class ExpressionFunction2 : public Expression
{
protected:
  ExpressionPtr left, right;
  ExpressionFunction2(const ExpressionPtr &l, const ExpressionPtr &r) : left(l), right(r) {}

public:
  virtual ~ExpressionFunction2() {}
  virtual std::string   name() const {return "operator";}
  virtual std::string   string() const override {return name()+"("+left->string()+", "+right->string()+")";}
  virtual ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const = 0;
  virtual void          usedVariables(const VariableList &varList, std::set<std::string> &usedName) const override {left->usedVariables(varList, usedName); right->usedVariables(varList, usedName);}
  virtual ExpressionPtr simplify(const VariableList &varList, Bool &resolved) const override;
  virtual ExpressionPtr derivative(const std::string &/*var*/) const override {throw(Exception("Derivative not defined for \""+name()+"\"."));}
  virtual ExpressionPtr clone()  const override {return create(left->clone(), right->clone());}
  virtual UInt priority() const override {return Expression::Priority::FUNCTION;}
};

/***********************************************/

// Template for functions with three arguments
class ExpressionFunction3 : public Expression
{
protected:
  ExpressionPtr arg1, arg2, arg3;
  ExpressionFunction3(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) : arg1(a1), arg2(a2), arg3(a3) {}

public:
  virtual ~ExpressionFunction3() {}
  virtual std::string   name() const {return "operator";}
  virtual std::string   string() const override {return name()+"("+arg1->string()+", "+arg2->string()+", "+arg3->string()+")";}
  virtual ExpressionPtr create(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) const = 0;
  virtual void          usedVariables(const VariableList &varList, std::set<std::string> &usedName) const override {arg1->usedVariables(varList, usedName); arg2->usedVariables(varList, usedName); arg3->usedVariables(varList, usedName);}
  virtual ExpressionPtr simplify(const VariableList &varList, Bool &resolved) const override;
  virtual ExpressionPtr derivative(const std::string &/*var*/) const override {throw(Exception("Derivative not defined for \""+name()+"\"."));}
  virtual ExpressionPtr clone()  const override {return create(arg1->clone(), arg2->clone(), arg3->clone());}
  virtual UInt priority() const override {return Expression::Priority::FUNCTION;}
};

/***********************************************/
/*** constants *********************************/
/***********************************************/

class ExpressionPi : public ExpressionFunction0
{
public:
  ExpressionPi() : ExpressionFunction0(PI) {}
  std::string   name()   const override {return "pi";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionPi>();}
};
FUNCLIST0(ExpressionPi)

inline ExpressionPtr exprPi() {return std::make_shared<ExpressionPi>();}

/***********************************************/

class ExpressionRho : public ExpressionFunction0
{
public:
  ExpressionRho() : ExpressionFunction0(RAD2DEG) {}
  std::string   name()   const override {return "rho";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionRho>();}
};
FUNCLIST0(ExpressionRho)

inline ExpressionPtr exprRho() {return std::make_shared<ExpressionRho>();}

/***********************************************/

class ExpressionNan : public ExpressionFunction0
{
public:
  ExpressionNan() : ExpressionFunction0(NAN_EXPR) {}
  std::string   name()   const override {return "nan";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionNan>();}
};
FUNCLIST0(ExpressionNan)

inline ExpressionPtr exprNan() {return std::make_shared<ExpressionNan>();}

/***********************************************/

class ExpressionLightVelovity : public ExpressionFunction0
{
public:
  ExpressionLightVelovity() : ExpressionFunction0(LIGHT_VELOCITY) {}
  std::string   name()   const override {return "c";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionLightVelovity>();}
};
FUNCLIST0(ExpressionLightVelovity)

/***********************************************/

class ExpressionGravitationalConstant : public ExpressionFunction0
{
public:
  ExpressionGravitationalConstant() : ExpressionFunction0(GRAVITATIONALCONSTANT) {}
  std::string   name()   const override {return "G";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionGravitationalConstant>();}
};
FUNCLIST0(ExpressionGravitationalConstant)

/***********************************************/

class ExpressionGM : public ExpressionFunction0
{
public:
  ExpressionGM() : ExpressionFunction0(DEFAULT_GM) {}
  std::string   name()   const override {return "GM";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionGM>();}
};
FUNCLIST0(ExpressionGM)

/***********************************************/

class ExpressionR : public ExpressionFunction0
{
public:
  ExpressionR() : ExpressionFunction0(DEFAULT_R) {}
  std::string   name()   const override {return "R";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionR>();}
};
FUNCLIST0(ExpressionR)

/***********************************************/
/*** basic arithmetic operations ***************/
/***********************************************/

class ExpressionNegative : public ExpressionFunction1
{
public:
  explicit ExpressionNegative(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   string() const override {return "-"+operand->string(priority());}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionNegative>(ob);}
  Double        evaluate(const VariableList &v) const override {return -operand->evaluate(v);}
  ExpressionPtr derivative(const std::string &var) const override {return -operand->derivative(var);}
  UInt priority() const override {return Expression::Priority::UNARY;}
};

inline ExpressionPtr operator-(const ExpressionPtr &ob) {return std::make_shared<ExpressionNegative>(ob);}

/***********************************************/

class ExpressionAdd : public ExpressionFunction2
{
public:
  ExpressionAdd(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l, r) {}
  std::string   string() const override {return left->string(priority()) + " + " + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionAdd>(l, r);}
  Double        evaluate(const VariableList &v) const override {return left->evaluate(v) + right->evaluate(v);}
  ExpressionPtr derivative(const std::string &var) const override {return left->derivative(var) + right->derivative(var);}
  UInt priority() const override {return Expression::Priority::ADDITIVE;}
};

inline ExpressionPtr operator+(const ExpressionPtr &l, const ExpressionPtr &r) {return std::make_shared<ExpressionAdd>(l,r);}

/***********************************************/

class ExpressionSub : public ExpressionFunction2
{
public:
  ExpressionSub(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l, r) {}
  std::string   string() const override {return left->string(priority()) + " - " + right->string(priority()+1);}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionSub>(l, r);}
  Double        evaluate(const VariableList &v) const override {return left->evaluate(v) - right->evaluate(v);}
  ExpressionPtr derivative(const std::string &var) const override {return left->derivative(var) - right->derivative(var);}
  UInt priority() const override {return Expression::Priority::ADDITIVE;}
};

inline ExpressionPtr operator-(const ExpressionPtr &l, const ExpressionPtr &r) {return std::make_shared<ExpressionSub>(l,r);}

/***********************************************/

class ExpressionMult : public ExpressionFunction2
{
public:
  ExpressionMult(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l, r) {}
  std::string   string() const override {return left->string(priority()) + "*" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionMult>(l, r);}
  Double        evaluate(const VariableList &v) const override {return left->evaluate(v) * right->evaluate(v);}
  ExpressionPtr derivative(const std::string &var) const override {return left->derivative(var)*right->clone() + left->clone()*right->derivative(var);}
  UInt priority() const override {return Expression::Priority::MULTIPLICATIVE;}
};

inline ExpressionPtr operator*(const ExpressionPtr &l, const ExpressionPtr &r) {return std::make_shared<ExpressionMult>(l,r);}

/***********************************************/

class ExpressionDiv : public ExpressionFunction2
{
public:
  ExpressionDiv(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l, r) {}
  std::string string() const override {return left->string(priority()) + "/" + right->string(priority()+1);}
  inline ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionDiv>(l, r);}
  inline Double        evaluate(const VariableList &v) const override {return left->evaluate(v) / right->evaluate(v);}
  inline ExpressionPtr derivative(const std::string &var) const override {return (left->derivative(var)*right->clone() - left->clone()*right->derivative(var))/(right->clone()^exprValue(2));}
  UInt priority() const override {return Expression::Priority::MULTIPLICATIVE;}
};

inline ExpressionPtr operator/(const ExpressionPtr &l, const ExpressionPtr &r) {return std::make_shared<ExpressionDiv>(l,r);}

/***********************************************/

class ExpressionPow : public ExpressionFunction2
{
public:
  ExpressionPow(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l, r) {}
  std::string   string() const override {return left->string(priority()) + "^" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionPow>(l, r);}
  Double        evaluate(const VariableList &v) const override {return pow(left->evaluate(v),right->evaluate(v));}
  ExpressionPtr derivative(const std::string &var) const override {return left->derivative(var)*right->clone()/left->clone() * (left->clone()^(right->clone()-exprValue(1)));} // this derivative is not completly correct!!!!
  UInt priority() const override {return Expression::Priority::EXPONENTIAL;}
};

inline ExpressionPtr operator^(const ExpressionPtr &l, const ExpressionPtr &r) {return std::make_shared<ExpressionPow>(l,r);}

/***********************************************/
/*** basic comparison operations ***************/
/***********************************************/

class ExpressionLessThan : public ExpressionFunction2
{
public:
  ExpressionLessThan(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator<";}
  std::string   string() const override {return left->string(priority()) + "<" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionLessThan>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) < right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::RELATION;}
};

/***********************************************/

class ExpressionLessEqualThan : public ExpressionFunction2
{
public:
  ExpressionLessEqualThan(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator<=";}
  std::string   string() const override {return left->string(priority()) + "<=" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionLessEqualThan>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) <= right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::RELATION;}
};

/***********************************************/

class ExpressionGreaterThan : public ExpressionFunction2
{
public:
  ExpressionGreaterThan(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator>";}
  std::string   string() const override {return left->string(priority()) + ">" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionGreaterThan>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) > right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::RELATION;}
};

/***********************************************/

class ExpressionGreaterEqualThan : public ExpressionFunction2
{
public:
  ExpressionGreaterEqualThan(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator>=";}
  std::string   string() const override {return left->string(priority()) + ">=" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionGreaterEqualThan>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) >= right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::RELATION;}
};

/***********************************************/

class ExpressionEqual : public ExpressionFunction2
{
public:
  ExpressionEqual(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator==";}
  std::string   string() const override {return left->string(priority()) + "==" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionEqual>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) == right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::EQUALITY;}
};

/***********************************************/

class ExpressionNotEqual : public ExpressionFunction2
{
public:
  ExpressionNotEqual(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator!=";}
  std::string   string() const override {return left->string(priority()) + "!=" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionNotEqual>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return left->evaluate(varList) != right->evaluate(varList) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::EQUALITY;}
};

/***********************************************/
/*** logical operators *************************/
/***********************************************/

class ExpressionLogicalNot : public ExpressionFunction1
{
public:
  explicit ExpressionLogicalNot(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name()   const override {return "operator!";}
  std::string   string() const override {return "!"+operand->string(priority());}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionLogicalNot>(ob);}
  Double        evaluate(const VariableList &v) const override {return operand->evaluate(v) == 0.0 ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::UNARY;}
};

/***********************************************/

class ExpressionLogicalAnd : public ExpressionFunction2
{
public:
  ExpressionLogicalAnd(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator&&";}
  std::string   string() const override {return left->string(priority()) + "&&" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionLogicalAnd>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return (left->evaluate(varList)!=0.0 &&  right->evaluate(varList)!=0.0) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::LOGICAL_AND;}
};

/***********************************************/

class ExpressionLogicalOr : public ExpressionFunction2
{
public:
  ExpressionLogicalOr(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name()   const override {return "operator||";}
  std::string   string() const override {return left->string(priority()) + "||" + right->string(priority());}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionLogicalOr>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return (left->evaluate(varList)!=0.0 || right->evaluate(varList)!=0.0) ? 1.0 : 0.0;}
  UInt priority() const override {return Expression::Priority::LOGICAL_OR;}
};

/***********************************************/

class ExpressionIfThenElse : public ExpressionFunction3
{
public:
  ExpressionIfThenElse(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) : ExpressionFunction3(a1,a2,a3) {}
  std::string   name() const override {return "if";}
  ExpressionPtr create(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) const override {return std::make_shared<ExpressionIfThenElse>(a1, a2, a3);}
  Double        evaluate(const VariableList &varList) const override {return arg1->evaluate(varList) != 0.0 ? arg2->evaluate(varList) : arg3->evaluate(varList);}
};
FUNCLIST3(ExpressionIfThenElse)

/***********************************************/

class ExpressionIsNan : public ExpressionFunction1
{
public:
  explicit ExpressionIsNan(const ExpressionPtr &a1) : ExpressionFunction1(a1) {}
  std::string   name() const override {return "isnan";}
  ExpressionPtr create(const ExpressionPtr &a1) const override {return std::make_shared<ExpressionIsNan>(a1);}
  Double        evaluate(const VariableList &varList) const override {return std::isnan(operand->evaluate(varList));}
};
FUNCLIST1(ExpressionIsNan)

/***********************************************/
/*** functions *********************************/
/***********************************************/

class ExpressionSqrt : public ExpressionFunction1
{
public:
  explicit ExpressionSqrt(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "sqrt";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionSqrt>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::sqrt(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)/sqrt(operand->clone());}
};
FUNCLIST1(ExpressionSqrt)

inline ExpressionPtr sqrt(const ExpressionPtr &ob) {return std::make_shared<ExpressionSqrt>(ob);}

/***********************************************/

class ExpressionExp : public ExpressionFunction1
{
public:
  explicit ExpressionExp(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "exp";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionExp>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::exp(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)*exp(operand->clone());}
};
FUNCLIST1(ExpressionExp)

inline ExpressionPtr exp(const ExpressionPtr &ob) {return std::make_shared<ExpressionExp>(ob);}

/***********************************************/

class ExpressionLog : public ExpressionFunction1
{
public:
  explicit ExpressionLog(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "log";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionExp>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::log(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)^exprValue(-1);}
};
FUNCLIST1(ExpressionLog)

inline ExpressionPtr log(const ExpressionPtr &ob) {return std::make_shared<ExpressionLog>(ob);}

/***********************************************/

class ExpressionSin : public ExpressionFunction1
{
public:
  explicit ExpressionSin(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "sin";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionSin>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::sin(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)*cos(operand->clone());}
};
FUNCLIST1(ExpressionSin)

inline ExpressionPtr sin(const ExpressionPtr &ob) {return std::make_shared<ExpressionSin>(ob);}

/***********************************************/

class ExpressionCos : public ExpressionFunction1
{
public:
  explicit ExpressionCos(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "cos";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionCos>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::cos(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return -operand->derivative(var)*sin(operand->clone());}
};
FUNCLIST1(ExpressionCos)

inline ExpressionPtr cos(const ExpressionPtr &ob) {return std::make_shared<ExpressionCos>(ob);}

/***********************************************/

class ExpressionTan : public ExpressionFunction1
{
public:
  explicit ExpressionTan(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "tan";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionTan>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::tan(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)/(cos(operand->clone())^exprValue(2));}
};
FUNCLIST1(ExpressionTan)

inline ExpressionPtr tan(const ExpressionPtr &ob) {return std::make_shared<ExpressionTan>(ob);}

/***********************************************/

class ExpressionAsin : public ExpressionFunction1
{
public:
  explicit ExpressionAsin(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "asin";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionAsin>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::asin(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)/sqrt(exprValue(1) - ((operand->clone())^exprValue(2)));}
};
FUNCLIST1(ExpressionAsin)

inline ExpressionPtr asin(const ExpressionPtr &ob) {return std::make_shared<ExpressionAsin>(ob);}

/***********************************************/

class ExpressionAcos : public ExpressionFunction1
{
public:
  explicit ExpressionAcos(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "acos";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionAcos>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::acos(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return -operand->derivative(var)/sqrt(exprValue(1) - ((operand->clone())^exprValue(2)));}
};
FUNCLIST1(ExpressionAcos)

inline ExpressionPtr acos(const ExpressionPtr &ob) {return std::make_shared<ExpressionAcos>(ob);}

/***********************************************/

class ExpressionAtan : public ExpressionFunction1
{
public:
  explicit ExpressionAtan(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "atan";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionAtan>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::atan(operand->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return operand->derivative(var)/(exprValue(1) + ((operand->clone())^exprValue(2)));}
};
FUNCLIST1(ExpressionAtan)

inline ExpressionPtr atan(const ExpressionPtr &ob) {return std::make_shared<ExpressionAtan>(ob);}

/***********************************************/

class ExpressionAbs : public ExpressionFunction1
{
public:
  explicit ExpressionAbs(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "abs";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionAbs>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::fabs(operand->evaluate(varList));}
};
FUNCLIST1(ExpressionAbs)

/***********************************************/

class ExpressionRound : public ExpressionFunction1
{
public:
  explicit ExpressionRound(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "round";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionRound>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::round(operand->evaluate(varList));}
};
FUNCLIST1(ExpressionRound)

/***********************************************/

class ExpressionCeil : public ExpressionFunction1
{
public:
  explicit ExpressionCeil(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "ceil";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionCeil>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::ceil(operand->evaluate(varList));}
};
FUNCLIST1(ExpressionCeil)

/***********************************************/

class ExpressionFloor : public ExpressionFunction1
{
public:
  explicit ExpressionFloor(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "floor";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionFloor>(ob);}
  Double        evaluate(const VariableList &varList) const override {return std::floor(operand->evaluate(varList));}
};
FUNCLIST1(ExpressionFloor)

/***********************************************/

class ExpressionDeg2Rad : public ExpressionFunction1
{
public:
  explicit ExpressionDeg2Rad(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "deg2rad";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionDeg2Rad>(ob);}
  Double        evaluate(const VariableList &varList) const override {return DEG2RAD * operand->evaluate(varList);}
};
FUNCLIST1(ExpressionDeg2Rad)

/***********************************************/

class ExpressionRad2Deg : public ExpressionFunction1
{
public:
  explicit ExpressionRad2Deg(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "rad2deg";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionRad2Deg>(ob);}
  Double        evaluate(const VariableList &varList) const override {return RAD2DEG * operand->evaluate(varList);}
};
FUNCLIST1(ExpressionRad2Deg)

/***********************************************/
/***********************************************/

class ExpressionAtan2 : public ExpressionFunction2
{
public:
  ExpressionAtan2(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name() const override {return "atan2";}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionAtan2>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return std::atan2(left->evaluate(varList), right->evaluate(varList));}
  ExpressionPtr derivative(const std::string &var) const override {return (left->derivative(var)*right->clone() - left->clone()*right->derivative(var))/((left->clone()^exprValue(2))+(right->clone()^exprValue(2)));}
};
FUNCLIST2(ExpressionAtan2)

/***********************************************/

class ExpressionMin : public ExpressionFunction2
{
public:
  ExpressionMin(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name() const override {return "min";}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionMin>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return std::min(left->evaluate(varList), right->evaluate(varList));}
};
FUNCLIST2(ExpressionMin)

/***********************************************/

class ExpressionMax : public ExpressionFunction2
{
public:
  ExpressionMax(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name() const override {return "max";}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionMax>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return std::max(left->evaluate(varList), right->evaluate(varList));}
};
FUNCLIST2(ExpressionMax)

/***********************************************/

class ExpressionMod : public ExpressionFunction2
{
public:
  ExpressionMod(const ExpressionPtr &l, const ExpressionPtr &r) : ExpressionFunction2(l,r) {}
  std::string   name() const override {return "mod";}
  ExpressionPtr create(const ExpressionPtr &l, const ExpressionPtr &r) const override {return std::make_shared<ExpressionMod>(l, r);}
  Double        evaluate(const VariableList &varList) const override {return std::fmod(left->evaluate(varList), right->evaluate(varList));}
};
FUNCLIST2(ExpressionMod)

/***********************************************/
/*** time conversions **************************/
/***********************************************/

class ExpressionNow : public ExpressionFunction0
{
public:
  ExpressionNow()
  {
    std::time_t tt = std::time(nullptr);
    std::tm     t  = *std::localtime(&tt);
    value = date2time(t.tm_year+1900, t.tm_mon+1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec).mjd();
  }
  std::string   name()   const override {return "now";}
  ExpressionPtr create() const override {return std::make_shared<ExpressionNow>();}
};
FUNCLIST0(ExpressionNow)

/***********************************************/

class ExpressionDate2mjd : public ExpressionFunction3
{
public:
  ExpressionDate2mjd(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) : ExpressionFunction3(a1,a2,a3) {}
  std::string   name() const override {return "date2mjd";}
  ExpressionPtr create(const ExpressionPtr &a1, const ExpressionPtr &a2, const ExpressionPtr &a3) const override {return std::make_shared<ExpressionDate2mjd>(a1, a2, a3);}
  Double        evaluate(const VariableList &varList) const override {return date2time(static_cast<UInt>(arg1->evaluate(varList)), static_cast<UInt>(arg2->evaluate(varList)), static_cast<UInt>(arg3->evaluate(varList))).mjd();}
};
FUNCLIST3(ExpressionDate2mjd)

/***********************************************/

class ExpressionGps2Utc : public ExpressionFunction1
{
public:
  explicit ExpressionGps2Utc(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "gps2utc";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionGps2Utc>(ob);}
  Double        evaluate(const VariableList &varList) const override {return timeGPS2UTC(mjd2time(operand->evaluate(varList))).mjd();}
};
FUNCLIST1(ExpressionGps2Utc)

/***********************************************/

class ExpressionUtc2Gps : public ExpressionFunction1
{
public:
  explicit ExpressionUtc2Gps(const ExpressionPtr &ob) : ExpressionFunction1(ob) {}
  std::string   name() const override {return "utc2gps";}
  ExpressionPtr create(const ExpressionPtr &ob) const override {return std::make_shared<ExpressionUtc2Gps>(ob);}
  Double        evaluate(const VariableList &varList) const override {return timeUTC2GPS(mjd2time(operand->evaluate(varList))).mjd();}
};
FUNCLIST1(ExpressionUtc2Gps)

/***********************************************/

class ExpressionDayOfYear : public ExpressionFunction1
{
public:
  explicit ExpressionDayOfYear(const ExpressionPtr &a1) : ExpressionFunction1(a1) {}
  std::string   name() const override {return "dayofyear";}
  ExpressionPtr create(const ExpressionPtr &a1) const override {return std::make_shared<ExpressionDayOfYear>(a1);}
  Double        evaluate(const VariableList &varList) const override {return static_cast<Double>(mjd2time(operand->evaluate(varList)).dayOfYear());}
};
FUNCLIST1(ExpressionDayOfYear)

/***********************************************/

class ExpressionDecimalYear : public ExpressionFunction1
{
public:
  explicit ExpressionDecimalYear(const ExpressionPtr &a1) : ExpressionFunction1(a1) {}
  std::string   name() const override {return "decimalyear";}
  ExpressionPtr create(const ExpressionPtr &a1) const override {return std::make_shared<ExpressionDecimalYear>(a1);}
  Double        evaluate(const VariableList &varList) const override {return mjd2time(operand->evaluate(varList)).decimalYear();}
};
FUNCLIST1(ExpressionDecimalYear)

/***********************************************/
/***********************************************/
/***********************************************/

std::string ExpressionValue::string() const
{
  std::stringstream ss;
  ss<<value;
  return ss.str();
}

/***********************************************/

std::string Expression::string(UInt aprio) const
{
  return (priority()<aprio) ? "("+string()+")" : string();
}

/***********************************************/

Double ExpressionVar::evaluate(const VariableList &varList) const
{
  auto variable = varList.find(name);
  if(!variable)
    throw(Exception("unknown variable: "+name));
  return variable->evaluate(varList);
}

/***********************************************/

void ExpressionVar::usedVariables(const VariableList &varList, std::set<std::string> &usedName) const
{
  usedName.insert(name);
  auto variable = varList.find(name);
  if(variable)
    variable->usedVariables(varList, usedName);
}

/***********************************************/

ExpressionPtr ExpressionVar::simplify(const VariableList &varList, Bool &resolved) const
{
  auto variable = varList.find(name);
  if(!variable)
  {
    resolved = FALSE;
    return clone();
  }
  try
  {
    resolved = TRUE;
    return exprValue(variable->evaluate(varList));
  }
  catch(...)
  {
    resolved = FALSE;
    return clone();
  }
}

/***********************************************/

ExpressionPtr ExpressionFunction1::simplify(const VariableList &varList, Bool &resolved) const
{
  ExpressionPtr expr = create(operand->simplify(varList, resolved));
  if(resolved)
    return exprValue(expr->evaluate(varList));
  return expr;
}

/***********************************************/

ExpressionPtr ExpressionFunction2::simplify(const VariableList &varList, Bool &resolved) const
{
  Bool resolved1, resolved2;
  ExpressionPtr expr = create(left->simplify(varList, resolved1), right->simplify(varList, resolved2));
  resolved = resolved1 && resolved2;
  if(resolved)
    return exprValue(expr->evaluate(varList));
  return expr;
}

/***********************************************/

ExpressionPtr ExpressionFunction3::simplify(const VariableList &varList, Bool &resolved) const
{
  Bool resolved1, resolved2, resolved3;
  ExpressionPtr expr = create(arg1->simplify(varList, resolved1), arg2->simplify(varList, resolved2), arg3->simplify(varList, resolved3));
  resolved = (resolved1 && resolved2) && resolved3;
  if(resolved)
    return exprValue(expr->evaluate(varList));
  return expr;
}

/***********************************************/
/***********************************************/

class Tokenizer
{
  std::string text;
  UInt pos, posStart;
  Bool unused;

public:
  enum Type   {END, SINGLE, VALUE, NAME, OPERATOR};
  enum OperatorType   {NONE, ADD, SUB, MULT, DIV, POW, EQUAL, NOT_EQUAL, LESS, GREATER, LESS_EQUAL, GREATER_EQUAL, LOGICAL_AND, LOGICAL_OR};
  Type        type;
  OperatorType opType;
  char        c;
  Double      value;
  std::string name;

  explicit Tokenizer(const std::string &text);

  void get();
  void putBack() {unused = TRUE;}
  std::string infoString(const std::string &info) const;
};

/***********************************************/

Tokenizer::Tokenizer(const std::string &txt) : text(txt), pos(0), posStart(0), unused(FALSE)
{
  get();
  putBack();
}

/***********************************************/

void Tokenizer::get()
{
  if(unused)
  {
    unused=FALSE;
    return;
  }

  name.clear();
  c     = ' ';
  value = 0.0;
  type  = Tokenizer::END;
  opType = Tokenizer::NONE;

  // eat whitespace
  while((pos<text.size()) && isspace(text.at(pos)))
    pos++;

  if(pos>=text.size())
  {
    type = Tokenizer::END;
    return;
  }

  posStart = pos;
  c = text.at(pos);

  // operator? (+,-,..)
  if(!(isalnum(c)||(c=='_')))
  {
    type = Tokenizer::SINGLE;
    if(c == '+') opType = Tokenizer::ADD;
    if(c == '-') opType = Tokenizer::SUB;
    if(c == '*') opType = Tokenizer::MULT;
    if(c == '/') opType = Tokenizer::DIV;
    if(c == '^') opType = Tokenizer::POW;
    if(c == '>')
    {
      opType = Tokenizer::GREATER;
      if(pos+1<text.size() && text.at(pos+1) == '=')
      {
        opType = Tokenizer::GREATER_EQUAL;
        pos++;
      }
    }
    if(c == '<')
    {
      opType = Tokenizer::LESS;
      if(pos+1<text.size() && text.at(pos+1) == '=')
      {
        opType = Tokenizer::LESS_EQUAL;
        pos++;
      }
    }
    if((c == '=') && (pos+1<text.size() && text.at(pos+1) == '='))
    {
      opType = Tokenizer::EQUAL;
      pos++;
    }
    if((c == '!') && (pos+1<text.size() && text.at(pos+1) == '='))
    {
      opType = Tokenizer::NOT_EQUAL;
      pos++;
    }
    if((c == '&') && (pos+1<text.size() && text.at(pos+1) == '&'))
    {
      opType = Tokenizer::LOGICAL_AND;
      pos++;
    }
    if((c == '|') && (pos+1<text.size() && text.at(pos+1) == '|'))
    {
      opType = Tokenizer::LOGICAL_OR;
      pos++;
    }
    pos++;
    return;
  }

  // number?
  if(std::isdigit(c))
  {
    type  = Tokenizer::VALUE;
    const char *ptr1 = text.c_str()+pos;
    char *ptr2;
    value = std::strtod(ptr1, &ptr2);
    pos += ptr2-ptr1;
    return;
  }

  // name
  type  = Tokenizer::NAME;
  while((pos<text.size()) && (isalnum(text.at(pos))||(text.at(pos)=='_')))
    name.push_back(text.at(pos++));
}

/***********************************************/

std::string Tokenizer::infoString(const std::string &info) const
{
  std::stringstream ss;
  ss<<"'"<<text<<"', col="<<posStart<<", "<<info<<std::endl
    <<std::setfill(' ')<<std::setw(posStart+2)<<std::right<<'^';
  return ss.str();
}

/***********************************************/

static ExpressionPtr parse(Tokenizer &token, UInt priority);
static ExpressionPtr operand(Tokenizer &token);

/***********************************************/

ExpressionPtr operand(Tokenizer &token)
{
  token.get();
  if(token.type==Tokenizer::SINGLE)
  {
    switch(token.c)
    {
    case '+':
      return  parse(token, Expression::Priority::UNARY);
    case '-':
      return -parse(token, Expression::Priority::UNARY);
    case '!':
      return std::make_shared<ExpressionLogicalNot>(parse(token, Expression::Priority::UNARY));
    case '(':
      ExpressionPtr expr = parse(token, 0);
      token.get();
      if((token.type!=Tokenizer::SINGLE)||(token.c!=')'))
        throw(Exception(token.infoString("expected ')'")));
      return expr;
    }
    throw(Exception(token.infoString("unknown token")));
  }
  else if(token.type==Tokenizer::VALUE)
  {
    return exprValue(token.value);
  }
  else if(token.type==Tokenizer::NAME)
  {
    const std::string name = token.name;
    token.get();
    if((token.type!=Tokenizer::SINGLE)||(token.c!='('))
    {
      token.putBack();
      return exprVar(name); // new variable
    }

    // is function/constant without arguments?
    for(UInt i=0; i<FunctionList::func0List.size(); i++)
      if(name == FunctionList::func0List.at(i)->name())
      {
        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=')'))
          throw(Exception(token.infoString("expected ')'")));
        return FunctionList::func0List.at(i)->create();
      }

    // is Function with one argument?
    for(UInt i=0; i<FunctionList::func1List.size(); i++)
      if(name == FunctionList::func1List.at(i)->name())
      {
        ExpressionPtr expr = FunctionList::func1List.at(i)->create(parse(token, 0));
        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=')'))
          throw(Exception(token.infoString("expected ')'")));
        return expr;
      }

    // is Function with two arguments?
    for(UInt i=0; i<FunctionList::func2List.size(); i++)
      if(name == FunctionList::func2List.at(i)->name())
      {
        ExpressionPtr ob1 = parse(token, 0);

        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=','))
          throw(Exception(token.infoString("expected ','")));

        ExpressionPtr ob2 = parse(token, 0);

        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=')'))
          throw(Exception(token.infoString("expected ')'")));

        return FunctionList::func2List.at(i)->create(ob1, ob2);
      }

    // is Function with three arguments?
    for(UInt i=0; i<FunctionList::func3List.size(); i++)
      if(name == FunctionList::func3List.at(i)->name())
      {
        ExpressionPtr ob1 = parse(token, 0);

        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=','))
          throw(Exception(token.infoString("expected ','")));

        ExpressionPtr ob2 = parse(token, 0);

        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=','))
          throw(Exception(token.infoString("expected ','")));

        ExpressionPtr ob3 = parse(token, 0);

        token.get();
        if((token.type!=Tokenizer::SINGLE)||(token.c!=')'))
          throw(Exception(token.infoString("expected ')'")));

        return FunctionList::func3List.at(i)->create(ob1, ob2, ob3);
      }

    throw(Exception(token.infoString("unknown function")));
  }

  throw(Exception(token.infoString("expected operand")));
}

/***********************************************/

ExpressionPtr parse(Tokenizer &token, UInt priority)
{
  UInt priorityNew = Expression::Priority::NONE;
  ExpressionPtr ob = operand(token);

  do
  {
    token.get();
    if(token.type != Tokenizer::SINGLE)
      break;
    switch(token.opType)
    {
    case Tokenizer::ADD:
      priorityNew = Expression::Priority::ADDITIVE;
      if(priorityNew > priority)
        ob = ob + parse(token, priorityNew);
      break;
    case Tokenizer::SUB:
      priorityNew = Expression::Priority::ADDITIVE;
      if(priorityNew > priority)
        ob = ob - parse(token, priorityNew);
      break;
    case Tokenizer::MULT:
      priorityNew = Expression::Priority::MULTIPLICATIVE;
      if(priorityNew > priority)
        ob = ob * parse(token, priorityNew);
      break;
    case Tokenizer::DIV:
      priorityNew = Expression::Priority::MULTIPLICATIVE;
      if(priorityNew > priority)
        ob = ob / parse(token, priorityNew);
      break;
    case Tokenizer::POW:
      priorityNew = Expression::Priority::EXPONENTIAL;
      if(priorityNew > priority)
        ob = ob ^ parse(token, priorityNew);
      break;
    case Tokenizer::GREATER:
      priorityNew = Expression::Priority::RELATION;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionGreaterThan>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::LESS:
      priorityNew = Expression::Priority::RELATION;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionLessThan>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::LESS_EQUAL:
      priorityNew = Expression::Priority::RELATION;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionLessEqualThan>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::GREATER_EQUAL:
      priorityNew = Expression::Priority::RELATION;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionGreaterEqualThan>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::EQUAL:
      priorityNew = Expression::Priority::EQUALITY;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionEqual>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::NOT_EQUAL:
      priorityNew = Expression::Priority::EQUALITY;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionNotEqual>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::LOGICAL_AND:
      priorityNew = Expression::Priority::LOGICAL_AND;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionLogicalAnd>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::LOGICAL_OR:
      priorityNew = Expression::Priority::LOGICAL_OR;
      if(priorityNew > priority)
        ob = std::make_shared<ExpressionLogicalOr>(ob, parse(token, priorityNew));
      break;
    case Tokenizer::NONE:
      priorityNew = Expression::Priority::NONE;
      break;
    }
  }
  while(priority<priorityNew);

  token.putBack();
  if(!ob)
    throw(Exception(token.infoString("unknown error")));
  return ob;
}

/***********************************************/

ExpressionPtr Expression::parse(const std::string &text)
{
  try
  {
    Tokenizer token(text);
    ExpressionPtr expr = ::parse(token, 0);
    if(token.type!=Tokenizer::END)
      throw(Exception(token.infoString("unknown token")));
    return expr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

ExpressionVariable::ExpressionVariable(const std::string &name)
  : _name(name), status(UNDEFINED) {}

ExpressionVariable::ExpressionVariable(const std::string &name, Double _value)
  : _name(name), status(VALUE), value(_value) {}

ExpressionVariable::ExpressionVariable(const std::string &name, const std::string &_text)
  : _name(name), status(TEXT), text(_text) {}

ExpressionVariable::ExpressionVariable(const std::string &name, ExpressionPtr _expr)
  : _name(name), status(EXPRESSION), expr(_expr) {if(!expr) throw(Exception("ExpressionVariable: Null-Pointer."));}

ExpressionVariable::ExpressionVariable(const ExpressionVariable &x)
  : _name(x._name), status(x.status), text(x.text), value(x.value) {if(x.expr) expr = x.expr->clone();}


/***********************************************/

ExpressionVariable &ExpressionVariable::operator=(const ExpressionVariable &x)
{
  _name  = x._name;
  status = x.status;
  text   = x.text;
  value  = x.value;
  expr   = (x.expr) ? x.expr->clone() : x.expr;
  return *this;
}

/***********************************************/

std::string ExpressionVariable::getText() const
{
  try
  {
    if(status == VALUE)
    {
      std::stringstream ss;
      ss<<value;
      return ss.str();
    }

    if(status == EXPRESSION)
      return expr->string();

    if(status == UNDEFINED)
      return "{_undefined_}";

    // status == TEXT
    return text;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+"')", e)
  }
}

/***********************************************/

void ExpressionVariable::parseVariableName(const VariableList &varList)
{
  try
  {
    if(status != TEXT)
      throw(Exception("must contain a variable name ('name [= expr]')"));

    Bool resolved;
    auto text_ = StringParser::parse(name(), this->text, varList, resolved);

    // parse variable name
    auto pos     = text_.find('=');
    auto nameStr = text_.substr(0, pos);
    if(!resolved)
      nameStr = StringParser::parse(name(), nameStr, varList, resolved);
    if(!resolved)
      throw(Exception("must contain a variable name ('name [= expr]')"));
    Tokenizer token(nameStr);
    token.get();
    if(token.type!=Tokenizer::NAME)
      throw(Exception("must contain a variable name ('name [= expr]')"));
    _name = token.name; // set new name
    token.get();
    if(token.type != Tokenizer::END)
      throw(Exception(token.infoString("unknown token")));

    if(pos == std::string::npos) // only variable name
    {
      status = VALUE;
      value  = 0.;
      return;
    }

    // expression after '='
    this->text = text_.substr(pos+1, std::string::npos);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/

void ExpressionVariable::usedVariables(const VariableList &varList, std::set<std::string> &usedName)
{
  try
  {
    if(status == UNDEFINED)
      throw(Exception("undefined variable"));
    if(status == VALUE)
      return;
    if(status == TEXT)
    {
      Bool resolved;
      std::string text_ = StringParser::parse(name(), this->text, varList, resolved);
      if(!resolved)
        throw(Exception("unresolved variables"));
      expr = Expression::parse(text_);
      status = EXPRESSION;
    }
    if(status == EXPRESSION)
      expr->usedVariables(varList, usedName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/

void ExpressionVariable::simplify(const VariableList &varList)
{
  try
  {
    if(status == TEXT)
    {
      Bool resolved;
      std::string text_ = StringParser::parse(name(), this->text, varList, resolved);
      if(!resolved)
        return;
      expr = Expression::parse(text_);
      status = EXPRESSION;
    }

    if(status == EXPRESSION)
    {
      Bool resolved;
      expr = expr->simplify(varList, resolved);
      if(!resolved)
        return;
      value  = expr->evaluate(varList);
      status = VALUE;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/

ExpressionVariablePtr ExpressionVariable::derivative(const std::string &varName, const VariableList &varList) const
{
  try
  {
    if(status == UNDEFINED)
      throw(Exception("undefined variable"));
    if(status == VALUE)
      return std::make_shared<ExpressionVariable>(name(), 0.0);
    if(status == EXPRESSION)
      return std::make_shared<ExpressionVariable>(name(), expr->derivative(varName));
    // status == TEXT
    Bool resolved;
    std::string text_ = StringParser::parse(name(), this->text, varList, resolved);
    if(!resolved)
      throw(Exception("unresolved variables"));
    return std::make_shared<ExpressionVariable>(name(), Expression::parse(text_)->derivative(varName));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/

Double ExpressionVariable::evaluate(const VariableList &varList) const
{
  try
  {
    if(status == CIRCULAR)
      throw(Exception("circular expression"));
    if(status == VALUE)
      return value;
    if(status == EXPRESSION)
      return expr->evaluate(varList);

    Bool resolved;
    std::string text_ = StringParser::parse(name(), this->text, varList, resolved);
    if(!resolved)
      throw(Exception("unresolved variables"));

    const Status stat = status;
    const_cast<ExpressionVariable*>(this)->status = CIRCULAR;
    Double d = Expression::parse(text_)->evaluate(varList);
    const_cast<ExpressionVariable*>(this)->status = stat;
    return d;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/

std::string ExpressionVariable::getParsedText(const VariableList &varList, Bool &resolved) const
{
  try
  {
    if(status == VALUE)
    {
      std::stringstream ss;
      ss<<std::setprecision(std::numeric_limits<decltype(value)>::digits10)<<value;
      return ss.str();
    }

    if(status == EXPRESSION)
      return expr->string();

    if(status == UNDEFINED)
    {
      resolved = FALSE;
      return "{_undefined_}";
    }

    if(status == CIRCULAR)
      throw(Exception("is circular defined"));

    // status == TEXT
    const_cast<ExpressionVariable*>(this)->status = CIRCULAR;
    auto result = StringParser::parse(name(), text, varList, resolved);
    const_cast<ExpressionVariable*>(this)->status = TEXT;
    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("Expression("+name()+" = '"+getText()+"')", e)
  }
}

/***********************************************/
/***********************************************/

VariableList::VariableList(const VariableList &x) : map(x.map)
{
  for(auto &p : map)
    p.second = std::make_shared<ExpressionVariable>(*p.second); // clone data
}

/***********************************************/

VariableList &VariableList::operator=(const VariableList &x)
{
  map = x.map;
  for(auto &p : map)
    p.second = std::make_shared<ExpressionVariable>(*p.second); // clone data
  return *this;
}

/***********************************************/

VariableList &VariableList::operator+=(const VariableList &x)
{
  for(const auto &var : x)
    addVariable(var.second);
  return *this;
}

/***********************************************/

ExpressionVariablePtr VariableList::operator[](const std::string &name)
{
  try
  {
    auto var = find(name);
    if(!var)
      return addVariable(std::make_shared<ExpressionVariable>(name));
    return var;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ExpressionVariablePtr VariableList::find(const std::string &name) const
{
  try
  {
    auto iter = map.find(name);
    if(iter != map.end())
      return iter->second;
    // Hack to read old constant definition without brackets
    if(name == "pi")  return std::make_shared<ExpressionVariable>("pi",  PI);
    if(name == "rho") return std::make_shared<ExpressionVariable>("rho", RAD2DEG);
    if(name == "nan") return std::make_shared<ExpressionVariable>("nan", NAN_EXPR);
    return ExpressionVariablePtr(nullptr);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ExpressionVariablePtr VariableList::addVariable(ExpressionVariablePtr var)
{
  try
  {
    var = std::make_shared<ExpressionVariable>(*var); // copy
    map[var->name()] = var;
    return var;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addVariable(const std::string &name, VariableList &varList)
{
  addVariable(std::make_shared<ExpressionVariable>(name), varList);
}

/***********************************************/

void addVariable(const std::string &name, Double value, VariableList &varList)
{
  addVariable(std::make_shared<ExpressionVariable>(name, value), varList);
}

/***********************************************/

void addVariable(const std::string &name, const std::string &text, VariableList &varList)
{
  addVariable(std::make_shared<ExpressionVariable>(name, text), varList);
}

/***********************************************/

void addVariable(ExpressionVariablePtr var, VariableList &varList)
{
  varList.addVariable(var);
}

/***********************************************/
