/***********************************************/
/**
* @file stringParser.h
*
* @brief string manipulation parser
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2016-03-26
*
*/
/***********************************************/

#ifndef __GROOPS_STRINGPARSER__
#define __GROOPS_STRINGPARSER__

// Latex documentation
#ifdef DOCSTRING_Parser
static const char *docstringParserText = R"(
\subsection{Text parser}\label{general.parser:text}
Before the mathematical expression parser evaluates the expression, a simple text parser is applied.
The text parser is used for all input fields (also file names). It scans the text for terms like
\verb|{variable}| and replaces it by the text content of the \verb|variable| defined in the global section.

The text parser also evaluates terms in the form \verb|{expression:format}| and replaces it by a formatted
output. The \verb|format| contains the text to be written as output.
It can contain embedded format specifiers that are replaced by the value of the expression
and formatted as requested (also multiple times). In the following, the resulting formatted output is given in the
brackets for an expression with the example value of 57493.8:
\begin{itemize}
\item \verb|%i|: Integer [57494]
\item \verb|%f|: Decimal floating point [57493.800000]
\item \verb|%e|: Scientific notation [5.749380e+04]
\item \verb|%g|: Use the shortest representation: \verb|%e| or \verb|%f| [57493.8]
\item \verb|%c|: Interpret number as ASCII character
\item \verb|%%|: Write a single literal \verb|%| character
\end{itemize}
The following specifiers interpret the value of the expression as MJD (modified Julian date):
\begin{itemize}
\item \verb|%y|: Four digit year [2016]
\item \verb|%Y|: Two digit year [16]
\item \verb|%m|: Month [04]
\item \verb|%d|: Day of month [15]
\item \verb|%H|: Hour [19]
\item \verb|%M|: Minute [12]
\item \verb|%S|: Second [00]
\item \verb|%D|: Date (same as \verb|%y-%m-%d|) [2016-04-15]
\item \verb|%T|: Time (same as \verb|%H-%M-%S|) [19-12-00]
\item \verb|%W|: GPS week [1892]
\item \verb|%w|: Day of GPS week (0..6) [5]
\item \verb|%O|: Day of year (1..366)
\end{itemize}
The format can be specified further with \verb|%[width][.precision]specifier|,
where \verb|[width]| is the minimum number of characters to be printed.
If the value to be printed is shorter than this number, the result is padded with blank spaces
(or zeros if \verb|[width]| starts with a zero).
The \verb|[.precision]| defines the number of digits after the period (for \verb|%g| the number of
significant digits instead).

Example:
Two variables \config{time}=\verb|57493+19/24+12/1440| and \config{satellite}=\verb|swarm| are
set in the global section. The \config{inputfile}=\verb|data/{time:%y}/{satellite}_{time:%D}.dat|
is expanded to \verb|"data/2016/swarm_2016-04-15.dat"|.

Example:
The variable \config{x}=\verb|3+5| is set in the global section.
The expression \config{number}=\verb|2*x| is evaluated by the expression parser to \verb|=16|.
In contrast if we use brackets like in \config{number}=\verb|2*{x}| the expression is first evaluated
by the text parser to \verb|"2*3+5"| and the expression parser now gives the result \verb|=11|.
)";
#endif

/***********************************************/

#include "base/importStd.h"
#include "expressionParser.h"

/***** CLASS ***********************************/

/** @brief string manipulation.
* @ingroup parserGroup */
namespace StringParser
{
  /** @brief string manipulation.
  * {variable} is replaced by the content of the variable list.
  * {expression:format} is replaced by the format string. Before the expression is evaluated and
  * the result is inserted and formatted at every format identifier:
  * - %%: the single % character.
  * - %c: character.
  * - %i: integer.
  * - %f, %e, %g: float.
  * - %y, %m, %d: year, month, day.
  * - %H, %M, %S: hour, minute, second.
  * - %Y: two digit year.
  * - %D: Date: yyyy-mm-dd.
  * - %T: Time: hh-mm-ss.
  * - %O: day of year.
  * - %W, %w: GPS week, day of week.
  *
  * Example {51544.5:%D_T%T} -> '2000-01-01_T12-00-00'. */
  std::string parse(const std::string &name, const std::string &text, const VariableList &varList, Bool &resolved);

  /** @brief string manipulation.
  * Convenience function.
  * An expception is thrown if the string cannot resolved completly. */
  std::string parse(const std::string &text, const VariableList &varList);
}

/***********************************************/

#endif /* __GROOPS__ */
