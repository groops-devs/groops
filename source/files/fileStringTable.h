/***********************************************/
/**
* @file fileStringTable.h
*
* @brief Read/write table of strings.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-03
*
*/
/***********************************************/

#ifndef __GROOPS_FILESTRINGTABLE__
#define __GROOPS_FILESTRINGTABLE__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_StringList
static const char *docstringStringList = R"(
White space separated list of strings.
Comments are allowed and all the text from the character '\verb|#|' to the end of the line is ignored.
Strings containing white spaces or the '\verb|#|' character must be set in quotes ('\verb|""|').

\begin{verbatim}
# IGSR3 stations
abmf
abpo
ade1
adis
ajac
albh
algo
alic
alrt
amc2
aoml
areq
arev
artu
asc1
\end{verbatim}
)";
#endif

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_StringTable
static const char *docstringStringTable = R"(
White space separated table of strings in row and columns.
Additional columns in a row may represent alternatives, e.g. for core stations in a GNSS network.
Comments are allowed and all the text from the character '\verb|#|' to the end of the line is ignored.
Strings containing white spaces or the '\verb|#|' character must be set in quotes  ('\verb|""|').

\begin{verbatim}
# core network with alternative stations
artu mdvj mdvo nril
asc1 sthl
bahr bhr1 yibl nama
chat chti auck
chpi braz ufpr savo
ckis nium
coco xmis dgar dgav
cro1 scub abmf lmmf aoml
daej taej suwn osn1
darw kat1 tow2 alic
dav1 maw1
drao albh will holb nano
fair whit
glps guat
gode godz usno usn3
goug
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "inputOutput/fileName.h"

/** @addtogroup filesGroup */
/// @{

/***** FUNCTIONS *******************************/

/** @brief Write into a stringList file. */
inline void writeFileStringList(const FileName &fileName, const std::vector<std::string> &x);

/** @brief Write into a stringTable file. */
void writeFileStringTable(const FileName &fileName, const std::vector<std::vector<std::string>> &x);

/** @brief Read from a stringList file. */
inline void readFileStringList(const FileName &fileName, std::vector<std::string> &x);

/** @brief Read from a stringTable file. */
void readFileStringTable(const FileName &fileName, std::vector<std::vector<std::string>> &x);

/// @}

/***********************************************/
/***** INLINES ********************************/
/***********************************************/

inline void writeFileStringList(const FileName &fileName, const std::vector<std::string> &x)
{
  try
  {
    std::vector<std::vector<std::string>> table(x.size());
    for(UInt i=0; i<x.size(); i++)
      table.at(i) = {x.at(i)};
    writeFileStringTable(fileName, table);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void readFileStringList(const FileName &fileName, std::vector<std::string> &x)
{
  try
  {
    std::vector<std::vector<std::string>> table;
    readFileStringTable(fileName, table);
    for(const auto &line : table)
      x.insert(x.end(), line.begin(), line.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
