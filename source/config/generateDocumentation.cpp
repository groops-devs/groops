/***********************************************/
/**
* @file generateDocumentation.cpp
*
* @brief Generates a latex user documentation file.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2008-10-14
*/
/***********************************************/

#define DOCSTRING_Doodson
#define DOCSTRING_GnssType
#define DOCSTRING_Parser

#include "base/import.h"
#include "base/string.h"
#include "base/doodson.h"
#include "base/gnssType.h"
#include "parser/expressionParser.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "inputOutput/file.h"
#include "inputOutput/logging.h"
#include "files/fileFormatRegister.h"
#include "config/configRegister.h"
#include "programs/program.h"
#include "generateDocumentation.h"

/***********************************************/

class _ClassDoodson : public SchemaClass
{
public:
  std::string typeName() const {return "doodson";}
  void registerConfigSchema(Config &/*config*/) const {}
  void generateDocumentation(Documentation &documentation) const {documentation.writeText(docstringDoodson);}
};
static _ClassDoodson _classDoodson;

/***********************************************/

class _ClassGnssType : public SchemaClass
{
public:
  std::string typeName() const {return "gnssType";}
  void registerConfigSchema(Config &/*config*/) const {}
  void generateDocumentation(Documentation &documentation) const {documentation.writeText(docstringGnssType);}
};
static _ClassGnssType _classGnssType;

/***********************************************/

class Documentation::Token
{
public:
  enum Type {TEXT, PARAGRAPH, NEWLINE, ESCAPE, BRACKET, COMMAND, EQUATION, EQUATIONINLINE};
  Type        type;
  std::string text;
  std::vector<std::string> options;

  Token(Type type_, const std::string &text_) : type(type_), text(text_) {}
};

/***********************************************/

std::vector<Documentation::Token> Documentation::tokenize(const std::string &text)
{
  try
  {
    std::vector<Token> tokens;
    if(text.empty())
      return tokens;

    std::string::size_type pos = 0, posOld = 0;
    for(;;)
    {
      pos = text.find_first_of("\n%{$\\", pos);

      // single return -> ignore
      if((pos != std::string::npos) && (text.at(pos) == '\n') && !((pos+1 < text.size()) && (text.at(pos+1) == '\n')))
      {
        pos++;
        continue;
      }

      // save old text before
      if(pos > posOld)
        tokens.push_back(Token(Token::TEXT, String::replaceAll(text.substr(posOld, pos-posOld), {{"~", " "}, {"<", "&lt;"}})));
      if(pos == std::string::npos)
        break;

      if(text.at(pos) == '%')       // comment
      {
        pos = text.find("\n", pos);
      }
      else if(text.at(pos) == '\n') // paragraph
      {
        pos += 2;
        if(tokens.size() && (tokens.back().type != Token::PARAGRAPH))
          tokens.push_back(Token(Token::PARAGRAPH, ""));
      }
      else if(text.at(pos) == '{') // bracket
      {
        auto posTmp = ++pos;
        Int depth  = 1;
        do
        {
          posTmp = text.find_first_of("{}", posTmp);
          if(posTmp == std::string::npos)
            throw(Exception("missing '}'"));
          depth += (text.at(posTmp++) == '{') ? +1 : -1;
        }
        while(depth>0);

        // special case \begin{equation}
        if((text.substr(pos, posTmp-pos-1) == "equation") && (tokens.back().text == "begin"))
        {
          pos = posTmp;
          posTmp = text.find(R"(\end{equation})", pos);
          if(posTmp == std::string::npos)
            throw(Exception("expected '\\end{equation}'"));
          tokens.back() = Token(Token::EQUATION, text.substr(pos, posTmp-pos));
          pos = posTmp + 15;
        }
        else
        {
          tokens.push_back(Token(Token::BRACKET, text.substr(pos, posTmp-pos-1)));
          pos = posTmp;
        }
      }
      else if(text.at(pos) == '$') // EQUATIONINLINE
      {
        auto posTmp = text.find('$', ++pos);
        if(posTmp == std::string::npos)
          throw(Exception("missing '$'"));
        tokens.push_back(Token(Token::EQUATIONINLINE, text.substr(pos, posTmp-pos)));
        pos = posTmp+1;
      }
      else if(text.at(pos) == '\\')
      {
        pos++;
        if(pos >= text.size())
          throw(Exception("text ended unexpectedly after '\\'"));
        if(text.at(pos) == '[') // EQUATION
        {
          auto posTmp = text.find(R"(\])", ++pos);
          if(posTmp == std::string::npos)
            throw(Exception("missing '\\]'"));
          tokens.push_back(Token(Token::EQUATION, text.substr(pos, posTmp-pos)));
          pos = posTmp+2;
        }
        else if(text.at(pos) == '\\')        // NEWLINE
          tokens.push_back(Token(Token::NEWLINE, text.substr(pos++, 1)));
        else if(!std::isalpha(text.at(pos))) // ESCPAE
          tokens.push_back(Token(Token::ESCAPE, text.substr(pos++, 1)));
        else // COMMAND
        {
          posOld = pos;
          while((pos < text.size()) && std::isalpha(text.at(pos))) pos++;
          tokens.push_back(Token(Token::COMMAND, text.substr(posOld, pos-posOld)));
          while((pos < text.size()) && std::isspace(text.at(pos))) pos++;
          if(tokens.back().text == "verb")
          {
            auto separator = text.at(pos++);
            auto posTmp = text.find(separator, pos);
            if(posTmp == std::string::npos)
              throw(Exception("missing end of \\verb"s+separator));
            tokens.push_back(Token(Token::BRACKET, String::replaceAll(text.substr(pos, posTmp-pos), "<", "&lt;")));
            pos = posTmp+1;
          }
        }
      }
      else
        throw(Exception("unexpected character"));

      posOld = pos;
    } // for(;;)

    return tokens;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

class Documentation::TableLine
{
public:
  std::string name, type, annotation;
  Bool        optional, unbounded;
  UInt        depth;
};

/***********************************************/

void Documentation::generateTableLine(XmlNodePtr xmlNode, UInt depth, std::vector<TableLine> &tableLines)
{
  try
  {
    TableLine tableLine;
    tableLine.depth = depth;

    // name
    XmlAttrPtr attr = xmlNode->getAttribute("name");
    if(attr)
      tableLine.name = attr->getText();

    // optional
    tableLine.optional = FALSE;
    if(xmlNode->getAttribute("minOccurs"))
      tableLine.optional = TRUE;

    // unbounded
    tableLine.unbounded = FALSE;
    if(xmlNode->getAttribute("maxOccurs"))
      tableLine.unbounded = TRUE;

    // comment
    XmlNodePtr xmlChild = xmlNode->getChild("xs:annotation");
    if(xmlChild)
      xmlChild = xmlChild->getChild("xs:documentation");
    if(xmlChild)
      tableLine.annotation = xmlChild->getText();

    // type
    attr = xmlNode->getAttribute("type");
    if(attr)
      tableLine.type = attr->getText();

    if(xmlNode->findChild("xs:complexType"))
      xmlNode  = xmlNode->getChild("xs:complexType");
    XmlNodePtr sequence = xmlNode->getChild("xs:sequence");
    XmlNodePtr choice   = xmlNode->getChild("xs:choice");

    // setType
    if(choice)
      tableLine.type = "choice";
    if(sequence)
      tableLine.type = "sequence";

    tableLines.push_back(tableLine);

    // sequence
    if(sequence)
      while(sequence->hasChildren())
        generateTableLine(sequence->getNextChild(), depth+1, tableLines);

    // choice
    if(choice)
      while(choice->hasChildren())
        generateTableLine(choice->getNextChild(), depth+1, tableLines);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Documentation::TableLine> Documentation::generateTable(Config &config)
{
  try
  {
    std::vector<TableLine> tableLines;
    XmlNodePtr rootNode = config.table();
    if(!rootNode->hasChildren())
      return tableLines;
    while(rootNode->hasChildren())
      generateTableLine(rootNode->getNextChild(), 0, tableLines);
    return tableLines;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Documentation::isClassType(const std::string &type)
{
  try
  {
    const std::vector<SchemaClass*> classList = SchemaClass::classList();
    return std::find_if(classList.begin(), classList.end(), [&](auto c) {return c->typeName() == type;}) != classList.end();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/
/***** LaTex ***********************************/
/***********************************************/

class DocumentationLatex : public Documentation
{
public:
  OutFile file;

  DocumentationLatex(const FileName &fileName);
 ~DocumentationLatex() {}

  std::string escaping(const std::string &text);

  void writeText(const std::string &text) override {file<<text<<std::endl;}
  void writeConfigTable(Config &config) override;
};

/***********************************************/

DocumentationLatex::DocumentationLatex(const FileName &fileName)
{
  try
  {
    file.open(fileName);
    file<<"% auto generated by GROOPS"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string DocumentationLatex::escaping(const std::string &text)
{
  return String::replaceAll(text, {{"\\", "\\textbackslash{}"}, {"{", "\\{"}, {"}", "\\}"}, {"<", "$<$"}, {">", "$>$"},
                                   {"&", "\\&"}, {"%", "\\%"}, {"$", "\\$"}, {"#", "\\#"}, {"_", "\\_"},
                                   {"~", "\\textasciitilde{}"}, {"^", "\\textasciicircum{}"}});
}

/***********************************************/

void DocumentationLatex::writeConfigTable(Config &config)
{
  try
  {
    std::vector<TableLine> tableLines = generateTable(config);
    if(!tableLines.size())
      return;

    file<<std::endl;
    file<<R"(\keepXColumns)"<<std::endl;
    file<<R"(\begin{tabularx}{\textwidth}{N T A})"<<std::endl;
    file<<R"(\hline)"<<std::endl;
    file<<R"(Name & Type & Annotation\\)"<<std::endl;
    file<<R"(\hline)"<<std::endl;
    for(UInt i=0; i<tableLines.size(); i++)
    {
      std::string type = tableLines.at(i).type;
      if(type == "gnssType")
        type = R"(\hyperref[gnssType]{gnssType})"; // set link
      else if(isClassType(type))
        type = R"(\hyperref[)"+type+"]{"+type.substr(0, type.rfind("Type"))+"}"; // set link

      file<<R"(\hfuzz=500pt)";
      for(UInt k=1; k<tableLines.at(i).depth; k++)
        file<<R"(\quad)";
      if(tableLines.at(i).depth > 0)
        file<<R"(\includegraphics[width=1em]{connector.pdf})";
      if(tableLines.at(i).optional && tableLines.at(i).unbounded)        file<<R"(\includegraphics[width=1em]{element-unbounded.pdf}~)";
      else if(tableLines.at(i).optional && !tableLines.at(i).unbounded)  file<<R"(\includegraphics[width=1em]{element.pdf}~)";
      else if(!tableLines.at(i).optional && tableLines.at(i).unbounded)  file<<R"(\includegraphics[width=1em]{element-mustset-unbounded.pdf}~)";
      else if(!tableLines.at(i).optional && !tableLines.at(i).unbounded) file<<R"(\includegraphics[width=1em]{element-mustset.pdf}~)";
      file<<escaping(tableLines.at(i).name)<<R"( & \hfuzz=500pt )"<<type<<R"( & \hfuzz=500pt )"<<escaping(tableLines.at(i).annotation)<<R"(\\)"<<std::endl;
    }
    file<<R"(\hline)"<<std::endl;
    file<<R"(\end{tabularx})"<<std::endl;
    file<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Documentation::writeLatex(const FileName &directoryName)
{
  try
  {
    // section parser
    // --------------
    {
      DocumentationLatex docu(directoryName.append("general.parser.auto.tex"));
      docu.writeText(docstringParserExpression);
      docu.writeText(docstringParserText);
      docu.writeText(docstringParserDataVariables);
    }

    // File description
    // ----------------
    {
      DocumentationLatex docu(directoryName.append("general.fileFormat.auto.tex"));
      std::vector<FileFormat*> fileFormatList = FileFormat::list();
      FileFormat::sort(fileFormatList);
      for(UInt i=0; i<fileFormatList.size(); i++)
      {
        docu.file<<"\\subsection{"<<fileFormatList.at(i)->title()<<"}\\label{general.fileFormat:"<<fileFormatList.at(i)->type()<<"}";
        docu.writeText(fileFormatList.at(i)->documentation());
        docu.file<<"\n%==================================\n";
      }
    }

    // Program description
    // -------------------
    {
      DocumentationLatex docu(directoryName.append("programs.auto.tex"));
      std::vector<Program::Program*> programList = Program::Program::programList();
      Program::Program::sortList(programList);
      Program::Tags tag = Program::NONE;
      for(UInt i=0; i<programList.size(); i++)
      {
        if(programList.at(i)->tags().at(0) != tag)
        {
          tag = programList.at(i)->tags().at(0);
          docu.file<<"\\section{Programs: "<<Program::tagStrings[tag]<<"}"<<std::endl;
        }

        Config config;
        programList.at(i)->run(config, Parallel::selfCommunicator());
        docu.file<<"\\subsection{"<<programList.at(i)->name()<<"}\\label{"<<programList.at(i)->name()<<"}";
        docu.writeText(programList.at(i)->documentation());
        docu.writeConfigTable(config);
        if(!programList.at(i)->isSingleProcess())
          docu.file<<"This program is \\reference{parallelized}{general.parallelization}."<<std::endl;
        docu.file<<"\\clearpage\n%==================================\n";
      }
    }

    // Classes description
    // -------------------
    {
      DocumentationLatex docu(directoryName.append("classes.auto.tex"));
      std::vector<SchemaClass*> classList = SchemaClass::classList();
      SchemaClass::sort(classList);
      for(UInt i=0; i<classList.size(); i++)
      {
        classList.at(i)->generateDocumentation(docu);
        docu.file<<"\\clearpage\n%==================================\n";
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** HTML ************************************/
/***********************************************/

class DocumentationHtml : public Documentation
{
public:
  FileName          directoryName;
  std::stringstream ss;
  std::stringstream ssSearchDisplayText; //  display text for search results
  std::stringstream ssSearchTokens; // description text tokens
  std::stringstream ssSearchMiscTokens; // config table, tags, ...

  DocumentationHtml(const FileName &directoryName_) : directoryName(directoryName_) {}
 ~DocumentationHtml() {}

  void writeText(const std::string &text) override;
  void writeConfigTable(Config &config) override;
};

/***********************************************/

void DocumentationHtml::writeText(const std::string &text)
{
  try
  {
    std::vector<Token> token = tokenize(text);
    std::string searchTokens;

    Bool isItemOpen = FALSE;
    for(UInt i=0; i<token.size(); i++)
    {
      // lambda function
      // --------------
      auto reference = [](const std::string &text) -> std::string
      {
        auto pos = text.find(":");
        if(pos != std::string::npos)
          return text.substr(0, pos)+".html#"+text.substr(pos+1);
        return text+".html";
      };
      // --------------

      // lambda function
      // --------------
      auto checkBrackets = [&](UInt count=1)
      {
        for(UInt k=0; k<count; k++)
          if(token.at(++i).type != Token::BRACKET)
            throw(Exception("expected '{'"));
      };
      // --------------

      // lambda function
      // --------------
      auto checkLabel = [&]() -> std::string
      {
        if((token.at(i+1).type != Token::COMMAND) || (token.at(i+1).text != "label"))
          return "";
        i++;
        checkBrackets();
        std::string refFile, refAnchor = token.at(i).text;
        auto pos = token.at(i).text.find(":");
        if(pos != std::string::npos)
        {
          refFile   = token.at(i).text.substr(0, pos);
          refAnchor = token.at(i).text.substr(pos+1);
        }
        return R"( id=")"+refAnchor+R"(")";
      };
      // --------------

      // lambda function
      // --------------
      auto addToSearchTokens = [&](const std::string &token)
      {
        searchTokens += (token+" ");
      };

      std::string html;
      Bool addToSearchDisplayText = TRUE;

      const std::string &text = token.at(i).text;
      switch(token.at(i).type)
      {
        case Token::TEXT:           html = text; addToSearchTokens(text); break;
        case Token::EQUATIONINLINE: html = "$"+text+"$"; break;
        case Token::EQUATION:       html = "\\["+ text +"\\]"; break;
        case Token::ESCAPE:         html = text; break;
        case Token::PARAGRAPH:      html = "</p><p>"; break;
        case Token::NEWLINE:
        {
          logWarning<<"Ignored newline '\\\\'"<<Log::endl;
          break;
        }
        case Token::BRACKET:
        {
          html = R"(<strong style="color: red;">{)" + text + "}</strong>";
          logWarning<<"Unexpected {"<<token.at(i).text<<"}"<<Log::endl;
          break;
        }
        case Token::COMMAND:
        {
          const std::string &cmd = token.at(i).text;
          if(cmd == "ldots")            {ss<<"&hellip; ";}
          else if(cmd == "section")     {checkBrackets();  std::string title = token.at(i).text; html = R"(<h1)" + checkLabel() + ">" + title + R"(</h1><p>)"; addToSearchDisplayText = FALSE; }
          else if(cmd == "subsection")  {checkBrackets();  std::string title = token.at(i).text; html = R"(<h2)" + checkLabel() + ">" + title + R"(</h2><p>)";}
          else if(cmd == "emph")        {checkBrackets();  html = R"(<em>)" + token.at(i).text + R"(</em>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "textbf")      {checkBrackets();  html = R"(<b>)" + token.at(i).text + R"(</b>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "verb")        {checkBrackets();  html = R"(<code>)" + token.at(i).text + R"(</code>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "url")         {checkBrackets();  html = R"(<a href=")" + token.at(i).text + R"(" target="_blank">)" + token.at(i).text + R"(</a>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "href")        {checkBrackets(2); html = R"(<a href=")" + token.at(i-1).text + R"(">)" + token.at(i).text + R"(</a>)";}
          else if(cmd == "eqref")       {checkBrackets();  html = R"(\eqref{)" + token.at(i).text + R"(})";}
          else if(cmd == "program")     {checkBrackets();  html = R"(<a class="groops-program" href=")" + reference(token.at(i).text) + R"(">)" + token.at(i).text + R"(</a>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "config")      {checkBrackets();  html = R"(<strong class="groops-config-element">)" + token.at(i).text + R"(</strong>)"; addToSearchTokens(token.at(i).text); }
          else if(cmd == "configClass") {checkBrackets(2); html = R"(<a class="groops-class" href=")" + reference(token.at(i).text) + R"(">)" + token.at(i-1).text + R"(</a>)"; addToSearchTokens(token.at(i-i).text);}
          else if(cmd == "configFile")  {checkBrackets(2); html = R"(<a class="groops-class" href=")" + reference("fileFormat_"+token.at(i).text) + R"(">)" + token.at(i-1).text + R"(</a>)"; addToSearchTokens(token.at(i-i).text);}
          else if(cmd == "file")        {checkBrackets(2); html = R"(<a class="groops-file" href=")" + reference("fileFormat_"+token.at(i).text) + R"(">)" + token.at(i-1).text + R"(</a>)"; addToSearchTokens(token.at(i-1).text); }
          else if(cmd == "reference")   {checkBrackets(2); html = R"(<a class="groops-ref" href=")" + reference(token.at(i).text) + R"(">)" + token.at(i-1).text + R"(</a>)"; addToSearchTokens(token.at(i-1).text); }
          else if(cmd == "ref")         {checkBrackets();  html = R"(<a href=")" + reference(token.at(i).text) + R"(">)" + token.at(i).text + R"(</a>)";}
          else if(cmd == "fig")
          {
            addToSearchDisplayText = FALSE;
            checkBrackets(5);
            html = R"(<figure><img class="figure" style="width:)" + (100*String::toDouble(token.at(i-3).text))%"%i"s + R"(%;" )"
                 + R"(src="../figures/)" + token.at(i-2).text + R"(.png" alt=")" + token.at(i-2).text + R"(">)"
                 + R"(<figcaption class="center">Figure: )" + token.at(i).text + R"(</figcaption></figure>)";
          }
          else if(cmd == "input")
          {
            checkBrackets();
            InFile infile(directoryName.append("..").append("latex").append(token.at(i).text+".tex"));
            std::string text((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
            writeText(text);
          }
          else if(cmd == "label")
          {
            checkBrackets();
            html = R"(<a id=")" + token.at(i).text + R"("></a>)";
            logWarning<<"Label at unexpected position: \\label{"<<token.at(i).text<<"} -> ignored"<<Log::endl;
          }
          else if(cmd == "begin")
          {
            checkBrackets();
            if(token.at(i).text == "itemize")
             html = "<ul>\n";
            else if(token.at(i).text == "enumerate")
              html = "<ol>\n";
            else if(token.at(i).text == "verbatim")
              html = "<pre>";
            else if(token.at(i).text == "equation")
              html = R"(\begin{)" + token.at(i).text + "}";
            else
            {
              html = R"(\begin{)" + token.at(i).text + "}";
              logWarning<<"Unknown Latex command: \\begin{"<<token.at(i).text<<"}"<<Log::endl;
            }
          }
          else if(cmd == "end")
          {
            checkBrackets();
            if(token.at(i).text == "itemize")
            {
              if(isItemOpen)
                html = "</li></ul>\n";
              else
                html = "</ul>";
              isItemOpen = FALSE;
            }
            else if(token.at(i).text == "enumerate")
            {
              if(isItemOpen)
                html = "</li></ol>\n";
              else
                html = "</ol>";
              isItemOpen = FALSE;
            }
            else if(token.at(i).text == "verbatim")
              html = "</pre>";
            else if(token.at(i).text == "equation")
              html = R"(\end{)" + token.at(i).text + "}";
            else
            {
              html = R"(\end{)" + token.at(i).text + "}";
              logWarning<<"Unknown Latex command: \\end{"<<token.at(i).text<<"}"<<Log::endl;
            }
          } // "end"
          else if(cmd == "item")
          {
            if(isItemOpen)
              html = "</li><li>\n";
            else
              html = "<li>";
            isItemOpen = TRUE;
          } // "item"
          else
          {
            html = R"(<strong style="color: red;">\)" + cmd + "</strong>";
            logWarning<<"Unknown Latex command: \\"<<cmd<<Log::endl;
          }
          break;
        }
      }
      ss<<html;
      if(addToSearchDisplayText) ssSearchDisplayText<<html;
    } // for(i)
    ss<<R"(</p>)"<<std::endl;
    ssSearchTokens<<searchTokens;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("text='"+text+"'", e)
  }
}

/***********************************************/

void DocumentationHtml::writeConfigTable(Config &config)
{
  try
  {
    std::vector<TableLine> tableLines = generateTable(config);
    if(!tableLines.size())
      return;

    ss<<R"(<table class="table table-hover">)"<<std::endl;
    ss<<R"(<tr class="table-primary"><th>Name</th><th>Type</th><th>Annotation</th></tr>)"<<std::endl;
    for(UInt i=0; i<tableLines.size(); i++)
    {
      if(i%2)
        ss<<R"(<tr class="">)";
      else
        ss<<R"(<tr class="table-light">)";
      ss<<R"(<td class="m-0"><div class="h-100 config-tree depth-)"<<tableLines.at(i).depth<<R"("><div class="h-100 config )";
      if(tableLines.at(i).optional && tableLines.at(i).unbounded)        ss<<R"(optional-unbounded)";
      else if(tableLines.at(i).optional && !tableLines.at(i).unbounded)  ss<<R"(optional)";
      else if(!tableLines.at(i).optional && tableLines.at(i).unbounded)  ss<<R"(mustset-unbounded)";
      else if(!tableLines.at(i).optional && !tableLines.at(i).unbounded) ss<<R"(mustset)";
      ss<<R"(">)"<<tableLines.at(i).name<<R"(</div></div></td>)";
      if(tableLines.at(i).type == "gnssType")
        ss<<R"(<td><a href="gnssType.html">gnssType</a></td>)";
      else if(isClassType(tableLines.at(i).type))
        ss<<R"(<td><a href=")"<<tableLines.at(i).type<<R"(.html">)"<<tableLines.at(i).type.substr(0,tableLines.at(i).type.rfind("Type"))<<R"(</a></td>)";
      else
        ss<<R"(<td>)"<<tableLines.at(i).type<<R"(</td>)";
      ss<<R"(<td>)"<<tableLines.at(i).annotation<<R"(</td>)";
      ss<<R"(</tr>)"<<std::endl;
      ssSearchMiscTokens<<tableLines.at(i).name<<" "<<tableLines.at(i).type<<" "<<tableLines.at(i).annotation<<" ";
    }
    ss<<R"(</table>)"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Documentation::writeHtml(const FileName &directoryName)
{
  try
  {
    // Load Template
    // -------------
    InFile templateFile(directoryName.append("static").append("template.html"));
    std::string templateHtml((std::istreambuf_iterator<char>(templateFile)), std::istreambuf_iterator<char>());

    OutFile documentationSearchIndex(directoryName.append("documentationSearchIndex.js"));
    documentationSearchIndex<<"var documents = {"<<std::endl;

    // lambda function
    // --------------
    auto writeHtmlFile = [&](const FileName &fileName, const std::string &name, const std::string &content)
    {
      OutFile htmlFile(directoryName.append(fileName));
      std::string page = templateHtml;
      const std::string tagTitle   = "[[[title]]]";
      const std::string tagContent = "[[[content]]]";
      page.replace(page.find(tagTitle),   tagTitle.length(),   name);
      page.replace(page.find(tagContent), tagContent.length(), content);
      htmlFile<<page;
    };
    // --------------

    // lambda function
    // --------------
    auto writeSearchIndex = [&](const std::string &fileName, const std::string &name, const std::string &search, const std::string &table, const std::string &displayText)
    {
      documentationSearchIndex<<"'"<<fileName<<"': { 'name': '"<<name<<"', 'key': '"<<fileName;
      documentationSearchIndex<<"', 'description': '" <<String::trim(String::replaceAll(search,      {{"\n", " "}, {"\\", "\\\\"}, {"\'", "\\'"}}));
      documentationSearchIndex<<"', 'config_table': '"<<String::trim(String::replaceAll(table,       {{"\n", " "}, {"\\", "\\\\"}, {"\'", "\\'"}}));
      documentationSearchIndex<<"', 'display_text': '"<<String::trim(String::replaceAll(displayText, {{"\n", " "}, {"\\", "\\\\"}, {"\'", "\\'"}}));
      documentationSearchIndex<<"'},"<<std::endl;
    };
    // --------------

    // lambda function
    // --------------
    auto writeLatexFile = [&](const std::string &fileName, const std::string &title)
    {
      InFile infile(directoryName.append("..").append("latex").append(fileName+".tex"));
      std::string text((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
      DocumentationHtml doc(directoryName);
      doc.writeText(text);
      writeHtmlFile(fileName+".html", title, doc.ss.str());
      writeSearchIndex(fileName, fileName, doc.ssSearchTokens.str(), doc.ssSearchMiscTokens.str(), doc.ssSearchDisplayText.str());
    };
    // --------------

    std::vector<Program::Program*> programList = Program::Program::programList();
    Program::Program::sortList(programList);

    std::vector<SchemaClass*> classList = SchemaClass::classList();
    SchemaClass::sort(classList);

    std::vector<FileFormat*> fileFormatList = FileFormat::list();
    FileFormat::sort(fileFormatList);

    // Index
    // -----
    {
      std::stringstream content;
      content<<R"(<h1>GROOPS Documentation and Cookbook</h1>)"<<std::endl;
      content<<R"(<p>GROOPS is a software toolkit for gravity field recovery, GNSS processing, and statistical analysis of time series and spatial data.</p>)"<<std::endl;
      content<<R"(<div class="container">)"<<std::endl;
      content<<R"(<div class="row">
                    <div class="col-sm top-buffer">
                          <h5 class="card-title">General Information</h5>
                            <p class="card-text">File formats and basic concepts for using GROOPS and its GUI.</p>
                          <ul class="list-group list-group-flush">
                            <li><a href="general.configFiles.html">Config files</a></li>
                            <li><a href="general.parser.html">Parser</a></li>
                            <li><a href="general.loopsAndConditions.html">Loops and conditions</a></li>
                            <li><a href="general.constants.html">Constants and the setting file</a></li>
                            <li><a href="general.parallelization.html">Parallelization</a></li>
                            <li><a href="general.gui.html">Graphical User Interface (GUI)</a></li>
                            <li><a href="general.fileFormat.html">File formats</a></li>
                          </ul>
                    </div>
                    <div class="col-sm top-buffer">
                          <h5 class="card-title">Mathematic Fundamentals</h5>
                            <p class="card-text">Description of methods used throughout GROOPS.</p>
                          <ul class="list-group list-group-flush">
                            <li><a href="fundamentals.robustLeastSquares.html">Robust least squares adjustment</a></li>
                            <li><a href="fundamentals.basisSplines.html">Basis splines</a></li>
                            <li><a href="fundamentals.autoregressiveModel.html">Autoregressive models</a></li>
                          </ul>
                    </div>
                  </div>)"<<std::endl;
      content<<R"(<div class="row">
                    <div class="col-sm top-buffer">
                          <h5 class="card-title">Cookbook</h5>
                            <p class="card-text">Recipes to explore the GROOPS feature set.</p>
                          <ul class="list-group list-group-flush">
                            <li><a href="cookbook.instrument.html">Instrument data handling</a></li>
                            <li><a href="cookbook.gnssNetwork.html">GNSS satellite orbit determination and station network analysis</a></li>
                            <li><a href="cookbook.gnssPpp.html">GNSS precise point positioning (PPP)</a></li>
                            <li><a href="cookbook.kinematicOrbit.html">Kinematic orbit determination of LEO satellites</a></li>
                            <li><a href="cookbook.gravityFieldPod.html">Gravity field determination from POD data</a></li>
                            <li><a href="cookbook.gravityFieldGrace.html">GRACE gravity field recovery</a></li>
                            <li><a href="cookbook.regionalGeoid.html">Regional geoid determination</a></li>
                          </ul>
                    </div>
                    <div class="col-sm top-buffer">
                          <h5 class="card-title">Program/Class Reference</h5>
                          <p class="card-text">Reference and documentation for all programs and classes in GROOPS.</p>
                          <ul class="list-group list-group-flush">
                            <li><a href="programType.html">List of programs</a></li>
                            <li><a href="classes.html">List of classes</a></li>
                          </ul>
                    </div>
                  </div>)"<<std::endl;
      content<<R"(</div>)"<<std::endl;
      writeHtmlFile("index.html", "Overview", content.str());
    }

    // Translate Latex files
    // ---------------------
    writeLatexFile("general.configFiles",              "Config files");
    writeLatexFile("general.parser",                   "Parser");
    writeLatexFile("general.loopsAndConditions",       "Loops and conditions");
    writeLatexFile("general.gui",                      "Graphical User Interface (GUI)");
    writeLatexFile("general.constants",                "Constants and the setting file");
    writeLatexFile("general.parallelization",          "Parallelization");
    writeLatexFile("fundamentals.robustLeastSquares",  "Robust least squares adjustment");
    writeLatexFile("fundamentals.basisSplines",        "Basis splines");
    writeLatexFile("fundamentals.autoregressiveModel", "Autoregressive model");
    writeLatexFile("cookbook.instrument",              "Instrument data handling");
    writeLatexFile("cookbook.gnssNetwork",             "GNSS satellite orbit determination and station network analysis");
    writeLatexFile("cookbook.gnssPpp",                 "GNSS precise point positioning (PPP)");
    writeLatexFile("cookbook.kinematicOrbit",          "Kinematic orbit determination of LEO satellites");
    writeLatexFile("cookbook.gravityFieldPod",         "Gravity field determination from POD data");
    writeLatexFile("cookbook.gravityFieldGrace",       "GRACE gravity field recovery");
    writeLatexFile("cookbook.regionalGeoid",           "Regional geoid determination");

    // Search page
    // -----------
    {
      std::stringstream content;
      content<<R"(<h1>Search</h1>)"<<std::endl;
      content<<R"(<p id="search-meta-info"></p>)"<<std::endl;
      content<<R"(<div id="searchResults"></div>)"<<std::endl;
      content<<R"(<script type="text/javascript" src="documentationSearchIndex.js"></script>)"<<std::endl;
      content<<R"(<script type="text/javascript" src="static/searchtools.js"></script>)"<<std::endl;
      content<<R"(<script>window.onload = startSearch();</script>)"<<std::endl;
      writeHtmlFile("search.html", "Search", content.str());
    }

    // Progams overview page
    // ---------------------
    {
      std::stringstream content;
      content<<R"(<h1>Programs</h1>)"<<std::endl;
      content<<R"(<p>This reference manual details programs included in GROOPS, describing what they are and what they do.
      For usage examples see the cookbook in the <a href="index.html">documentation overview</a>.</p>)"<<std::endl;
      Program::Tags tag = Program::NONE;
      for(UInt i=0; i<programList.size(); i++)
      {
        if(programList.at(i)->tags().at(0) != tag)
        {
          if(tag != Program::NONE)
            content<<R"(</ul></p>)"<<std::endl;
          tag = programList.at(i)->tags().at(0);
          content<<R"(<h2>)"<<Program::tagStrings[tag]<<R"(</h2>)"<<std::endl;
          content<<R"(<p><ul>)"<<std::endl;
        }
        std::string name = programList.at(i)->name();
        content<<R"(<li><a href=")"<<name<<R"(.html">)"<<name<<R"(</a></li>)"<<std::endl;
      }
      content<<R"(</ul></p>)"<<std::endl;
      writeHtmlFile("programType.html", "Programs", content.str());
    }

    // Program pages
    // -------------
    for(UInt i=0; i<programList.size(); i++)
    {
      const std::string name = programList.at(i)->name();
      DocumentationHtml doc(directoryName);
      doc.ss<<R"(<h1>)"<<name<<R"(</h1><p>)"<<std::endl;
      doc.writeText(programList.at(i)->documentation());
      Config config;
      programList.at(i)->run(config, Parallel::selfCommunicator());
      doc.writeConfigTable(config);
      if(!programList.at(i)->isSingleProcess())
        doc.ss<<R"(This program is <a class="groops-ref" href="general.parallelization.html">parallelized</a>.)"<<std::endl;
      writeHtmlFile(name+".html", name, doc.ss.str());
      writeSearchIndex(name, name, doc.ssSearchTokens.str(), doc.ssSearchMiscTokens.str(), doc.ssSearchDisplayText.str());
    }

    // Classes overview page
    // ---------------------
    {
      std::stringstream content;
      content<<R"(<h1>Classes</h1>)"<<std::endl;
      content<<R"(<p>This reference manual details classes included in GROOPS, describing what they are and what they do.
      For usage examples see the cookbook in the <a href="index.html">documentation overview</a>.</p>)"<<std::endl;
      content<<R"(<p><ul>)"<<std::endl;
      for(UInt i=0; i<classList.size(); i++)
      {
        std::string type = classList.at(i)->typeName();
        if(type == "gnssType")
          content<<R"(<li><a href="gnssType.html">gnssType</a></li>)"<<std::endl;
        else
          content<<R"(<li><a href=")"<<type<<R"(.html">)"<<type.substr(0, type.rfind("Type"))<<R"(</a></li>)"<<std::endl;
      }
      content<<R"(</ul></p>)"<<std::endl;
      writeHtmlFile("classes.html", "Classes", content.str());
    }

    // Class pages
    // -----------
    for(UInt i=0; i<classList.size(); i++)
    {
      std::string pageContent = templateHtml;
      const std::string name = classList.at(i)->typeName();
      DocumentationHtml doc(directoryName);
      classList.at(i)->generateDocumentation(doc);
      writeHtmlFile(name+".html", name, doc.ss.str());
      writeSearchIndex(name, name/*.substr(0, name.find("Type"))*/, doc.ssSearchTokens.str(), doc.ssSearchMiscTokens.str(), doc.ssSearchDisplayText.str());
    }

    // File format overview page
    // -------------------------
    {
      const std::string fileName = "general.fileFormat";
      const std::string title    = "File Formats";
      InFile infile(directoryName.append("..").append("latex").append(fileName+".tex"));
      std::string text((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
      DocumentationHtml doc(directoryName);
      doc.writeText(text);

      std::stringstream content;
      content<<R"(<p><ul>)"<<std::endl;
      for(UInt i=0; i<fileFormatList.size(); i++)
      {
        const std::string name = fileFormatList.at(i)->title();
        content<<R"(<li><a href="fileFormat_)"<<fileFormatList.at(i)->type()<<R"(.html">)"<<name<<R"(</a></li>)"<<std::endl;
      }
      content<<R"(</ul></p>)"<<std::endl;

      writeHtmlFile(fileName+".html", title, doc.ss.str() + content.str());
      writeSearchIndex(fileName, fileName, doc.ssSearchTokens.str(), doc.ssSearchMiscTokens.str(), doc.ssSearchDisplayText.str());
    }

    // File format pages
    // -----------------
    for(UInt i=0; i<fileFormatList.size(); i++)
    {
      const std::string name = fileFormatList.at(i)->title();
      DocumentationHtml doc(directoryName);
      doc.ss<<R"(<h1>)"<<name<<R"( (file format)</h1><p>)"<<std::endl;
      doc.writeText(fileFormatList.at(i)->documentation());
      writeHtmlFile("fileFormat_"+fileFormatList.at(i)->type()+".html", name, doc.ss.str());
      writeSearchIndex("fileFormat_"+fileFormatList.at(i)->type(), name, doc.ssSearchTokens.str(), doc.ssSearchMiscTokens.str(), doc.ssSearchDisplayText.str());
    }

    documentationSearchIndex<<"};"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Documentation::write(const FileName &directoryName)
{
  try
  {
    writeLatex(directoryName.append("latex"));
    writeHtml(directoryName.append("html"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
