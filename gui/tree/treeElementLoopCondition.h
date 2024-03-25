/***********************************************/
/**
* @file treeElementLoopCondition.h
*
* @brief Loop or condition of an element.
*
* @author Torsten Mayer-Guerr
* @date 2024-01-22
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTLOOPCONDITION__
#define __GROOPSGUI__TREEELEMENTLOOPCONDITION__

#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElementChoice.h"

/***** CLASS ***********************************/

class TreeElementLoopCondition : public TreeElementChoice
{
  Q_OBJECT

  QString _name, _label;

public:
  TreeElementLoopCondition(Tree *tree, TreeElement *affectedElement, XsdElementPtr xsdElement, XmlNodePtr xmlNode)
    : TreeElementChoice(tree, nullptr, xsdElement, "", xmlNode, !xmlNode, false), affectedElement(affectedElement)
    {
      _name  = (xsdElement == tree->xsdElementLoop) ? "loop" : "condition";
      _label = (xsdElement == tree->xsdElementLoop) ? "_localLoop_" : "_localCondition_";
    }

  TreeElement *affectedElement;

  QString name()         const override {return _name;}
  QString label()        const override {return _label;}
  bool optional()        const override {return true;}
  bool unbounded()       const override {return false;}
  bool canSetLoop()      const override {return false;}
  bool canSetCondition() const override {return false;}
  bool canDisabled()     const override {return false;}
  bool canRename()       const override {return false;}
};

/***********************************************/

#endif
