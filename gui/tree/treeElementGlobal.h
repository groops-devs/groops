/***********************************************/
/**
* @file treeElementGlobal.h
*
* @brief The global element.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-07-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTGLOBAL__
#define __GROOPSGUI__TREEELEMENTGLOBAL__

#include "base/importGroops.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementGlobal : public TreeElementComplex
{
  Q_OBJECT

public:
  TreeElementGlobal(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                    XmlNodePtr xmlNode);
  virtual ~TreeElementGlobal() override;

/** @brief Generate XML-tree.
* recursively called for all children. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief calls @a element->newLink(elementInGlobal) for all elements in global. */
void informAboutGlobalElements(TreeElement *element) const;

/** @brief Returns list of global types. */
XsdComplexPtr xsdComplexGlobalTypes() const;

/** @brief Returns list of element names for all child elements in global. */
QStringList getChildrenNames() const;

/** @brief create elements for a choice with @a index. */
void createChildrenElements(int index, XmlNodePtr xmlNode=XmlNodePtr(nullptr));

// ========================================================

// Add/Remove elements
// -------------------
public:
/** @brief Is it possible to insert an element with @a type before @a beforeElement?
* New elements can be added directly to the global section. */
virtual Bool canAddChild(TreeElement *beforeElement, const QString &type) const override;

/** @brief Is it possible to remove this element from tree?
* All elements except the add element in the global section can be deleted. */
virtual Bool canRemoveChild(TreeElement *element) const override;

// ========================================================

/** @brief can @a element be set in the global section? */
Bool canSetGlobal(TreeElement *element) const;

/** @brief Creates a new @a element in the global section.
* Is undoable.
* @return success? */
Bool setGlobal(TreeElement *element);

/** @brief add child before @a beforeElement.
* Is undoable.
* the new element is selected.
* @return success */
virtual Bool addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, bool moved = false) override;
virtual Bool addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, QString &label, bool moved = false);

/** @brief remove child @a element.
* Is undoable.
* if successful the element is deleted.
* @return success */
virtual Bool removeChild(TreeElement *element) override;
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTGLOBAL__ */
