#ifndef QUANTIFIERTREE_HPP
#define QUANTIFIERTREE_HPP

#include <list>
#include "variable.hpp"
#include "formula.hpp"

class QuantifierTreeFormula;

class QuantifierTreeNode : virtual public QuantifiedVariablesManipulator {
public:
    virtual ~QuantifierTreeNode() = default;
    virtual void localise() = 0; // TODO ake parametre a return???
    void pushExistVar(Variable var);
    void pushUnivVar(Variable var);
    virtual QuantifierTreeFormula* getFormula(Cudd &mgr) = 0;
    //virtual void renameVar(Variable oldVar, Variable newVar) = 0;
};

class QuantifierTreeFormula : public QuantifierTreeNode, public Formula {
public:
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    void localise();
    QuantifierTreeFormula* getFormula(Cudd &mgr);
};

class QuantifierTree : public QuantifierTreeNode {
private:
    std::list<QuantifierTreeNode*> children;
    bool isConj; // TODO change to operator, will have to learn how to do that tho

    VariableSet supportSet;

    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove);

    //void changeChildWithFormula(std::list<QuantifierTreeNode*>::iterator childToChange, Cudd &mgr);
public:
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    ~QuantifierTree();

    void localise();
    VariableSet getSupportSet() override;
    std::ostream& print(std::ostream& out) const override;

    QuantifierTreeFormula* getFormula(Cudd &mgr);
    //void renameVar(Variable oldVar, Variable newVar);

    void addChild(QuantifierTreeNode *child);
};


#endif