#ifndef QUANTIFIERTREE_HPP
#define QUANTIFIERTREE_HPP

#include <list>
#include "variable.hpp"
#include "formula.hpp"


class QuantifierTreeNode {
public:    
    virtual VariableSet getSupportSet() = 0;
    virtual void localise() = 0; // TODO ake parametre a return???
    virtual void pushExistVar(Variable var) = 0;
    virtual void pushUnivVar(Variable var) = 0;
    virtual Formula* getFormula(Cudd &mgr) = 0;
    //virtual void renameVar(Variable oldVar, Variable newVar) = 0;
};

class QuantifierTree : public QuantifierTreeNode, public QuantifiedVariablesManipulator {
private:
    std::list<QuantifierTreeNode*> children;
    bool isConj; // TODO change to operator, will have to learn how to do that tho

    VariableSet supportSet;

    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove);
public:
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    ~QuantifierTree();

    void localise();
    void pushExistVar(Variable var);
    void pushUnivVar(Variable var);
    VariableSet getSupportSet();

    Formula* getFormula(Cudd &mgr);
    //void renameVar(Variable oldVar, Variable newVar);

    void addChild(QuantifierTreeNode *child);
};

class QuantifierTreeVariable : public QuantifierTreeNode, public Variable {
    bool isExistential = false;
    bool isUniversal = false;
    QuantifiedVariablesManager *qvMgr;
public:
    QuantifierTreeVariable(int id, Cudd &mgr, QuantifiedVariablesManager &qvMgr);
    VariableSet getSupportSet();
    void localise();
    void pushExistVar(Variable var);
    void pushUnivVar(Variable var);
    Formula* getFormula(Cudd &mgr);
};

#endif