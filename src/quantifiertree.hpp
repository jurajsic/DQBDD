#ifndef DQBF_BDD_QUANTIFIERTREE_HPP
#define DQBF_BDD_QUANTIFIERTREE_HPP

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
    virtual QuantifierTreeFormula* changeToFormula(Cudd &mgr) = 0;
    virtual void negate() = 0;
    //virtual void renameVar(Variable oldVar, Variable newVar) = 0;
};

class QuantifierTreeFormula : public QuantifierTreeNode, public Formula {
public:
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    void localise();
    QuantifierTreeFormula* changeToFormula(Cudd &mgr);
    void negate();
};

class QuantifierTree : public QuantifierTreeNode {
private:
    std::list<QuantifierTreeNode*> children;
    bool isConj; // TODO change to operator, will have to learn how to do that tho


    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove);

    std::ostream& print(std::ostream& out) const override;

    //void changeChildWithFormula(std::list<QuantifierTreeNode*>::iterator childToChange, Cudd &mgr);
public:
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr);
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManipulator &qvManipulator);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    ~QuantifierTree();

    void localise();

    /**
     * @brief Changes this instance of QuantifierTree to the instance of Formula,
     * it is assumed that this instance was created by new, after calling this function 
     * the pointer to the QuantifierTree is invalidated (do NOT delete it!)
     * 
     * @param mgr The Cudd manager used for creating matrix of Formula
     * @return Pointer to the resulting instance of Formula, needs to be deleted after it was used
     */
    QuantifierTreeFormula* changeToFormula(Cudd &mgr);
    void negate();
    //void renameVar(Variable oldVar, Variable newVar);

    void addChild(QuantifierTreeNode *child);
};


#endif