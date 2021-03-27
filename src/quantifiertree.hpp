/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DQBDD_QUANTIFIERTREE_HPP
#define DQBDD_QUANTIFIERTREE_HPP

#include <list>
#include <ostream>

#include "DQBDDvariable.hpp"
#include "DQBDDformula.hpp"

class QuantifierTreeFormula;

/**
 * @brief Base class for nodes in quantifier trees
 */
class QuantifierTreeNode : virtual public QuantifiedVariablesManipulator {
protected:
    /**
     * @brief the set of universal variables which are outside this subformula
     * or in dependency set of existential variable which is outside this subformula
     */
    VariableSet uVarsOutsideThisSubtree = { };
    // the set of universal variables in support set or in dependecy set of ex. var in support set
    VariableSet uVarsSupportSet = {};

public:
    QuantifierTreeNode(QuantifiedVariablesManager &qvmgr);
    virtual ~QuantifierTreeNode() = default;
    /**
     * @brief pushes quantifiers inside the subtree rooted in this node
     * 
     * @param uVarsOutsideThisSubtree - universal variables which are in uVarsSupportSet
     * of nodes outside this subtree (but only those which are in this node's
     *  uVarsSupportSet); used during localisation for disjunction
     */
    virtual void localise(const VariableSet &uVarsOutsideThisSubtree) = 0;
    // pushes the existential variable var into this node
    void pushExistVar(Variable var);
    // pushes the universal variable var into this node
    void pushUnivVar(Variable var);
    /**
     * @brief Changes this instance of node into a formula
     */
    virtual QuantifierTreeFormula* changeToFormula(Cudd &mgr) = 0;
    // negates this node
    virtual void negate() = 0;

    virtual VariableSet const &getUVarsSupportSet();
};

/**
 * @brief A node that is also a formula (used for terminal nodes)
 */
class QuantifierTreeFormula : public QuantifierTreeNode, public Formula {
private:
    void computeSupportSets();
public:
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    void localise(const VariableSet &uVarsOutsideThisSubtree) override;
    QuantifierTreeFormula* changeToFormula(Cudd &) override;
    void negate() override;

    VariableSet const &getSupportSet() override;
    VariableSet const &getUVarsSupportSet() override;
    //void addToUVarsOutsideThisSubtree(const VariableSet &varsToAdd) override;
};

/**
 * @brief A node that represents a root of a quantifier subtree
 */
class QuantifierTree : public QuantifierTreeNode {
private:
    // the children of root
    std::list<QuantifierTreeNode*> children;

    // if isConj=true, the root has assigned conjuction, otherwise disjunction
    bool isConj;

    /**
     * @brief Adds another child to the list of children (assuming *this is a root)
     * 
     * If child has the same operator as this quantifier tree and does not contain a quantifier prefix,
     * the children of this child are added instead, while child is deleted.
     * 
     * @param child - the root of subtree of the child to add
     */
    void addChild(QuantifierTreeNode *child);

    /**
     * @brief Removes from one ordered list another where both are assumed to be ordered by the ordering in of children list
     * 
     * @return true if something was removed
     */
    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove);

    /**
     * This mapping maps each variable that we will want to push during localisation
     * and which is not possible to push to every child (for conjuction this means all exist
     * vars, for disjunction this is all univ vars and those exist vars which do not fulfill
     * the conditions for pushing into every child). This means that for var v
     * childrenToCombineMapping[v] will contain all children which either contain v or contain
     * some exist var dependent on v (this is only if v is univ var). At the beginning, this 
     * mappind does not contain any variables, they have to be added by using functions
     * addExistVarsToChildrenToCombineMapping and addUnivVarsToChildrenToCombineMapping.
     * Function pushVarsWithCombining is then used to create new children and pushing variables.
     */
    std::unordered_map<Variable, std::list<QuantifierTreeNode*>> childrenToCombineMapping;
    //void initialiseChildrenToCombineMapping(const VariableSet &varsToCombine);
    /**
     * Adds each exist var y from eVarsToAdd to childrenToCombineMapping, i.e.
     * childrenToCombineMapping[y] will contain all children which contain y. Later,
     * this will be used to create a new child with these children to which y will be pushed.
     */
    void addExistVarsToChildrenToCombineMapping(const VariableSet &eVarsToAdd);
    /**
     * Adds each univ var x from uVarsToAdd to childrenToCombineMapping, i.e.
     * childrenToCombineMapping[x] will contain all children which either contain x
     * or some exist var y which depends on x.
     */
    void addUnivVarsToChildrenToCombineMapping(const VariableSet &uVarsToAdd);
    /**
     * Iteratively creates for each var v a new child combining children from 
     * childrenToCombineMapping[v] to which v is pushed starting from the 
     * var v which has the smallest children containing it.
     * It assumes that it starts with only exist vars that we want to push, univ
     * vars that can be pushed are also added and pushed in this function.
     * For conjunction, it will push every pushable univ var (pushable = it does
     * not depend on any exist vars in this node of quantifier tree) to each child,
     * for disjunction it will combine children containing them.
     */
    void pushVarsWithCombining();
    /**
     * @brief For each exist checks if it is pushable according to conditions for disjunction and pushes those that are
     * 
     * It is assumed that childrenToCombineMapping were initialised with addExistVarsToChildrenToCombineMapping(getExistVars())
     * and that it is used only for disjunction.
     */
    void pushExistVarsSeparately(const VariableSet &uVarsOutsideThisSubtree);

    /* Keeps for each (original) child the universal variables that are outside the subtree rooted in the child 
     * (but only those that occur in the subtree). Used for checking conditions of pushing existential
     * variables separately for disjunction and also to call localise for a child. For disjunction, this
     * is computed during localisation for each original child and then it can be reused during the calls 
     * to localisation. For conjunction, we compute this only when we want to localise. Therefore, to
     * access uVarsOutsideChildSubtree[child] use the function getUVarsOutsideChildSubtree(child, ...).
     */
    std::unordered_map<QuantifierTreeNode*, VariableSet> uVarsOutsideChildSubtree;
    // to save on time, only access uVarsOutsideChildSubtree[child] trough this function, it will compute it only once
    VariableSet &getUVarsOutsideChildSubtree(QuantifierTreeNode* child, std::list<QuantifierTreeNode*> originalChildren, const VariableSet &uVarsOutsideThisSubtree);


    std::ostream& print(std::ostream& out) const override;

    /**
     * @brief Construct a new Quantifier Tree inside a tree
     * 
     * We assume that children used to be children of parent, siblings are the remaining children of parent,
     * after creating this tree, it should be added to parent as a child.
     */
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr, 
                        const std::list<QuantifierTreeNode*> &siblings, QuantifierTree *parent);

public:
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr);
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManipulator &qvManipulator);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    ~QuantifierTree();

    void localise(const VariableSet &uVarsOutsideThisSubtree) override;

    /**
     * @brief Changes this instance of QuantifierTree to the instance of Formula,
     * it is assumed that this instance was created by new, after calling this function 
     * the pointer to the QuantifierTree is invalidated (do NOT delete it!)
     * 
     * @param mgr The Cudd manager used for creating matrix of Formula
     * @return Pointer to the resulting instance of Formula, needs to be deleted after it was used
     */
    QuantifierTreeFormula* changeToFormula(Cudd &mgr) override;
    void negate() override;
};


#endif
