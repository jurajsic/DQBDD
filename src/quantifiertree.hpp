/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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

#include "dqbddvariable.hpp"
#include "dqbddformula.hpp"

namespace dqbdd {

class QuantifierTreeFormula;

/**
 * @brief Base class for nodes in quantifier trees
 */
class QuantifierTreeNode : virtual public QuantifiedVariablesManipulator {
protected:
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
    void pushExistVar(const Variable &var);
    // pushes the universal variable var into this node
    void pushUnivVar(const Variable &var);
    /**
     * @brief Changes this instance of node into a formula
     */
    virtual QuantifierTreeFormula* changeToFormula(Cudd &mgr) = 0;

    virtual VariableSet const &getUVarsSupportSet();

    virtual QuantifierTreeNode* getCopy() = 0;
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
    void localise(const VariableSet&) override;
    QuantifierTreeFormula* changeToFormula(Cudd &) override;
    QuantifierTreeNode* getCopy() override;

    VariableSet const &getSupportSet() override;
    VariableSet const &getUVarsSupportSet() override;
};

/**
 * @brief A connection for quantifier tree to its child
 */
class QuantifierTreeConnection {
private:
    // number of quantifier trees that have this connection as child
    unsigned numOfParents = 0;
public:
    QuantifierTreeNode *child;

    QuantifierTreeConnection(QuantifierTreeNode *child);

    /**
     * @brief Changes child to quantifier tree formula
     * 
     * @return the new child
     */
    QuantifierTreeFormula* changeChildToQuantifierTreeFormula(Cudd &mgr);


    void decreaseNumOfParents();
    void increaseNumOfParents();
    unsigned getNumOfParents();
};

/**
 * @brief A node that represents a root of a quantifier subtree
 */
class QuantifierTree : public QuantifierTreeNode {
private:
    // the connections with children of this tree
    std::list<QuantifierTreeConnection*> childrenConnections;

    // if isConj==true, this tree has assigned conjuction, otherwise disjunction
    bool isConj;

    /* if needsToLocalise==true, then in localise(), we need to check and push quantifiers
     * further down, otherwise it is assumed that nothing can be pushed (newly created
     * combined children during localisation have this attribute)
     */
    bool needsToLocalise = true;

    /**
     * @brief Adds another child to the list of children (assuming *this is a root)
     * 
     * @param child - the connection with the root of subtree of the child to add
     * @param collapseChildren - if true and the child has the same operator as this quantifier 
     * tree and does not contain a quantifier prefix, then the children of the child are added instead
     * @param changeNumOfParents if true, the newly added child will update its numOfParents (increments it by 1)
     */
    void addChild(QuantifierTreeConnection *childConnection, bool collapseChildren = true, bool changeNumOfParents = true);

    // decrements the number of parents for a child and possibly deletes the child and its childConnection (but only if this was the last parent)
    void deleteChildConnection(QuantifierTreeConnection *childConnection);

    /**
     * @brief Removes from one ordered list another where both are assumed to be ordered by the ordering of childrenConnections list
     * 
     * @return true if something was removed
     */
    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeConnection*> &listToRemoveFrom, std::list<QuantifierTreeConnection*> &listOfItemsToRemove);

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
    std::unordered_map<Variable, std::list<QuantifierTreeConnection*>> childrenToCombineMapping;
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
     * (but only those that occur in the subtree). Used for checking conditions of pushing existential variables 
     * separately for disjunction and also to call localise for a child. For disjunction, this is computed during 
     * localisation for each original child and then it can be reused during the calls to localisation. For 
     * conjunction, we compute this only when we want to localise and for newly created children, which allows 
     * to remove non-renamed copies of univ vars which should not be passed during localisation. Therefore, to 
     * access uVarsOutsideChildSubtree[child] use the function getUVarsOutsideChildSubtree(child, ...).
     */
    std::unordered_map<QuantifierTreeConnection*, VariableSet> uVarsOutsideChildSubtree;
    // to save on time, only access uVarsOutsideChildSubtree[child] trough this function, it will compute it only once
    VariableSet &getUVarsOutsideChildSubtree(QuantifierTreeConnection* childConnection, const std::list<QuantifierTreeConnection*> &childAndSiblingsConnections, const VariableSet &uVarsOutsideThisSubtree);


    std::ostream& print(std::ostream& out) const override;

public:
    /**
     * @brief Construct a new Quantifier Tree object
     * 
     * @param collapseChildren - should we collapse children (i.e. the children of children with the same operation are moved into this node)
     * @param changeNumOfParents - should we update the numOfParents of the newly added children (i.e. increment by 1)
     */
    QuantifierTree(bool isConj, std::list<QuantifierTreeConnection*> childrenConnections, QuantifiedVariablesManager &qvMgr, bool collapseChildren = true, bool changeNumOfParents = true);
    /**
     * @brief Construct a new Quantifier Tree object from some existing QuantifiedVariablesManipulator (i.e. it copies the quantifiers from the manipulator here)
     * 
     * @param qvManipulator - the QuantifiedVariablesManipulator whose quantifiers will be copied here
     * @param collapseChildren - should we collapse children (i.e. the children of children with the same operation are moved into this node)
     */
    QuantifierTree(bool isConj, std::list<QuantifierTreeConnection*> childrenConnections, QuantifiedVariablesManipulator &qvManipulator, bool collapseChildren = true);
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
    
    QuantifierTreeNode* getCopy() override;
};

} // namespace dqbdd

#endif
