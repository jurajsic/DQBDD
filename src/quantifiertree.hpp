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
#include <memory>

#include "DQBDDvariable.hpp"
#include "DQBDDformula.hpp"

// TODO maybe add needToRecomputeSupportSets for quantifierTree also, so that it is not computed
// for combine trees in disjunction, as it will not needed there (I think) supportSet and uVarSupportSet

/* TODO how to implement shared subtrees
 *     - change children into QuantifierTreeConnection
 *     - probably keep the number of parents???
 *     - the next thing is probably wrong, have to think it trough again:
 *     - add something like turnToNNF(??over which QuantifierTreeConnection??)
 *          - this function is called only on the subtree without quantifiers (only in root)
 *          - keep the number of parents which were already processed (if it turns to 0, we start pushing for the node)
 *          - i have to also keep the set of parents (or connections to parents) which are waiting for negation to be pushed further (i.e. they have negation on the parent-this node connection)
 *          - if "the number of parents which were already processed" == "the number of parents"
 *              - if "the set of parents which are waiting for negation to be pushed" == "the number of parents"
 *                  - we just delete this set and call pushNegationsFurther() for this node, as it needs to be pushed for every parent
 *              - else if "the set of parents which are waiting for negation to be pushed" > 0
 *                  - take the set of (connections to) parents that are waiting
 *                  - create new node with same children connections as this node (except we have to change 'parent' field of course)
 *                  - "the set of parents which are waiting for negation to be pushed" should be zero + the actual set should be empty for the new node
 *                  - update the connection from the set with the new node as a child, change isChildNegated to false (should be true here)
 *                  - delete this set and update number of parents for both nodes
 *                  - call pushNegationsFurther() for the newly created node
 *              - call turnToNNF for each child
 *          - else 
 *              - if isChildNegated==true in the connnection to parent which called turnToNNF(), change it to false and add it to the set of connections of parent that are waiting
 *              - decrement "the number of parents which were already processed"
 *     - pushNegationsFurther()
 *          - takes each children connection and changes the value of isChildNegated to its complement
 *     - update localisation, we need to check for number of parents during pushing, if it is more than one, we need to create new node
 *     - update changeToFormula() .... how?????, I need parents!!
 * 
 *     - maybe pushUVar,pushEVar should also be implemented in QuantifierTreeConnection where it would do the shit about splitting child if it has more than one parent??
 */

class QuantifierTreeNode;
class QuantifierTreeFormula;

class QuantifierTreeConnection {
private:
    // you sure you want shared?? what happens in changeToFormula?? won't there be cyclical ownership? -- how will i 
    //QuantifierTreeNode* parent;
    bool isNegated;
    std::shared_ptr<QuantifierTreeNode> node;

    // TODO name and what exactly needed
    bool hasParent = false;
    std::list<std::list<QuantifierTreeConnection>::iterator>::iterator iterToChildToParent;

    QuantifierTreeConnection() = delete;


    QuantifierTreeConnection(std::shared_ptr<QuantifierTreeNode> node, bool isNegated);
    QuantifierTreeConnection(Variable var, bool isNegated);

public:
    // TODO copy and assignment constructor should not be deleted? because of parent or something??
    //QuantifierTreeConnection(const QuantifierTreeConnection&) = delete;
    //QuantifierTreeConnection& operator=(const QuantifierTreeConnection&) = delete;

    // but move is alright, because it moves??
    //QuantifierTreeConnection& operator=(QuantifierTreeConnection&&) = default;

    static QuantifierTreeConnection andQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr);
    static QuantifierTreeConnection orQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr);
    static QuantifierTreeConnection varQT(bool isNegated, const Variable &var, const Cudd &mgr, QuantifiedVariablesManager &qvMgr);
    static QuantifierTreeConnection negatedQT(const QuantifierTreeConnection &operand);

    //bool isNegated() { return isNegated; }
    //std::shared_ptr<QuantifierTreeNode> getNode() { return node; }

    void changeNodeToFormula();
    void turnToNNF();
    void localise(const VariableSet &uVarsOutsideThisSubtree);

    //TODO add operations here

    // HMMMM, if we only keep this as children connections and only update the child in the connection, there is no problem, right??? no need to delete from children and we do not keep the set of parents
};

/**
 * @brief Base class for nodes in quantifier trees
 */
class QuantifierTreeNode : virtual public QuantifiedVariablesManipulator {
protected:

    // the connections of parents to this node
    std::list<std::weak_ptr<QuantifierTreeConnection>> parentsConnections;

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
    virtual std::shared_ptr<QuantifierTreeFormula> changeToFormula(Cudd &mgr) = 0;
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
    std::shared_ptr<QuantifierTreeFormula> changeToFormula(Cudd &) override;
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
    // the connections to children of the root
    std::list<std::shared_ptr<QuantifierTreeConnection>> childrenConnections;

    // if isConj==true, the root has assigned conjuction, otherwise disjunction
    bool isConj;

    /* if needsToLocalise==true, then in localise(), we need to check and push quantifiers
     * further down, otherwise it is assumed that nothing can be pushed (newly created
     * combined children during localisation have this attribute)
     */
    bool needsToLocalise = true;

    /**
     * @brief Adds another child to the list of children
     * 
     * If collapseChildren==true and the child is without quantifiers and has the same operator as this quantifier tree,
     * then the children of this child are added instead, while child is deleted.
     */
    void addChild(std::shared_ptr<QuantifierTreeConnection> childConnection, bool collapseChildren = true);

    /**
     * @brief Removes from one ordered list of subset of childrent another where both are assumed to be ordered by the ordering of children list
     * 
     * @return true if something was removed
     */
    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<std::weak_ptr<QuantifierTreeConnection>> &listToRemoveFrom, const std::list<std::weak_ptr<QuantifierTreeConnection>> &listOfItemsToRemove);

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
    std::unordered_map<Variable, std::list<std::weak_ptr<QuantifierTreeConnection>>> childrenToCombineMapping;
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
     * (but only those that occur in the subtree). Used for checking conditions of pushing existential variables 
     * separately for disjunction and also to call localise for a child. For disjunction, this is computed during 
     * localisation for each original child and then it can be reused during the calls to localisation. For 
     * conjunction, we compute this only when we want to localise and for newly created children, which allows 
     * to remove non-renamed copies of univ vars which should not be passed during localisation. Therefore, to 
     * access uVarsOutsideChildSubtree[child] use the function getUVarsOutsideChildSubtree(child, ...).
     */
    std::unordered_map<std::weak_ptr<QuantifierTreeConnection>, VariableSet> uVarsOutsideChildSubtree;
    // to save on time, only access uVarsOutsideChildSubtree[child] trough this function, it will compute it only once
    VariableSet &getUVarsOutsideChildSubtree(std::shared_ptr<QuantifierTreeConnection> child, const std::list<std::shared_ptr<QuantifierTreeConnection>> &childAndSiblings, const VariableSet &uVarsOutsideThisSubtree);


    std::ostream& print(std::ostream& out) const override;

    /**
     * @brief Construct a new Quantifier Tree inside a tree
     * 
     * We assume that children used to be children of parent, siblings are the remaining children of parent,
     * after creating this tree, it should be added to parent as a child.
     */
    //QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr, 
    //                    const std::list<QuantifierTreeNode*> &siblings, QuantifierTree *parent);

public:
    /**
     * @brief Construct a new Quantifier Tree object
     * 
     * @param collapseChildren - should we collapse children (i.e. the children of children with the same operation are moved into this node)
     */
    QuantifierTree(bool isConj, const std::list<std::shared_ptr<QuantifierTreeConnection>> &childrenConnections, QuantifiedVariablesManager &qvMgr, bool collapseChildren = true);
    /**
     * @brief Construct a new Quantifier Tree object from some existing QuantifiedVariablesManipulator (i.e. it copies the quantifiers from the manipulator here)
     * 
     * @param qvManipulator - the QuantifiedVariablesManipulator whose quantifiers will be copied here
     * @param collapseChildren - should we collapse children (i.e. the children of children with the same operation are moved into this node)
     */
    QuantifierTree(bool isConj, const std::list<std::shared_ptr<QuantifierTreeConnection>> &childrenConnections, QuantifiedVariablesManipulator &qvManipulator, bool collapseChildren = true);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    //~QuantifierTree();

    void localise(const VariableSet &uVarsOutsideThisSubtree) override;

    /**
     * @brief Changes this instance of QuantifierTree to the instance of Formula,
     * it is assumed that this instance was created by new, after calling this function 
     * the pointer to the QuantifierTree is invalidated (do NOT delete it!)
     * 
     * @param mgr The Cudd manager used for creating matrix of Formula
     * @return Pointer to the resulting instance of Formula, needs to be deleted after it was used
     */
    std::shared_ptr<QuantifierTreeFormula> changeToFormula(Cudd &mgr) override;
    void negate() override;
};


#endif
