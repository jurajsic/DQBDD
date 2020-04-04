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

#include "DQBDDvariable.hpp"
#include "DQBDDformula.hpp"

class QuantifierTreeFormula;

/**
 * @brief Base class for nodes in quantifier trees
 */
class QuantifierTreeNode : virtual public QuantifiedVariablesManipulator {
public:
    QuantifierTreeNode(QuantifiedVariablesManager &qvmgr);
    virtual ~QuantifierTreeNode() = default;
    // pushes quantifiers inside the subtree rooted in this node
    virtual void localise() = 0;
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
};

/**
 * @brief A node that is also a formula (used for terminal nodes)
 */
class QuantifierTreeFormula : public QuantifierTreeNode, public Formula {
public:
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    void localise() override;
    QuantifierTreeFormula* changeToFormula(Cudd &) override;
    void negate() override;

    VariableSet const &getSupportSet() override;
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
     * @brief Removes from one ordered list another where both are assumed to be ordered by the ordering in of children list
     * 
     * @return true if something was removed
     */
    bool removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove);

    std::ostream& print(std::ostream& out) const override;

public:
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr);
    QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManipulator &qvManipulator);
    QuantifierTree(const QuantifierTree&) = delete;
    QuantifierTree& operator=(const QuantifierTree&) = delete;

    ~QuantifierTree();

    void localise() override;

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

    /**
     * @brief Adds another child to the list of children
     * 
     * If child has the same operator as this quantifier tree and does not contain a quantifier prefix,
     * the children of this child are added instead, while child is deleted.
     */
    void addChild(QuantifierTreeNode *child);
};


#endif
