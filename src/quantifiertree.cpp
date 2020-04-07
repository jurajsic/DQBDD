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

#include <algorithm>

#include "quantifiertree.hpp"
#include "DQBDDexceptions.hpp"

QuantifierTreeNode::QuantifierTreeNode(QuantifiedVariablesManager &qvmgr) : QuantifiedVariablesManipulator(qvmgr) {}

void QuantifierTreeNode::pushExistVar(Variable var) {
    // only put var in this node if it is used in this node
    if (getSupportSet().contains(var)) {
        addExistVar(var);
    }
}

void QuantifierTreeNode::pushUnivVar(Variable var) {
    /* We do not need to rename universal variables,
     * because we assume all existential variables
     * depending on it in this subtree are not in 
     * different subtrees (which is always true if we
     * created this tree by localising, because only 
     * for disjunction we create two subtrees with
     * same existential variables, but for universal
     * variable we cannot push two copies inside the subtrees
     * for disjunction.
     */ 
    if (getSupportSet().contains(var)) { // only put var in this node if it is used in this node
        addUnivVar(var);
    } else { // otherwise we can basically eliminate it by removing dependencies 
        for (const Variable &dependentExistVar : getUnivVarDependencies(var)) {
            if (getSupportSet().contains(dependentExistVar)) { // we assume this is the only tree that contains dependentExistVar
                removeDependency(dependentExistVar,var);
            }
        }
    }
}

/*******************************************/
/*******************************************/
/************ QUANTIFIER TREE **************/
/*******************************************/
/*******************************************/ 

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr) : QuantifiedVariablesManipulator(qvMgr), QuantifierTreeNode(qvMgr), isConj(isConj) {
    supportSet = {};
    if (children.size() < 2) {
        throw DQBDDexception("You cannot create a quantifier tree with only one or zero operands");
    }
    for (QuantifierTreeNode *child : children) {
        addChild(child);
    }
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManipulator &qvManipulator) : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), isConj(isConj) {
    supportSet = {};
    if (children.size() < 2) {
        throw DQBDDexception("You cannot create a quantifier tree with only one or zero operands");
    }
    for (QuantifierTreeNode *child : children) {
        addChild(child);
    }
}

QuantifierTree::~QuantifierTree() {
    for (QuantifierTreeNode *child : children) {
        delete child;
    }
}

bool QuantifierTree::removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove) {
    auto childrenCopy = children;
    auto orderIter = childrenCopy.begin();
    auto listToRemoveFromIter = listToRemoveFrom.begin();
    auto listOfItemsToRemoveIter = listOfItemsToRemove.begin();

    bool somethingWasDeleted = false;

    while (orderIter != childrenCopy.end() || listToRemoveFromIter != listToRemoveFrom.end() || listOfItemsToRemoveIter != listOfItemsToRemove.end()) {
        if (*orderIter == *listToRemoveFromIter && *listToRemoveFromIter == *listOfItemsToRemoveIter) {
            auto iterToDelete = listToRemoveFromIter;
            ++listToRemoveFromIter;
            ++listOfItemsToRemoveIter;
            listToRemoveFrom.erase(iterToDelete);
            somethingWasDeleted = true;
        } else if (*orderIter == *listToRemoveFromIter) {
            ++listToRemoveFromIter;
        } else if (*orderIter == *listOfItemsToRemoveIter) {
            ++listOfItemsToRemoveIter;
        }

        ++orderIter;
    }
    return somethingWasDeleted;
}

void QuantifierTree::localise() {
    removeUnusedVars();
    if (isConj) {
        /****************************************************************************
         * For conjunction, first find for each exist var the set of children that  
         * contain it and merge them into new QuantifierTree with same operator and 
         * then push universal variables on which none of the leftover existential 
         * variables depend into children.
         ****************************************************************************/

        // for each existential variable, create a set of children which contain it
        std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingExistVar;
        for (QuantifierTreeNode *child : children) {
            VariableSet childSupportSet = child->getSupportSet();
            for (const Variable &eVar : getExistVars()) {
                if (childSupportSet.contains(eVar)) {
                    childrenContainingExistVar[eVar].push_back(child);
                }
            }
        }
        
        // push exist vars while it is possible
        while (!getExistVars().empty()) {
            // exist variable that has minimal number of occurences in children
            auto eVarToPushIter = std::min_element(getExistVars().begin(), getExistVars().end(), 
                        [&](Variable v1, Variable v2) { 
                            return (childrenContainingExistVar[v1].size() < childrenContainingExistVar[v2].size());
                        }
                     );

            std::list<QuantifierTreeNode*> &childrenToCombine = childrenContainingExistVar[*eVarToPushIter];

            // if we hit the exist variable that is in all children (meaning 
            // that all other exist variables are in all children too) we stop
            if (childrenToCombine.size() == children.size()) {
                break;
            }

            // if *eVarToPushIter is needed to be pushed to only one child, there is no need to create new children
            if (childrenToCombine.size() == 1) {
                (*childrenToCombine.begin())->pushExistVar(*eVarToPushIter);
                removeExistVar(*eVarToPushIter);
                continue;
            }
            
            // create a new child to which this existential variable is pushed
            QuantifierTree *newChild = new QuantifierTree(isConj, childrenToCombine, *qvMgr);
            newChild->pushExistVar(*eVarToPushIter);
            removeExistVar(*eVarToPushIter);

            // erase the children that were connected to create a new child
            // and add the new child if some children were connected
            for (const Variable &eVarToBePushed : getExistVars()) {
                if (removeFromOrderedListOtherOrderedListUsingChildrenOrder(childrenContainingExistVar[eVarToBePushed], childrenToCombine)) {
                    childrenContainingExistVar[eVarToBePushed].push_back(newChild);
                }
            }
            removeFromOrderedListOtherOrderedListUsingChildrenOrder(children, childrenToCombine);
            children.push_back(newChild);

            childrenContainingExistVar.erase(*eVarToPushIter);
        }
        
        // get all universal variables on which some leftover exist var depends (these vars cant be pushed)
        VariableSet dependentUnivVars; 
        for (const Variable &existVar : getExistVars()) {
            dependentUnivVars.insert(getExistVarDependencies(existVar).begin(), getExistVarDependencies(existVar).end());
        }

        // push possible univ vars (they do not depend on any exist var in current node)
        VariableSet univVars = getUnivVars();
        for (const Variable &univVarToPush : univVars) {
            if (!dependentUnivVars.contains(univVarToPush)) { // no leftover existential variable depends on univVarToPush --> can be pushed
                // push it into every child, pushUnivVar() will take care of deciding whether to keep it there or not
                for (QuantifierTreeNode *child : children) {
                    child->pushUnivVar(univVarToPush);
                }
                removeUnivVar(univVarToPush);
            }
        }
    } else {
        /***********************************************************************************
         * For disjunction, we can easily push every existential quantifier into children
         * and then we create for each univ var new subtrees from children which contain them
         * in which we push them.
         ***********************************************************************************/

        // push all existential variables into children
        VariableSet existVars = getExistVars();
        for (const Variable &existVarToPush : existVars) {
            for (QuantifierTreeNode *child : children) {
                /* There is no need to rename them, because universal vars on which this exist var 
                 * depends will either stay as ancestor/same level or will be pushed into some sibling 
                 * and then the dependency will be removed. Also, pushExistVar() takes care whether
                 * to actually push it inside or not.       
                 */
                child->pushExistVar(existVarToPush);
            }
            removeExistVar(existVarToPush);
        }

        // for each universal variable find the set of children that contain it or contain exist variable which depends
        // on it and merge them into new QuantifierTree with same operator and remove it if it was pushed
        std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingUnivVar;
        for (QuantifierTreeNode *child : children) {
            VariableSet childSupportSet = child->getSupportSet();
            for (const Variable &uVar : getUnivVars()) {
                // if the child contains uVar...
                if (childSupportSet.contains(uVar)) {
                    childrenContainingUnivVar[uVar].push_back(child);
                    continue;
                }
                // ... or some exist var that depends on uVar
                for (const Variable &eVar : getUnivVarDependencies(uVar)) {
                    if (childSupportSet.contains(eVar)) {
                        childrenContainingUnivVar[uVar].push_back(child);
                        continue;
                    }
                }
            }
        }
        
        // push all universal variables (except those on which some exist var from every child depends)
        while (!getUnivVars().empty()) {
            // univ variable that has minimal number of occurences in children
            auto uVarToPushIter = std::min_element(getUnivVars().begin(), getUnivVars().end(), 
                        [&](Variable v1, Variable v2) { 
                            return (childrenContainingUnivVar[v1].size() < childrenContainingUnivVar[v2].size());
                        }
                     );

            std::list<QuantifierTreeNode*> &childrenToCombine = childrenContainingUnivVar[*uVarToPushIter];

            // if we hit the univ variable that is in all children (meaning 
            // that all other univ variables are in all children too) we stop
            if (childrenToCombine.size() == children.size()) {
                break;
            }

            // if *uVarToPushIter is needed to be pushed to only one child, there is no need to create new child
            if (childrenToCombine.size() == 1) {
                (*childrenToCombine.begin())->pushUnivVar(*uVarToPushIter);
                removeUnivVar(*uVarToPushIter);
                continue;
            }

            // create a new child to which this universal variable is pushed
            QuantifierTree *newChild = new QuantifierTree(isConj, childrenToCombine, *qvMgr);
            newChild->pushUnivVar(*uVarToPushIter);
            removeUnivVar(*uVarToPushIter);

            // erase the children that were connected to create a new child
            // and add the new child if some children were deleted
            for (const Variable &uVarToBePushed : getUnivVars()) {
                if (removeFromOrderedListOtherOrderedListUsingChildrenOrder(childrenContainingUnivVar[uVarToBePushed], childrenToCombine)) {
                    childrenContainingUnivVar[uVarToBePushed].push_back(newChild);
                }
            }
            removeFromOrderedListOtherOrderedListUsingChildrenOrder(children, childrenToCombine);
            children.push_back(newChild);

            childrenContainingUnivVar.erase(*uVarToPushIter);
        }
    }

    // recursively call localise for each child
    for (QuantifierTreeNode *child : children) {
        child->localise();
    }
}

QuantifierTreeFormula* QuantifierTree::changeToFormula(Cudd &mgr) {
    BDD matrix;
    if (isConj) {
        matrix = mgr.bddOne();
    } else {
        matrix = mgr.bddZero();
    }

    auto childIter = children.begin();
    while (childIter != children.end()) {
        // get the formula from child
        QuantifierTreeFormula *childFormula = (*childIter)->changeToFormula(mgr);

        // remove the child from the list of children
        auto childToRemoveIter = childIter;
        ++childIter;
        children.erase(childToRemoveIter);

        switch (qvMgr->options.treeElimChoice)
        {
        case TreeElimChoice::None:
            break;
        case TreeElimChoice::Simple: {
            // eliminate possible existential variables
            // TODO: maybe also universal without dependencies
            VariableSet existVarsToEliminate = childFormula->getPossibleExistVarsToEliminate();
            while (existVarsToEliminate.size() != 0) {
                childFormula->eliminateExistVars(existVarsToEliminate);
                childFormula->removeUnusedVars();
                existVarsToEliminate = childFormula->getPossibleExistVarsToEliminate();
            }
            break;
        }
        case TreeElimChoice::All: {
            // eliminate all possible universal and existential variables
            childFormula->eliminatePossibleVars();
            break;
        }
        default: {
            throw DQBDDexception("Selected choice of variables to eliminate in tree is not supported.");
            break;
        }
        }
        
        // check if we can stop (the resulting formula after adding the child formula is simple)
        if (isConj) {
            matrix &= childFormula->getMatrix();
            if (matrix.IsZero()) {
                delete childFormula;
                break;
            }
        } else {
            matrix |= childFormula->getMatrix();
            if (matrix.IsOne()) {
                delete childFormula;
                break;
            }
        }

        // move the quantifier prefix up
        for (const Variable &eVar : childFormula->getExistVars()) {
            addExistVar(eVar);
        }
        for (const Variable &uVar : childFormula->getUnivVars()) {
            addUnivVar(uVar);
        }

        delete childFormula;
    }
    
    // delete leftover children
    while (childIter != children.end()) {
        auto childToRemoveIter = childIter;
        ++childIter;
        delete (*childToRemoveIter);
        children.erase(childToRemoveIter);
    }

    QuantifierTreeFormula* f = new QuantifierTreeFormula(mgr, *qvMgr);
    f->setMatrix(matrix);
    for (const Variable &uVar : getUnivVars()) {
        f->addUnivVar(uVar);
        //removeUnivVar(uVar);
    }
    for (const Variable &eVar : getExistVars()) {
        f->addExistVar(eVar);
        //removeExistVar(eVar);
    }
    delete this; // removing is done here
    f->removeUnusedVars();
    return f;
}

void QuantifierTree::negate() {
    isConj = !isConj;
    for (QuantifierTreeNode *child : children) {
        child->negate();
    }
}

void QuantifierTree::addChild(QuantifierTreeNode *child) {
    // add variables of the child to the support set here
    for (const Variable var : child->getSupportSet()) {
        supportSet.insert(var);
    }

    // check if this child is not quantifier tree with the same operator like here and empty quantifier prefix
    auto treeChild = dynamic_cast<QuantifierTree*>(child);
    if (treeChild != nullptr && treeChild->isConj == isConj 
                && treeChild->getExistVars().empty() && treeChild->getUnivVars().empty()) {
        // if it is, set its children as current children, not itself
        for (auto childOfChild : treeChild->children) {
            children.push_back(childOfChild);
        }
        treeChild->children.clear();

        // and finally delete it
        delete child;
    } else {
        // otherwise just add this child to children
        children.push_back(child);
    }
}

std::ostream& QuantifierTree::print(std::ostream& out) const {
    std::string op;
    if (isConj) {
        op = '&';
    } else {
        op = '|';
    }

    auto size = children.size();
    for (QuantifierTreeNode *child : children) {
        out << std::string("(") << *child << std::string(")");
        // print the last element without operator
        if (size != 1) {
            out << std::string(" ") << op << std::string(" ");
        }
        --size;
    }
    return out;
}

/*********************************************************/
/*********************************************************/
/************     QUANTIFIERTREEFORMULA     **************/
/*********************************************************/
/*********************************************************/

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr)
                            : QuantifiedVariablesManipulator(qvmgr), QuantifierTreeNode(qvmgr), Formula(mgr, qvmgr) {}

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator)
                            : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), Formula(mgr, qvManipulator) {}

void QuantifierTreeFormula::localise() {
    removeUnusedVars();
}

QuantifierTreeFormula* QuantifierTreeFormula::changeToFormula(Cudd &) {
    return this;
}

void QuantifierTreeFormula::negate() {
    setMatrix(!getMatrix());
}

VariableSet const& QuantifierTreeFormula::getSupportSet() {
    return Formula::getSupportSet();
}