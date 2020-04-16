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

//#include<iostream>

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
        VariableSet dependentExistVars = getUnivVarDependencies(var);
        for (const Variable &dependentExistVar : dependentExistVars) {
            //std::cout << dependentExistVar << " is in support set:" << getSupportSet().contains(dependentExistVar) << std::endl;
            if (getSupportSet().contains(dependentExistVar)) { // we assume this is the only tree that contains dependentExistVar
                removeDependency(dependentExistVar,var);
            }
        }
    }
}

VariableSet const& QuantifierTreeNode::getUVarsOutsideThisSubtree() const {
    return uVarsOutsideThisSubtree;
}

VariableSet const& QuantifierTreeNode::getUVarsSupportSet() {
    return uVarsSupportSet;
}

void QuantifierTreeNode::addToUVarsOutsideThisSubtree(const Variable &varToAdd) {
    addToUVarsOutsideThisSubtree(VariableSet{varToAdd});
}

/*******************************************/
/*******************************************/
/************ QUANTIFIER TREE **************/
/*******************************************/
/*******************************************/ 

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr) : QuantifiedVariablesManipulator(qvMgr), QuantifierTreeNode(qvMgr), isConj(isConj) {
    supportSet = {};
    uVarsSupportSet = {};
    if (children.size() < 2) {
        throw DQBDDexception((std::string("You cannot create a quantifier tree with ") + std::to_string(children.size()) + std::string(" operands")).c_str());
    }
    for (QuantifierTreeNode *child : children) {
        addChild(child);
    }
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManipulator &qvManipulator) : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), isConj(isConj) {
    supportSet = {};
    uVarsSupportSet = {};
    if (children.size() < 2) {
        throw DQBDDexception((std::string("You cannot create a quantifier tree with ") + std::to_string(children.size()) + std::string(" operands")).c_str());
    }
    for (QuantifierTreeNode *child : children) {
        addChild(child);
    }
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr, 
                        const std::list<QuantifierTreeNode*> &siblings, QuantifierTree *parent) :QuantifiedVariablesManipulator(qvMgr), QuantifierTreeNode(qvMgr), isConj(isConj) {
    supportSet = {};
    uVarsSupportSet = {};
    if (children.size() < 2) {
        throw DQBDDexception((std::string("You cannot create a quantifier tree with ") + std::to_string(children.size()) + std::string(" operands")).c_str());
    }
    
    // we do not need to update uVarsOutsideThisSubtree of children because of 
    // the assumption where this constructor can be used 
    // so we only set it for this node
    uVarsOutsideThisSubtree = parent->getUVarsOutsideThisSubtree();
    for (QuantifierTreeNode *sibling : siblings) {
        uVarsOutsideThisSubtree.insert(sibling->getUVarsSupportSet().begin(),
                                       sibling->getUVarsSupportSet().end());
    }

    for (QuantifierTreeNode *child : children) {
        supportSet.insert(child->getSupportSet().begin(),
                          child->getSupportSet().end());
        uVarsSupportSet.insert(child->getUVarsSupportSet().begin(),
                               child->getUVarsSupportSet().end());
        this->children.push_back(child);
    }
}

QuantifierTree::~QuantifierTree() {
    for (QuantifierTreeNode *child : children) {
        delete child;
    }
}

bool QuantifierTree::removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove) {
    //std::cout << "Deleting from list of size " << listToRemoveFrom.size() << " the list of size " << listOfItemsToRemove.size() << " using children order of size " << children.size() << std::endl;
    
    auto childrenCopy = children;
    auto orderIter = childrenCopy.begin();
    auto listToRemoveFromIter = listToRemoveFrom.begin();
    auto listOfItemsToRemoveIter = listOfItemsToRemove.begin();

    bool somethingWasDeleted = false;

    while (orderIter != childrenCopy.end() && listToRemoveFromIter != listToRemoveFrom.end() && listOfItemsToRemoveIter != listOfItemsToRemove.end()) {
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
        //if (listToRemoveFromIter == listToRemoveFrom.end()) std::cout << "we found end of listtoremovefrom" << std::endl;
        //if (listOfItemsToRemoveIter == listOfItemsToRemove.end()) std::cout << "we found end of listofitemtoremove" << std::endl;
        //if (orderIter == childrenCopy.end()) std::cout << "we found end of children" << std::endl;
    }
    return somethingWasDeleted;
}

void QuantifierTree::addPossibleUnivVarsToMapping(std::unordered_map<Variable, std::list<QuantifierTreeNode*>> &childrenToCombineMapping) {
    //std::cout << "Finding new univ variables...." << std::endl;
    
    // get all universal variables on which some leftover exist var depends (these vars cant be pushed)...
    VariableSet dependentUnivVars; 
    for (const Variable &existVar : getExistVars()) {
        dependentUnivVars.insert(getExistVarDependencies(existVar).begin(), getExistVarDependencies(existVar).end());
    }
    // ...add to it variables which are already in mapping...
    for (const auto &univVarInMappingWithMappedContent : childrenToCombineMapping) {
        dependentUnivVars.insert(univVarInMappingWithMappedContent.first);
    }
    // ...and take all universal variables which are not in it
    VariableSet notDependentUnivVars = getUnivVars().minus(dependentUnivVars);

    // for each universal variable that is not in dependentUnivVars get the list of children containing it or its dependency
    if (!notDependentUnivVars.empty()) {
        for (QuantifierTreeNode *child : children) {
            const VariableSet &childSupportSet = child->getSupportSet();
            for (const Variable &uVar : notDependentUnivVars) {
                // if the child contains uVar...
                if (childSupportSet.contains(uVar)) {
                    childrenToCombineMapping[uVar].push_back(child);
                    continue;
                }
                // ... or some exist var that depends on uVar
                for (const Variable &eVar : getUnivVarDependencies(uVar)) {
                    if (childSupportSet.contains(eVar)) {
                        childrenToCombineMapping[uVar].push_back(child);
                        continue;
                    }
                }
            }
        }
    }

    //std::cout << "....found" << std::endl;
}


void QuantifierTree::pushVarsWithCombining(std::unordered_map<Variable, std::list<QuantifierTreeNode*>> &childrenToCombineMapping, bool findNewUnivVars) {
    //std::cout << "We are planning to remove these with combining" << std::endl;
    //for (auto i : childrenToCombineMapping) {
    //    std::cout << "  variable " << i.first << " with " << i.second.size() << " children to combine." << std::endl;
    //}
    
    while (!childrenToCombineMapping.empty()) {
        //std::cout << "These are leftover" << std::endl;
        //for (auto i : childrenToCombineMapping) {
        //    std::cout << "  variable " << i.first << " with " << i.second.size() << " children to combine." << std::endl;
        //}

        // get variable that has minimal number of occurences in children
        auto varToPushAndItsChildren = std::min_element(childrenToCombineMapping.begin(), childrenToCombineMapping.end(), 
                    [](const std::pair<Variable, std::list<QuantifierTreeNode*>> &a, 
                        const std::pair<Variable, std::list<QuantifierTreeNode*>> &b) { 
                        return (a.second.size() < b.second.size());
                    }
                    );
        Variable varToPush = varToPushAndItsChildren->first;
        //std::cout << "Pushing the " << (isVarUniv(varToPush) ? "universal" : "existential") << " variable " << varToPush << " to a new combined child." << std::endl;
        std::list<QuantifierTreeNode*> childrenToCombine = varToPushAndItsChildren->second;
        //std::cout << "The size of children to combine right after getting it: " << childrenToCombine.size() << std::endl;

        // if we hit the variable that is in all children (meaning 
        // that all other variables are in all children too) we stop
        if (childrenToCombine.size() == children.size()) {
            break;
        }

        // if varToPush is needed to be pushed to only one child, there is no need to create new child
        if (childrenToCombine.size() == 1) {
            if (isVarUniv(varToPush)) {
                (*childrenToCombine.begin())->pushUnivVar(varToPush);
                removeUnivVar(varToPush);
            } else if (isVarExist(varToPush)) {
                (*childrenToCombine.begin())->pushExistVar(varToPush);
                removeExistVar(varToPush);
                if (findNewUnivVars) {
                    addPossibleUnivVarsToMapping(childrenToCombineMapping);
                }
            } else {
                throw DQBDDexception("Found variable which is neither universal or existential, this should not happen");
            }
            childrenToCombineMapping.erase(varToPush);
            continue;
        }

        //std::cout << "Finding vars with changed children...." << std::endl;

        // for each variable v that we are planning to push later, we need to update
        // childrenToCombineMapping[v] by removing from them the set of children
        // which we are going to combine and also save them into varsWithChangedChildren
        // if something was deleted to add the new child with combined children to it
        VariableSet varsWithChangedChildren = {};
        for (auto &varToBePushedWithChildren : childrenToCombineMapping) {
            Variable varToBePushed = varToBePushedWithChildren.first;
            //std::cout << "   maybe var " << varToBePushed << " has changed children" << std::endl;
            if (removeFromOrderedListOtherOrderedListUsingChildrenOrder(childrenToCombineMapping[varToBePushed], childrenToCombine)) {
                //std::cout << "   yes it does" << std::endl;
                varsWithChangedChildren.insert(varToBePushed);
            }
        }
        
        //std::cout << ".....found" << std::endl;

        // create a new child to which varToPush is pushed
        removeFromOrderedListOtherOrderedListUsingChildrenOrder(children, childrenToCombine);
        //std::cout << "The size of children to combine right before creating new child: " << childrenToCombine.size() << std::endl;
        QuantifierTree *newChild = new QuantifierTree(isConj, childrenToCombine, *qvMgr, children, this);
        if (isVarExist(varToPush)) {
            newChild->pushExistVar(varToPush);
            removeExistVar(varToPush);
        } else if (isVarUniv(varToPush)) {
            newChild->pushUnivVar(varToPush);
            removeUnivVar(varToPush);
        } else {
            throw DQBDDexception("Found variable which is neither universal or existential, this should not happen");
        }
        children.push_back(newChild);

        
        //std::cout << "...combined tree created" << std::endl;

        // add the new child to each var in varsWithChangedChildren
        for (auto varWithChangedChildren : varsWithChangedChildren) {
            childrenToCombineMapping[varWithChangedChildren].push_back(newChild);
        }

        childrenToCombineMapping.erase(varToPush);

        // at the end, if the pushed variable was existential, some universal variables can become pushable
        if (isVarExist(varToPush) && findNewUnivVars) {
            addPossibleUnivVarsToMapping(childrenToCombineMapping);
        }
    }
}


void QuantifierTree::localiseOR() {
    /************************************************************************************
     * For disjunction, we can push only some existential quantifiers into children, we
     * have to check if the children fulfill the conditions for it to be possible. We 
     * can however still create the new subtrees which contain them to which we can them,
     * but also we do it at the same time as we push univ vars, for those we also create 
     * new subtrees which contain them or some existential variable depending on it.
     ************************************************************************************/

    // for each existential variable map to the list of children containing it
    std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingVarOrDependency;
    for (QuantifierTreeNode *child : children) {
        const VariableSet &childSupportSet = child->getSupportSet();
        for (const Variable &eVar : getExistVars()) {
            if (childSupportSet.contains(eVar)) {
                childrenContainingVarOrDependency[eVar].push_back(child);
            }
        }
    }

    /* first, try pushing existential variables into children if they all fulfill the needed conditions,
        * where for exist. var y with dependency set D_y:
        *     (1) check if all children containing y have uVarsSupportSet/D_y pairwise disjoint
        *     (2) at most one child containing y has some variable from uVarsSupportSet/D_y in uVarsSupportSet
        *       of child that does not contain y or in uVarsOutsideThisSubtree
        */
    VariableSet existVars = getExistVars();
    for (const Variable &existVarToPush : existVars) {
        // if true, existVarToPush fulfills the conditions
        bool canBePushed = true;
        // true if one of the processed children was already problematic
        bool wasThereOneProblem = false;

        // first check the condition
        for (QuantifierTreeNode *child : childrenContainingVarOrDependency[existVarToPush]) {
            /* uVarsToCheck = uVarsSupportSet/D_y
                * - the set of universal variables in the child (or in a dependency set of a exist var in child)
                *   which are not in the dependency set of existVarToPush
                */
            VariableSet uVarsToCheck = child->getUVarsSupportSet().minus(getExistVarDependencies(existVarToPush));

            /* check if uVarsToCheck are outside this child subtree
                *     - if it violates condition (1), then it is true and it either sets wasThereOneProblem
                *       to true and some next child with which it intersects will also evaluate this to true
                *       but wasThereOneProblem will be already set to true OR this is the second child and 
                *       wasThereOneProblem is already true
                *     - if it violates condition (2), then it is true and we set wasThereOneProblem to true
                *       so we know this one child is the one possible child that can violate (2)
                *     - if it does not violate either, then this is false 
                */
            if ( !(child->getUVarsOutsideThisSubtree().intersect(uVarsToCheck).empty()) ) {
                if (wasThereOneProblem) {
                    canBePushed = false;
                    break;
                } else {
                    wasThereOneProblem = true;
                }
            }
        }

        // if it can be pushed, push it
        if (canBePushed) {
            //std::cout << "Pushing " << existVarToPush << " using the weird rule." << std::endl;
            for (QuantifierTreeNode *child : childrenContainingVarOrDependency[existVarToPush]) {
                /* There is no need to rename them, because universal vars on which this exist var 
                * depends will either stay as ancestor/same level or will be pushed into some sibling 
                * and then the dependency will be removed. Also, pushExistVar() takes care whether
                * to actually push it inside or not.       
                */
                child->pushExistVar(existVarToPush);
            }
            removeExistVar(existVarToPush);
            childrenContainingVarOrDependency.erase(existVarToPush);
        }
    }

    addPossibleUnivVarsToMapping(childrenContainingVarOrDependency);

    // push all possible leftover existential and universal vars into a new child created from children containing them
    pushVarsWithCombining(childrenContainingVarOrDependency, true);
}

void QuantifierTree::localiseAND() {
    /****************************************************************************
     * For conjunction, first find for each exist var the set of children that  
     * contain it and merge them into new QuantifierTree with same operator and 
     * then push universal variables on which none of the leftover existential 
     * variables depend into children.
     ****************************************************************************/

    // for each existential variable, create a set of children which contain it
    std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingExistVar;
    for (QuantifierTreeNode *child : children) {
        const VariableSet &childSupportSet = child->getSupportSet();
        for (const Variable &eVar : getExistVars()) {
            if (childSupportSet.contains(eVar)) {
                childrenContainingExistVar[eVar].push_back(child);
            }
        }
    }
    
    // push exist vars while it is possible
    pushVarsWithCombining(childrenContainingExistVar, false);
    
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
                //std::cout << "Pushing " << univVarToPush << " to child " << *child << std::endl;
                //std::cout << child->getSupportSet() << std::endl;
                child->pushUnivVar(univVarToPush);
            }
            removeUnivVar(univVarToPush);
        }
    }
}

void QuantifierTree::localise() {
    removeUnusedVars();
    //std::cout << "Localising at tree: " << *this << " with: " << std::endl
    //          << "   support set:   " << getSupportSet() << std::endl
    //          << "   u support set: " << getUVarsSupportSet() << std::endl
    //          << "   u outside:     " << getUVarsOutsideThisSubtree() << std::endl;
    if (isConj) {
        localiseAND();
    } else {
        localiseOR();
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
            bool somethingWasEliminated = false;
            do {
                // eliminate possible existential variables...
                somethingWasEliminated = childFormula->eliminatePossibleExistVars();
                // ...and universal without dependencies
                somethingWasEliminated |= childFormula->eliminateSimpleUnivVars();
            } while (somethingWasEliminated);
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


void QuantifierTree::addToUVarsOutsideThisSubtree(const VariableSet &varsToAdd) {
    for (auto child : children) {
        child->addToUVarsOutsideThisSubtree(varsToAdd);
    }

    uVarsOutsideThisSubtree.insert(varsToAdd.begin(), varsToAdd.end());
}

void QuantifierTree::addChild(QuantifierTreeNode *child) {
    /* in oldUVarsSupportSet are all universal variables in all the children
     * of this tree before adding new child, therefore all those will be
     * added to uVarsOutsideThisSubtree of the new child
     */
    VariableSet oldUVarsSupportSet = getUVarsSupportSet();

    // add variables of the child to the support sets here
    supportSet.insert(child->getSupportSet().begin(),
                    child->getSupportSet().end());
    uVarsSupportSet.insert(child->getUVarsSupportSet().begin(),
                        child->getUVarsSupportSet().end());

    // update others children uVarsOutsideThisSubtree with the universal
    // variables that are in support set of what will be newly added child
    for (auto oldChild : children) {
        oldChild->addToUVarsOutsideThisSubtree(child->getUVarsSupportSet());
    }

    // check if this child is not quantifier tree with the same operator like here and empty quantifier prefix
    auto treeChild = dynamic_cast<QuantifierTree*>(child);
    if (treeChild != nullptr && treeChild->isConj == isConj 
                && treeChild->getExistVars().empty() && treeChild->getUnivVars().empty()) {

        // if it is, set its children as current children, not itself
        for (auto childOfChild : treeChild->children) {
            children.push_back(childOfChild);
            // we only nned to add to uVarsOutsideThisSubtree the uVarsInSupportSet
            // of old children, because these children of the child have already
            // the uVarsInSupportSet of other children of the child in their uVarsOutsideThisSubtree
            childOfChild->addToUVarsOutsideThisSubtree(oldUVarsSupportSet);
        }

        // and finally delete it
        treeChild->children.clear();
        delete child;
    } else {
        // otherwise just add this child to children
        children.push_back(child);
        child->addToUVarsOutsideThisSubtree(oldUVarsSupportSet);
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
    out << std::string("(");
    for (QuantifierTreeNode *child : children) {
        out << *child;
        // print the last element without operator
        if (size != 1) {
            out << std::string(" ") << op << std::string(" ");
        }
        --size;
    }
    out << std::string(")");
    return out;
}

/*********************************************************/
/*********************************************************/
/************     QUANTIFIERTREEFORMULA     **************/
/*********************************************************/
/*********************************************************/

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr)
                            : QuantifiedVariablesManipulator(qvmgr), QuantifierTreeNode(qvmgr), Formula(mgr, qvmgr) 
                            {
                            }

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator)
                            : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), Formula(mgr, qvManipulator)
                            {
                            }

void QuantifierTreeFormula::localise() {
    removeUnusedVars();
}

QuantifierTreeFormula* QuantifierTreeFormula::changeToFormula(Cudd &) {
    return this;
}

void QuantifierTreeFormula::negate() {
    setMatrix(!getMatrix());
}

// computes supportSet and uVarsSupportSet
void QuantifierTreeFormula::computeSupportSets() {
    if (needsToRecomputeSupportSet()) { // only if we need to compute it though
        // this computes supportSet
        Formula::getSupportSet();
        // this computes uVarsSupportSet
        uVarsSupportSet.clear();
        for (const Variable &var : supportSet) {
            if (isVarUniv(var)) {
                uVarsSupportSet.insert(var);
            } else if (isVarExist(var)) {
                uVarsSupportSet.insert(getExistVarDependencies(var).begin(),
                                       getExistVarDependencies(var).end());
            }
        }
    }
}

VariableSet const& QuantifierTreeFormula::getSupportSet() {
    computeSupportSets();
    return supportSet;
}

VariableSet const& QuantifierTreeFormula::getUVarsSupportSet() {
    computeSupportSets();
    return uVarsSupportSet;
}


void QuantifierTreeFormula::addToUVarsOutsideThisSubtree(const VariableSet &varsToAdd) {
    uVarsOutsideThisSubtree.insert(varsToAdd.begin(), varsToAdd.end());
}