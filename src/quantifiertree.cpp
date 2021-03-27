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

#include<iostream>

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

// VariableSet const& QuantifierTreeNode::getUVarsOutsideThisSubtree() const {
//     return uVarsOutsideThisSubtree;
// }

VariableSet const& QuantifierTreeNode::getUVarsSupportSet() {
    return uVarsSupportSet;
}

// void QuantifierTreeNode::addToUVarsOutsideThisSubtree(const Variable &varToAdd) {
//     addToUVarsOutsideThisSubtree(VariableSet{varToAdd});
// }

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

// void QuantifierTree::addPossibleUnivVarsToMapping(std::unordered_map<Variable, std::list<QuantifierTreeNode*>> &childrenToCombineMapping) {
//     //std::cout << "Finding new univ variables...." << std::endl;
    
//     // get all universal variables on which some leftover exist var depends (these vars cant be pushed)...
//     VariableSet dependentUnivVars; 
//     for (const Variable &existVar : getExistVars()) {
//         dependentUnivVars.insert(getExistVarDependencies(existVar).begin(), getExistVarDependencies(existVar).end());
//     }
//     // ...add to it variables which are already in mapping...
//     for (const auto &univVarInMappingWithMappedContent : childrenToCombineMapping) {
//         dependentUnivVars.insert(univVarInMappingWithMappedContent.first);
//     }
//     // ...and take all universal variables which are not in it
//     VariableSet notDependentUnivVars = getUnivVars().minus(dependentUnivVars);

//     // for each universal variable that is not in dependentUnivVars get the list of children containing it or its dependency
//     if (!notDependentUnivVars.empty()) {
//         for (QuantifierTreeNode *child : children) {
//             const VariableSet &childSupportSet = child->getSupportSet();
//             for (const Variable &uVar : notDependentUnivVars) {
//                 // if the child contains uVar...
//                 if (childSupportSet.contains(uVar)) {
//                     childrenToCombineMapping[uVar].push_back(child);
//                     continue;
//                 }
//                 // ... or some exist var that depends on uVar
//                 for (const Variable &eVar : getUnivVarDependencies(uVar)) {
//                     if (childSupportSet.contains(eVar)) {
//                         childrenToCombineMapping[uVar].push_back(child);
//                         continue;
//                     }
//                 }
//             }
//         }
//     }

//     //std::cout << "....found" << std::endl;
// }

void QuantifierTree::addExistVarsToChildrenToCombineMapping(const VariableSet &eVarsToAdd) {
    for (QuantifierTreeNode *child : children) {
        const VariableSet &childSupportSet = child->getSupportSet();
        for (const Variable &eVar : eVarsToAdd) {
            if (childSupportSet.contains(eVar)) {
                // we add only those children which contain eVar
                childrenToCombineMapping[eVar].push_back(child);
            }
        }
    }
}

void QuantifierTree::addUnivVarsToChildrenToCombineMapping(const VariableSet &uVarsToAdd) {
    for (QuantifierTreeNode *child : children) {
        const VariableSet &childUVarsSupportSet = child->getUVarsSupportSet();
        for (const Variable &uVar : uVarsToAdd) {
            if (childUVarsSupportSet.contains(uVar)) {
                // we add only those children which either contain uVar,
                // or contain some eVar that depends on uVar
                childrenToCombineMapping[uVar].push_back(child);
            }
        }
    }
}

void QuantifierTree::pushVarsWithCombining() {//std::unordered_map<Variable, std::list<QuantifierTreeNode*>> &childrenToCombineMapping, bool findNewUnivVars) {
    //std::cout << "We are planning to remove these with combining" << std::endl;
    //for (auto i : childrenToCombineMapping) {
    //    std::cout << "  variable " << i.first << " with " << i.second.size() << " children to combine." << std::endl;
    //}

    // for each universal variable from getUnivVars(), this map keeps the number
    // of dependencies in dependency sets of exist vars from getExistVars()
    std::unordered_map<Variable, unsigned> numOfDependenciesForUnivVars;
    // the set of universal variables that can be pushed, i.e. on which no exist var depends 
    VariableSet pushableUnivVars = getUnivVars();
    for (const Variable &existVar : getExistVars()) {
        for (const Variable &dependentUnivVar : getExistVarDependencies(existVar)) {
            ++numOfDependenciesForUnivVars[dependentUnivVar];
            pushableUnivVars.erase(dependentUnivVar);
        }
    }

    // a heplful function which will either push or prepare to push univ vars based on whether we have conjunction or disjunction
    auto pushOrPrepareToPushUnivVars = [&] {
        if (isConj) {
            // for conjunction, we can just push all pushable univ vars to each child
            for (const Variable &univVarToPush : pushableUnivVars) {
                // push it into every child, pushUnivVar() will take care of deciding whether to keep it there or not
                for (QuantifierTreeNode *child : children) {
                    //std::cout << "Pushing " << univVarToPush << " to child " << *child << std::endl;
                    //std::cout << child->getSupportSet() << std::endl;
                    child->pushUnivVar(univVarToPush);
                }
                removeUnivVar(univVarToPush);
            }
        } else {
            // for disjunction, we can only push univ vars to combined child, 
            // therefore we add them to the mapping and we will push them later
            addUnivVarsToChildrenToCombineMapping(pushableUnivVars);
        }
        // we can clear pushableUnivVars because we processed them
        pushableUnivVars.clear();
    };
    // we also run this helpful function now, it will push (or prepare to push) all univ vars on which nothing depends
    pushOrPrepareToPushUnivVars();
    
    // we are going to push each variable v to the children containing it (which are saved in childrenToCombineMapping[v])
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

        // check if there are some new univ vars that can become pushable after we will push varToPush
        if (isVarExist(varToPush)) {
            for (const Variable &dependentUnivVar : getExistVarDependencies(varToPush)) {
                --numOfDependenciesForUnivVars[dependentUnivVar];
                if (numOfDependenciesForUnivVars[dependentUnivVar] == 0) {
                    // after pushing varToPush, no exist var will be dependent on 
                    // dependentUnivVar in this node of the quantifier tree
                    pushableUnivVars.insert(dependentUnivVar);
                }
            }
        }

        // if varToPush is needed to be pushed to only one child, there is no need to create new child
        if (childrenToCombine.size() == 1) {
            if (isVarUniv(varToPush)) {
                (*childrenToCombine.begin())->pushUnivVar(varToPush);
                removeUnivVar(varToPush);
            } else if (isVarExist(varToPush)) {
                (*childrenToCombine.begin())->pushExistVar(varToPush);
                removeExistVar(varToPush);
                // for exist var, we need to push (or prepare to push) univ vars which became pushable
                pushOrPrepareToPushUnivVars();
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
        for (const Variable &varWithChangedChildren : varsWithChangedChildren) {
            childrenToCombineMapping[varWithChangedChildren].push_back(newChild);
        }

        childrenToCombineMapping.erase(varToPush);

        // for existential variable, we need to push or prepare to push univ vars which became pushable
        if (isVarExist(varToPush)) {
            pushOrPrepareToPushUnivVars();
        }
    }
}

VariableSet &QuantifierTree::getUVarsOutsideChildSubtree(QuantifierTreeNode* child, std::list<QuantifierTreeNode*> originalChildren, const VariableSet &uVarsOutsideThisSubtree) {
    if (uVarsOutsideChildSubtree.find(child) == uVarsOutsideChildSubtree.end()) { // if uVars outside child was not computed yet
        
        // the uVars outside child should be those that are already outside this subtree (reduced only to uVars actually occuring in the child)
        const VariableSet &uVarsInChildSupportSet = child->getUVarsSupportSet();
        uVarsOutsideChildSubtree[child] = uVarsOutsideThisSubtree.intersect(uVarsInChildSupportSet);
        
        // and furthermore, it should be also those, that are in (original) siblings of child (again reduced to uVars actually occuring in the child)
        for (auto originalChild : originalChildren) {
            if (originalChild != child) {
                VariableSet uVarsToAdd = originalChild->getUVarsSupportSet().intersect(uVarsInChildSupportSet);
                uVarsOutsideChildSubtree[child].insert(uVarsToAdd.begin(), uVarsToAdd.end());
            }
        }

    }

    return uVarsOutsideChildSubtree[child];
}

void QuantifierTree::pushExistVarsSeparately(const VariableSet &uVarsOutsideThisSubtree) {

    if (isConj) {
        throw DQBDDexception("Trying to push existential variables as in disjunction for a conjuction node");
    }

    /**
     * to be able to push exist. var y with dependency set D_y to each child separately:
     *    (1) all children containing y have (uVarsSupportSet / D_y) pairwise disjoint
     *    (2) at most one child containing y has some variable from (uVarsSupportSet / D_y) in uVarsSupportSet
     *        of child that does not contain y or in uVarsOutsideThisSubtree
     **/
    VariableSet existVars = getExistVars();
    for (const Variable &existVarToPush : existVars) {
        // if true, existVarToPush fulfills the conditions
        bool canBePushed = true;
        // true if one of the processed children was already problematic
        bool wasThereOneProblem = false;

        // first check the condition
        for (QuantifierTreeNode *child : childrenToCombineMapping[existVarToPush]) {
            /* uVarsToCheck
             *   - the set of universal variables occuring outside the subtree rooted in child intersected with the 
             *     actual universal variables occuring in this subtree (either directly or in the dependency set of some
             *     exist var occuring in the subtree) -- this is in getUVarsOutsideChildSubtree[child] -- and from these 
             *     only those which are not in the dependency set of existVarToPush
             * 
             * also we can call getUVarsOutsideChildSubtree with children, because this function should be called before combining
             */
            VariableSet uVarsToCheck = getUVarsOutsideChildSubtree(child, children, uVarsOutsideThisSubtree).minus(getExistVarDependencies(existVarToPush));

            /* check if uVarsToCheck is empty
                *     - if it violates condition (1), then it is true and it either sets wasThereOneProblem
                *       to true and some next child with which it intersects will also evaluate this to true
                *       but wasThereOneProblem will be already set to true OR this is the second child and 
                *       wasThereOneProblem is already true
                *     - if it violates condition (2), then it is true and we set wasThereOneProblem to true
                *       so we know this one child is the one possible child that can violate (2), or this 
                *       could be a second child violating (2) and then wasThereOneProblem is already true
                *       and we cannot push it
                *     - if it does not violate either, then this is false 
                */
            if ( !(uVarsToCheck.empty()) ) {
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
            for (QuantifierTreeNode *child : childrenToCombineMapping[existVarToPush]) {
                /* There is no need to rename them, because universal vars on which this exist var 
                * depends will either stay as ancestor/same level or will be pushed into some sibling 
                * and then the dependency will be removed. Also, pushExistVar() takes care whether
                * to actually push it inside or not.       
                */
                child->pushExistVar(existVarToPush);
            }
            removeExistVar(existVarToPush);
            childrenToCombineMapping.erase(existVarToPush);
        }
    }
}



void QuantifierTree::localise(const VariableSet &uVarsOutsideThisSubtree) {
    /**
     * Localisation for conjunction and disjunction works similarly:
     *  1) disjunction starts with attempting to push existential quantifiers into children separately
     *  2) both use childrenToCombineMapping and the related functions to push quantifiers into new
     *     children created by combining only those children that contain the variable that is being
     *     pushed where:
     *        - conjunction pushes this way only existential variables, universal variables are pushed
     *          separately to original and newly created children (as deep as possible based on the
     *          dependencies, i.e. if y depends on x and y is in some child, then x can be pushed only
     *          into this child and not further)
     *        - disjunction pushes this way the leftover existential variables which were not pushed
     *          separately in the first step, and all the possible universal variables (again as 
     *          deep as possible as explained in the previous remark)
     *  3) calling localisation on the original children (i.e. not those created by previous step)
     *     because the previous step should assure that everything that is possible to be pushed
     *     was pushed as deep as it could (deep means between this node through the new children
     *     to the original children)
     * 
     * Example:
     * If the children are ch_1, ch_2, ch_3, ch_4 then:
     *   - If we have disjunction, we check the conditions for each existential variable whether
     *     it can be pushed separately and push them separately into ch_1, ch_2, ch_3, and ch_4.
     *   - We call pushVarsWithCombining() which will do the second step, where for example if
     *     y was not pushed in previous step and appears only in ch_2 and ch_3, then we end up with
     *     new child ch_y which will have as children ch_2 and ch_3, and this node will then have
     *     children ch_1, ch_4, and ch_y. This happens iteratively, so we can get another new child
     *     ch_x, created for variable x which appears in ch_1 and ch_y and then ch_x will have children
     *     ch_1 and ch_y, and this node will updates its children as ch_4 and ch_x. During this,
     *     all variables are pushed as deep as possible till the original children, so all pushed variables
     *     appear somewhere between this node and the original children ch_1, ch_2... (or also in the 
     *     original children).
     *   - Localisation is then called for the original children ch_1, ch_2, ch_3, and ch_4 where we need
     *     to update uVarsOutsideThisSubtree for each child.
     */

    removeUnusedVars();
    // if we do not have any quantifiers to push, we are done with localising for this subtree
    if (getUnivVars().empty() && getExistVars().empty()) {
        return;
    }

    //std::cout << "Localising at tree: " << *this << " with: " << std::endl
    //          << "   support set:   " << getSupportSet() << std::endl
    //          << "   u support set: " << getUVarsSupportSet() << std::endl
    //          << "   u outside:     " << uVarsOutsideThisSubtree << std::endl;

    // for each exist var y compute the children that contain y
    // and save them in childrenToCombineMapping[y]
    addExistVarsToChildrenToCombineMapping(getExistVars());

    // do the 1) step
    if (!isConj) {
        pushExistVarsSeparately(uVarsOutsideThisSubtree);
    }

    // make a copy of original children to call localisation only for them at the end
    std::list<QuantifierTreeNode*> originalChildren = children;

    // do the 2) step 
    pushVarsWithCombining();

    // recursively call localise for each (original) child - step 3)
    for (QuantifierTreeNode *child : originalChildren) {
        // TODO check if the subtree contains any quantifiers, if not we can skip and not compute getUVarsOutsideChildSubtree for the child (important for conjunction, for disjunction it is computed)

        VariableSet &uVarsOutsideThisChildSubtree = getUVarsOutsideChildSubtree(child, originalChildren, uVarsOutsideThisSubtree);
        
        /* TODO we need to check all copies of univ vars outside childSubtree, whether they are really outside or they are copies. 
         * This next is just dirty not complete fix (the fix just removes universal variables quantified in the child, because 
         * obviously the occurences of these variables outside subtree will be just copies which we did not rename).
         */
        child->localise(uVarsOutsideThisChildSubtree.minus(child->getUnivVars()));
        
        // to save on memory, we can delete uVarsOutsideChildSubtree[child] because we do not use it further
        uVarsOutsideChildSubtree.erase(child);
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


// void QuantifierTree::addToUVarsOutsideThisSubtree(const VariableSet &varsToAdd) {
//     for (auto child : children) {
//         child->addToUVarsOutsideThisSubtree(varsToAdd);
//     }

//     uVarsOutsideThisSubtree.insert(varsToAdd.begin(), varsToAdd.end());
// }

void QuantifierTree::addChild(QuantifierTreeNode *child) {
    /* in oldUVarsSupportSet are all universal variables in all the children
     * of this tree before adding new child, therefore all those will be
     * added to uVarsOutsideThisSubtree of the new child
     */
    // VariableSet oldUVarsSupportSet = getUVarsSupportSet();

    // add variables of the child to the support sets here
    supportSet.insert(child->getSupportSet().begin(),
                    child->getSupportSet().end());
    uVarsSupportSet.insert(child->getUVarsSupportSet().begin(),
                        child->getUVarsSupportSet().end());

    // update others children uVarsOutsideThisSubtree with the universal
    // variables that are in support set of what will be newly added child
    // for (auto oldChild : children) {
    //     oldChild->addToUVarsOutsideThisSubtree(child->getUVarsSupportSet());
    // }

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
            //childOfChild->addToUVarsOutsideThisSubtree(oldUVarsSupportSet);
        }

        // and finally delete it
        treeChild->children.clear();
        delete child;
    } else {
        // otherwise just add this child to children
        children.push_back(child);
        //child->addToUVarsOutsideThisSubtree(oldUVarsSupportSet);
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

void QuantifierTreeFormula::localise(const VariableSet &uVarsOutsideThisSubtree) {
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


// void QuantifierTreeFormula::addToUVarsOutsideThisSubtree(const VariableSet &varsToAdd) {
//     uVarsOutsideThisSubtree.insert(varsToAdd.begin(), varsToAdd.end());
// }