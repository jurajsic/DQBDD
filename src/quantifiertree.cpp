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

#include <algorithm>

#include "quantifiertree.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

QuantifierTreeNode::QuantifierTreeNode(QuantifiedVariablesManager &qvmgr) : QuantifiedVariablesManipulator(qvmgr) {}

void QuantifierTreeNode::pushExistVar(const Variable &var) {
    // only put var in this node if it is used in this node
    if (getSupportSet().contains(var)) {
        addExistVar(var);
    }
}

void QuantifierTreeNode::pushUnivVar(const Variable &var) {
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

VariableSet const& QuantifierTreeNode::getUVarsSupportSet() {
    return uVarsSupportSet;
}

/*******************************************/
/*******************************************/
/******* QUANTIFIER TREE CONNECTION ********/
/*******************************************/
/*******************************************/ 

QuantifierTreeConnection::QuantifierTreeConnection(QuantifierTreeNode *child) : child(child) {}

void QuantifierTreeConnection::decreaseNumOfParents() {
    --numOfParents;
}

void QuantifierTreeConnection::increaseNumOfParents() {
    ++numOfParents;
}

unsigned QuantifierTreeConnection::getNumOfParents() {
    return numOfParents;
}

/*******************************************/
/*******************************************/
/************ QUANTIFIER TREE **************/
/*******************************************/
/*******************************************/ 


QuantifierTreeFormula* QuantifierTreeConnection::changeChildToQuantifierTreeFormula(Cudd &mgr) {
    QuantifierTreeFormula* childFormula = child->changeToFormula(mgr);
    child = childFormula;
    return childFormula;
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeConnection*> childrenConnections, QuantifiedVariablesManager &qvMgr, bool collapseChildren, bool changeNumOfParents) : QuantifiedVariablesManipulator(qvMgr), QuantifierTreeNode(qvMgr), isConj(isConj) {
    supportSet = {};
    uVarsSupportSet = {};
    if (childrenConnections.size() < 2) {
        throw dqbddException((std::string("You cannot create a quantifier tree with ") + std::to_string(childrenConnections.size()) + std::string(" operands")).c_str());
    }
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        addChild(childConnection, collapseChildren, changeNumOfParents);
    }
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeConnection*> childrenConnections, QuantifiedVariablesManipulator &qvManipulator, bool collapseChildren) : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), isConj(isConj) {
    supportSet = {};
    uVarsSupportSet = {};
    if (childrenConnections.size() < 2) {
        throw dqbddException((std::string("You cannot create a quantifier tree with ") + std::to_string(childrenConnections.size()) + std::string(" operands")).c_str());
    }
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        addChild(childConnection, collapseChildren);
    }
}

QuantifierTree::~QuantifierTree() {
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        deleteChildConnection(childConnection);
    }
}

void QuantifierTree::deleteChildConnection(QuantifierTreeConnection *childConnection) {
    childConnection->decreaseNumOfParents();
    if (childConnection->getNumOfParents() == 0) {
        delete childConnection->child;
        delete childConnection;
    }
}

bool QuantifierTree::removeFromOrderedListOtherOrderedListUsingChildrenOrder(std::list<QuantifierTreeConnection*> &listToRemoveFrom, std::list<QuantifierTreeConnection*> &listOfItemsToRemove) {
    //std::cout << "Deleting from list of size " << listToRemoveFrom.size() << " the list of size " << listOfItemsToRemove.size() << " using children order of size " << children.size() << std::endl;
    
    auto childrenConnectionsCopy = childrenConnections;
    auto orderIter = childrenConnectionsCopy.begin();
    auto listToRemoveFromIter = listToRemoveFrom.begin();
    auto listOfItemsToRemoveIter = listOfItemsToRemove.begin();

    bool somethingWasDeleted = false;

    while (orderIter != childrenConnectionsCopy.end() && listToRemoveFromIter != listToRemoveFrom.end() && listOfItemsToRemoveIter != listOfItemsToRemove.end()) {
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

void QuantifierTree::addExistVarsToChildrenToCombineMapping(const VariableSet &eVarsToAdd) {
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        //QuantifierTreeNode *child = childConnection->child;
        const VariableSet &childSupportSet = childConnection->child->getSupportSet();
        for (const Variable &eVar : eVarsToAdd) {
            if (childSupportSet.contains(eVar)) {
                // we add only those children which contain eVar
                childrenToCombineMapping[eVar].push_back(childConnection);
            }
        }
    }
}

void QuantifierTree::addUnivVarsToChildrenToCombineMapping(const VariableSet &uVarsToAdd) {
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        const VariableSet &childUVarsSupportSet = childConnection->child->getUVarsSupportSet();
        for (const Variable &uVar : uVarsToAdd) {
            if (childUVarsSupportSet.contains(uVar)) {
                // we add only those children which either contain uVar,
                // or contain some eVar that depends on uVar
                childrenToCombineMapping[uVar].push_back(childConnection);
            }
        }
    }
}

void QuantifierTree::pushVarsWithCombining() {
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
            if (getUnivVars().contains(dependentUnivVar)) {
                ++numOfDependenciesForUnivVars[dependentUnivVar];
                pushableUnivVars.erase(dependentUnivVar);
            }
        }
    }

    // a heplful function which will either push or prepare to push univ vars based on whether we have conjunction or disjunction
    auto pushOrPrepareToPushUnivVars = [this, &pushableUnivVars] {
        if (this->isConj) {
            // for conjunction, we can just push all pushable univ vars to each child
            for (const Variable &univVarToPush : pushableUnivVars) {
                // push it into every child, pushUnivVar() will take care of deciding whether to keep it there or not
                for (QuantifierTreeConnection *childConnection : this->childrenConnections) {
                    //std::cout << "Pushing " << univVarToPush << " to child " << *child << std::endl;
                    //std::cout << child->getSupportSet() << std::endl;
                    childConnection->child->pushUnivVar(univVarToPush);
                }
                this->removeUnivVar(univVarToPush);
            }
        } else {
            // for disjunction, we can only push univ vars to combined child, 
            // therefore we add them to the mapping and we will push them later
            this->addUnivVarsToChildrenToCombineMapping(pushableUnivVars);
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
                    [](const std::pair<Variable, std::list<QuantifierTreeConnection*>> &a, 
                        const std::pair<Variable, std::list<QuantifierTreeConnection*>> &b) { 
                        return (a.second.size() < b.second.size());
                    }
                    );
        Variable varToPush = varToPushAndItsChildren->first;
        //std::cout << "Pushing the " << (isVarUniv(varToPush) ? "universal" : "existential") << " variable " << varToPush << " to a new combined child." << std::endl;
        std::list<QuantifierTreeConnection*> childrenToCombine = varToPushAndItsChildren->second;
        //std::cout << "The size of children to combine right after getting it: " << childrenToCombine.size() << std::endl;

        // if we hit the variable that is in all children (meaning that all other variables are in all children too) we stop
        if (childrenToCombine.size() == childrenConnections.size()) {
            break;
        }

        // check if there are some new univ vars that can become pushable after we will push varToPush
        if (isVarExist(varToPush)) {
            for (const Variable &dependentUnivVar : getExistVarDependencies(varToPush)) {
                if (getUnivVars().contains(dependentUnivVar)) {
                    --numOfDependenciesForUnivVars[dependentUnivVar];
                    if (numOfDependenciesForUnivVars[dependentUnivVar] == 0) {
                        // after pushing varToPush, no exist var will be dependent on 
                        // dependentUnivVar in this node of the quantifier tree
                        pushableUnivVars.insert(dependentUnivVar);
                    }
                }
            }
        }

        // if varToPush is needed to be pushed to only one child, there is no need to create new child
        if (childrenToCombine.size() == 1) {
            if (isVarUniv(varToPush)) {
                (*childrenToCombine.begin())->child->pushUnivVar(varToPush);
                removeUnivVar(varToPush);
            } else if (isVarExist(varToPush)) {
                (*childrenToCombine.begin())->child->pushExistVar(varToPush);
                removeExistVar(varToPush);
                // for exist var, we need to push (or prepare to push) univ vars which became pushable
                pushOrPrepareToPushUnivVars();
            } else {
                throw dqbddException("Found variable which is neither universal or existential, this should not happen");
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
        removeFromOrderedListOtherOrderedListUsingChildrenOrder(childrenConnections, childrenToCombine);
        //std::cout << "The size of children to combine right before creating new child: " << childrenToCombine.size() << std::endl;
        QuantifierTree *newChild = new QuantifierTree(isConj, childrenToCombine, *qvMgr, false, false);
        // as nothing is possible to be pushed further in the new child, we can set this to false
        newChild->needsToLocalise = false;
        if (isVarExist(varToPush)) {
            newChild->pushExistVar(varToPush);
            removeExistVar(varToPush);
        } else if (isVarUniv(varToPush)) {
            newChild->pushUnivVar(varToPush);
            removeUnivVar(varToPush);
        } else {
            throw dqbddException("Found variable which is neither universal or existential, this should not happen");
        }
        QuantifierTreeConnection *newChildConnection = new QuantifierTreeConnection(newChild);
        childrenConnections.push_back(newChildConnection);
        newChildConnection->increaseNumOfParents();

        
        //std::cout << "...combined tree created" << std::endl;

        // add the new child to each var in varsWithChangedChildren
        for (const Variable &varWithChangedChildren : varsWithChangedChildren) {
            childrenToCombineMapping[varWithChangedChildren].push_back(newChildConnection);
        }

        childrenToCombineMapping.erase(varToPush);

        // for existential variable, we need to push or prepare to push univ vars which became pushable
        if (isVarExist(varToPush)) {
            pushOrPrepareToPushUnivVars();
        }
    }
}

VariableSet &QuantifierTree::getUVarsOutsideChildSubtree(QuantifierTreeConnection* childConnection, const std::list<QuantifierTreeConnection*> &childAndSiblingsConnections, const VariableSet &uVarsOutsideThisSubtree) {
    if (uVarsOutsideChildSubtree.find(childConnection) == uVarsOutsideChildSubtree.end()) { // if uVars outside child was not computed yet
        
        // the uVars outside child should be those that are already outside this subtree (reduced only to uVars actually occuring in the child)
        const VariableSet &uVarsInChildSupportSet = childConnection->child->getUVarsSupportSet();
        uVarsOutsideChildSubtree[childConnection] = uVarsOutsideThisSubtree.intersect(uVarsInChildSupportSet);
        
        // and furthermore, it should be also those, that are in (original) siblings of child (again reduced to uVars actually occuring in the child)
        for (QuantifierTreeConnection *childOrSiblingConnection : childAndSiblingsConnections) {
            if (childOrSiblingConnection != childConnection) {
                VariableSet uVarsToAdd = childOrSiblingConnection->child->getUVarsSupportSet().intersect(uVarsInChildSupportSet);
                uVarsOutsideChildSubtree[childConnection].insert(uVarsToAdd.begin(), uVarsToAdd.end());
            }
        }

    }

    return uVarsOutsideChildSubtree[childConnection];
}

void QuantifierTree::pushExistVarsSeparately(const VariableSet &uVarsOutsideThisSubtree) {

    if (isConj) {
        throw dqbddException("Trying to push existential variables as in disjunction for a conjuction node");
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
        for (QuantifierTreeConnection *childConnection : childrenToCombineMapping[existVarToPush]) {
            /* uVarsToCheck
             *   - the set of universal variables occuring outside the subtree rooted in child intersected with the 
             *     actual universal variables occuring in this subtree (either directly or in the dependency set of some
             *     exist var occuring in the subtree) -- this is in getUVarsOutsideChildSubtree[child] -- and from these 
             *     only those which are not in the dependency set of existVarToPush
             * 
             * also we can call getUVarsOutsideChildSubtree with children, because this function should be called before combining
             */
            const VariableSet &uVarsToCheck = getUVarsOutsideChildSubtree(childConnection, childrenConnections, uVarsOutsideThisSubtree).minus(getExistVarDependencies(existVarToPush));

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
            for (QuantifierTreeConnection *childConnection : childrenToCombineMapping[existVarToPush]) {
                /* There is no need to rename them, because universal vars on which this exist var 
                * depends will either stay as ancestor/same level or will be pushed into some sibling 
                * and then the dependency will be removed. Also, pushExistVar() takes care whether
                * to actually push it inside or not.       
                */
                childConnection->child->pushExistVar(existVarToPush);
            }
            removeExistVar(existVarToPush);
            childrenToCombineMapping.erase(existVarToPush);
        }
    }
}



void QuantifierTree::localise(const VariableSet &uVarsOutsideThisSubtree) {
    /**
     * Localisation for conjunction and disjunction works similarly:
     *  1) Disjunction starts with attempting to push existential quantifiers into children separately
     *  2) Both use childrenToCombineMapping and the related functions to push quantifiers into new
     *     children created by combining only those children that contain the variable that is being
     *     pushed where:
     *        - conjunction pushes this way only existential variables, universal variables are pushed
     *          separately to original and newly created children (as deep as possible based on the
     *          dependencies, i.e. if y depends on x and y is in some child, then x can be pushed only
     *          into this child and not further)
     *        - disjunction pushes this way the leftover existential variables which were not pushed
     *          separately in the first step, and all the possible universal variables (again as 
     *          deep as possible as explained in the previous remark)
     *  3) Localising only on the original children (i.e. not those created by previous step)
     *     because the previous step should assure that everything that is possible to be pushed
     *     was pushed as deep as it could (deep means between this node through the new children
     *     to the original children). However localise() function is called for conjunction on
     *     the new children, but for these, localisation is turned off, we only use that function
     *     to pass on uVarsOutsideSubtree in such a way, that copies of universal variables (because
     *     we do not rename variables, some copies can have same name) occuring outside each subtree 
     *     are not passed along.
     * 
     * Example:
     * If the children are ch_1, ch_2, ch_3, ch_4 then:
     *   - If we have disjunction, we check the conditions for each existential variable whether
     *     it can be pushed separately and push them separately into ch_1, ch_2, ch_3, and ch_4.
     *   - We call pushVarsWithCombining() which will do the second step, where for example if
     *     y was not pushed in previous step and appears only in ch_2 and ch_3, then we end up with
     *     new child ch_y which will have as children ch_2 and ch_3, and the current node will then have
     *     children ch_1, ch_4, and ch_y. This happens iteratively, so we can get another new child
     *     ch_x, created for variable x which appears in ch_1 and ch_y and then ch_x will have children
     *     ch_1 and ch_y, and this node will update its children to ch_4 and ch_x. During this,
     *     all variables are pushed as deep as possible between this node and the original children (included), 
     *     so all pushed variables appear somewhere between this node and the original children ch_1, ch_2, 
     *     ch_3, ch_4 (or also in the original children).
     *   - Localisation is then called for the original children ch_1, ch_2, ch_3, and ch_4 for disjunction
     *     and for ch_4 and ch_x (i.e. the new children) for conjunction. However, as quantifiers are pushed
     *     as deeply as possible, proper localisation will occur also only for ch_1, ch_2, ch_3, and ch_4, 
     *     localise() for ch_x and ch_y will only be used, to compute uVarsOutsideThisSubtree in such a way
     *     that that copies of universal variables occuring outside each subtree are not passed along.
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

    // replace all children that have multiple parents with a new copy, so that this quantifier tree is their only parent (so we can push quantifiers without conflicts)
    for (auto childrenConnectionsIter = childrenConnections.begin(); childrenConnectionsIter != childrenConnections.end(); ++childrenConnectionsIter) {
        QuantifierTreeConnection *currentChildConnection = *childrenConnectionsIter;
        if (currentChildConnection->getNumOfParents() > 1) {
            QuantifierTreeConnection *newChildConnection = new QuantifierTreeConnection(currentChildConnection->child->getCopy());
            (*childrenConnectionsIter) = newChildConnection;
            currentChildConnection->decreaseNumOfParents();
            newChildConnection->increaseNumOfParents();
        }
    }

    if (needsToLocalise) {
        // for each exist var y compute the children that contain y
        // and save them in childrenToCombineMapping[y]
        addExistVarsToChildrenToCombineMapping(getExistVars());

        // do the 1) step
        if (!isConj) {
            pushExistVarsSeparately(uVarsOutsideThisSubtree);
        }
    }

    // for disjunction, localisation will be called only on original children (and we need a copy here, before pushVarsWithCombining() changes children)
    std::list<QuantifierTreeConnection*> childrenToCallLocaliseOn;
    if (!isConj) {
        childrenToCallLocaliseOn = childrenConnections;
    }

    if (needsToLocalise) {
        // do the 2) step 
        pushVarsWithCombining();
    }

    // for conjunction, we want to call localisation on the newly created children, on which localisation 
    // will not occur (as quantifiers should be already maximally pushed for them), however, it will
    // compute uVarsOutside them in such a way, to exclude univ vars which are just non-renamed 'copies'
    if (isConj) {
        childrenToCallLocaliseOn = childrenConnections;
    }

    // recursively call localise for each (original) child - step 3)
    for (QuantifierTreeConnection *childConnection : childrenToCallLocaliseOn) {
        QuantifierTreeNode *child = childConnection->child;
        // if child does not have any quantifier to push, we can skip localisation, possibly saving on computing uVarsOutsideChildSubtree (for conjunction)
        if (child->getUnivVars().empty() && child->getExistVars().empty()) {
            // to save on memory, we can delete uVarsOutsideChildSubtree[child] because we do not use it further (here only important for disjunction)
            uVarsOutsideChildSubtree.erase(childConnection);
            continue;
        }

        VariableSet &uVarsOutsideThisChildSubtree = getUVarsOutsideChildSubtree(childConnection, childrenToCallLocaliseOn, uVarsOutsideThisSubtree);
        
        /* This part is important probably only for conjunction: 
         * To not pass along non-renamed copies of univ vars, we remove those universal quantifiers, which are quantified 
         * in the child. This removes all possible copies from siblings and because we assume, that localise() for this node
         * was called without copies, this should accomplish it. Also, for conjunction, we are doing this also for newly created
         * children, so it should really do it (for disjunction it is not important, no univ vars copies can be created, that is
         * why we can just skip to original children).
         */
        child->localise(uVarsOutsideThisChildSubtree.minus(child->getUnivVars()));
        
        // to save on memory, we can delete uVarsOutsideChildSubtree[child] because we do not use it further
        uVarsOutsideChildSubtree.erase(childConnection);
    }
}

QuantifierTreeFormula* QuantifierTree::changeToFormula(Cudd &mgr) {
    BDD matrix;
    if (isConj) {
        matrix = mgr.bddOne();
    } else {
        matrix = mgr.bddZero();
    }

    auto childConnectionIter = childrenConnections.begin();
    while (childConnectionIter != childrenConnections.end()) {
        QuantifierTreeConnection *childConnection = *childConnectionIter;

        // remove the child from the list of children
        childConnectionIter = childrenConnections.erase(childConnectionIter);

        // get the formula from child
        QuantifierTreeFormula *childFormula = childConnection->changeChildToQuantifierTreeFormula(mgr);

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
            throw dqbddException("Selected choice of variables to eliminate in tree is not supported.");
            break;
        }
        }
        
        // check if we can stop (i.e. the resulting formula after adding the child formula is simple)
        if (isConj) {
            matrix &= childFormula->getMatrix();
            if (matrix.IsZero()) {
                deleteChildConnection(childConnection);
                break;
            }
        } else {
            matrix |= childFormula->getMatrix();
            if (matrix.IsOne()) {
                deleteChildConnection(childConnection);
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

        /* here we assume that if child had any quantifiers in prefix, then this tree was its only parent,
         * therefore, deleteChildConnection will actually delete the child and its prefix
         */
        deleteChildConnection(childConnection);
    }
    
    // delete leftover children
    while (childConnectionIter != childrenConnections.end()) {
        deleteChildConnection(*childConnectionIter);
        childConnectionIter = childrenConnections.erase(childConnectionIter);
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

void QuantifierTree::addChild(QuantifierTreeConnection *childConnection, bool collapseChildren, bool changeNumOfParents) {
    QuantifierTreeNode *child = childConnection->child;

    // add variables of the child to the support sets here
    supportSet.insert(child->getSupportSet().begin(),
                    child->getSupportSet().end());
    uVarsSupportSet.insert(child->getUVarsSupportSet().begin(),
                        child->getUVarsSupportSet().end());

    // collapsing children
    auto treeChild = dynamic_cast<QuantifierTree*>(child);
    if (collapseChildren                        // if we want to collapse children...
            && treeChild != nullptr             // ...and the child is QuantifierTree...
            && treeChild->isConj == isConj      // ...with the same operation...
            // ...and with empty quantifier prefix...
            && treeChild->getExistVars().empty() && treeChild->getUnivVars().empty()) {

        // ...then we collapse children, i.e. set its children as current children, not itself...
        for (QuantifierTreeConnection *childOfChildConnection : treeChild->childrenConnections) {
            childrenConnections.push_back(childOfChildConnection);
            if (changeNumOfParents) {
                childOfChildConnection->increaseNumOfParents();
            }
        }

        // ...and finally delete it
        //treeChild->children.clear();
        //delete child;
    } else {
        // otherwise just add this child to children
        childrenConnections.push_back(childConnection);
        if (changeNumOfParents) {
            childConnection->increaseNumOfParents();
        }
    }
}

std::ostream& QuantifierTree::print(std::ostream& out) const {
    std::string op;
    if (isConj) {
        op = '&';
    } else {
        op = '|';
    }

    auto size = childrenConnections.size();
    out << std::string("(");
    for (QuantifierTreeConnection *childConnection : childrenConnections) {
        out << *(childConnection->child);
        // print the last element without operator
        if (size != 1) {
            out << std::string(" ") << op << std::string(" ");
        }
        --size;
    }
    out << std::string(")");
    return out;
}

QuantifierTreeNode* QuantifierTree::getCopy() {
    if (!(getUnivVars().empty() && getExistVars().empty())) {
        throw dqbddException("Trying to copy quantifier tree with some variables in prefix");
    }

    QuantifierTree *newTree = new QuantifierTree(isConj, childrenConnections, (*qvMgr), false);
    return newTree;
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

void QuantifierTreeFormula::localise(const VariableSet&) {
    removeUnusedVars();
}

QuantifierTreeFormula* QuantifierTreeFormula::changeToFormula(Cudd &) {
    return this;
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

QuantifierTreeNode* QuantifierTreeFormula::getCopy() {
    if (!(getUnivVars().empty() && getExistVars().empty())) {
        throw dqbddException("Trying to copy quantifier tree formula with some variables in prefix");
    }

    QuantifierTreeFormula *newFormula = new QuantifierTreeFormula(getMgr(), (*qvMgr));
    newFormula->setMatrix(getMatrix());
    return newFormula;
}

} // namespace dqbdd
