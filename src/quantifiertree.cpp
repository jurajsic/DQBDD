#include <algorithm>
#include "quantifiertree.hpp"


void QuantifierTreeNode::pushExistVar(Variable var) {
    if (getSupportSet().contains(var)) {
        addExistVar(var);
    }
}

// TODO no need for renaming new univ vars???
void QuantifierTreeNode::pushUnivVar(Variable var) {
    // no renaming??? should work
    VariableSet supportSet = getSupportSet();
    if (supportSet.contains(var)) {
        addUnivVar(var);
    } else { // we delete dependencies if var is not in support set??? hopefully no problem
        for (const Variable &dependentExistVar : getUnivVarDependencies(var)) {
            if (supportSet.contains(dependentExistVar)) { // we assume this is the only tree that contains dependentExistVar, shoudl work
                removeDependency(dependentExistVar,var);
            }
        }
    }
}

QuantifierTree::QuantifierTree(bool isConj, std::list<QuantifierTreeNode*> children, QuantifiedVariablesManager &qvMgr) : QuantifiedVariablesManipulator(qvMgr), isConj(isConj) {
    supportSet = {};
    for (QuantifierTreeNode *child : children) {
        addChild(child);
    }
}

QuantifierTree::~QuantifierTree() {
    for (QuantifierTreeNode *child : children) {
        delete child;
    }
}

// TODO implement ME!!!!
// removes from one list another one where both are in the same order (the order is based on the order of children)
// returns true if something was removed
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
        // TODO for each exist var find set of children that contain it 
        // and merge them into new QuantifierTree with same operator
        std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingExistVar;
        for (QuantifierTreeNode *child : children) {
            VariableSet supportSet = child->getSupportSet();
            for (const Variable &eVar : getExistVars()) {
                if (supportSet.contains(eVar)) {
                    childrenContainingExistVar[eVar].push_back(child);
                }
            }
        }
        
        // push exist vars while it is possible
        while (!getExistVars().empty()) {
            // find exist variable that has minimal number of occurences in children
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
        
        // find all univ vars that can be pushed (they do not depend on any exist var in current node)
        VariableSet dependentUnivVars; // universal variables on which some exist var depends
        for (const Variable &existVar : getExistVars()) {
            dependentUnivVars.insert(getExistVarDependencies(existVar).begin(), getExistVarDependencies(existVar).end());
        }

        // TODO delete dependencies for children that do not contain univ var
        // push possible univ vars
        VariableSet univVars = getUnivVars();
        for (const Variable &univVarToPush : univVars) {
            if (!dependentUnivVars.contains(univVarToPush)) { // no leftover existential variable depends on univVarToPush --> can be pushed
                for (QuantifierTreeNode *child : children) {
                    // TODO for push universal variables the function pushVar should rename it
                    child->pushUnivVar(univVarToPush);
                }
                removeUnivVar(univVarToPush);
            }
        }
    } else {
        VariableSet existVars = getExistVars();
        for (const Variable &existVarToPush : existVars) {
            for (QuantifierTreeNode *child : children) {
                // there is no need to rename them, because universal vars on which this exist var depends will either stay as ancestor/same level
                // or will be pushed into some sibling and then we can say that this exist var is not there????
                child->pushExistVar(existVarToPush);
            }
            // TODO cannot remove because it will also be removed from dependencyMap???
            removeExistVar(existVarToPush);
        }

        // for each universal variable find set of children that contain it + all exist variables which depend
        // on it and merge them into new QuantifierTree with same operator and remove it if it was pushed
        std::unordered_map<Variable,std::list<QuantifierTreeNode*>> childrenContainingUnivVar;
        for (QuantifierTreeNode *child : children) {
            VariableSet supportSet = child->getSupportSet();
            for (const Variable &uVar : getUnivVars()) {
                // if the child contains uVar...
                if (supportSet.contains(uVar)) {
                    childrenContainingUnivVar[uVar].push_back(child);
                    continue;
                }
                // ... or some exist var that depends on uVar
                for (const Variable &eVar : getUnivVarDependencies(uVar)) {
                    if (supportSet.contains(eVar)) {
                        childrenContainingUnivVar[uVar].push_back(child);
                        continue;
                    }
                }
            }
        }
        
        while (!getUnivVars().empty()) {
            // find univ variable that has minimal number of occurences in children
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

            // if *uVarToPushIter is needed to be pushed to only one child, there is no need to create new children
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

    for (QuantifierTreeNode *child : children) {
        child->localise();
    }
}

QuantifierTreeFormula* QuantifierTree::changeToFormula(Cudd &mgr) {
    auto childIter = children.begin();
    BDD matrix;
    if (isConj) {
        matrix = mgr.bddOne();
    } else {
        matrix = mgr.bddZero();
    }
    while (childIter != children.end()) {
        QuantifierTreeFormula *childFormula = (*childIter)->changeToFormula(mgr);

        // remove the child 
        auto childToRemoveIter = childIter;
        ++childIter;
        // TODO maybe change to shared_ptr so I don't have to this shit
        //delete (*childToRemoveIter);
        children.erase(childToRemoveIter);

        // TODO implement this
        childFormula->eliminatePossibleVars();
        
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

        // move all exist vars up (univ vars are not moved, because they were eliminated)
        for (const Variable &eVar : childFormula->getExistVars()) {
            addExistVar(eVar);
            // no need for this because we delete childFormula after and it is done automatically (and it wouldnt work because childFormula->getExistVars() returns reference)
            //childFormula->removeExistVar(eVar);
        }

        delete childFormula;
    }
    
    // delete leftover children
    while (childIter != children.end()) {
        // remove the child 
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
    children.push_back(child);
    for (const Variable var : child->getSupportSet()) {
        supportSet.insert(var);
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
                            : QuantifiedVariablesManipulator(qvmgr), Formula(mgr, qvmgr) {}

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator)
                            : QuantifiedVariablesManipulator(qvManipulator), Formula(mgr, qvManipulator) {}

void QuantifierTreeFormula::localise() {
    removeUnusedVars();
}

QuantifierTreeFormula* QuantifierTreeFormula::changeToFormula(Cudd &mgr) {
    return this;
}

void QuantifierTreeFormula::negate() {
    setMatrix(!getMatrix());
}