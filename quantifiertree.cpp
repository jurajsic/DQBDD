#include <algorithm>
#include "quantifiertree.hpp"

void QuantifierTree::renameVar(Variable oldVar, Variable newVar) {
    for (QuantifierTreeNode *child : children) {
        child->renameVar(oldVar, newVar);
    }
}

void QuantifierTree::pushExistVar(Variable var) {
    addExistVar(var);
}

// TODO no need for renaming new univ vars???
void QuantifierTree::pushUnivVar(Variable var) {
    Variable newUniVar = var.newVarAtSameLevel();
    addUnivVar(newUniVar);
    VariableSet supportSet = getSupportSet();
    for (Variable eVar : getUnivVarDependencies(var)) {
        if (supportSet.contains(eVar)) {
            removeDependency(eVar, var);
            addDependency(eVar,newUniVar);
        }
    }
    renameVar(var, newUniVar);
}

// TODO implement ME!!!!
// removes from one list another one where both are in the same order
// returns true if something was removed
bool removeFromOrderedListOtherOrderedList(std::list<QuantifierTreeNode*> &listToRemoveFrom, std::list<QuantifierTreeNode*> &listOfItemsToRemove) {
    return fuck;
}

void QuantifierTree::localise() {
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

/*
        std::list<Variable> eVarsToBePushed;
        // the existential variables that are in all children stays here, we do not create new children for them
        for (const Variable &eVar : getExistVars()) {
            if (childrenContainingExistVar[eVar].size() == children.size()) {
                childrenContainingExistVar.erase(eVar);
            } else {
                eVarsToBePushed.push_back(eVar);
            }
        }
*/
        
        while (!getExistVars().empty()) {
            // find exist variable that has minimal number of occurences in children
            auto eVarToPushIter = std::min_element(getExistVars().begin(), getExistVars().end(), 
                        [&](Variable v1, Variable v2) { 
                            return (childrenContainingExistVar[v1].size() < childrenContainingExistVar[v2].size());
                        }
                     );

            // if we hit the exist variable that is in all children (meaning 
            // that all other exist variables are in all children too) we stop
            if (childrenContainingExistVar[*eVarToPushIter].size() == children.size()) {
                break;
            }

            // if *eVarToPushIter is needed to be pushed to only one child, there is no need to create new children
            if (childrenContainingExistVar[*eVarToPushIter].size() == 1) {
                (*childrenContainingExistVar[*eVarToPushIter].begin())->pushExistVar(*eVarToPushIter);
                removeExistVar(*eVarToPushIter);
                continue;
            }

            std::list<QuantifierTreeNode*> &childrenToCombine = childrenContainingExistVar[*eVarToPushIter];
            
            // create a new child to which this existential variable is pushed
            QuantifierTree *newChild = new QuantifierTree(true, childrenToCombine);
            newChild->pushExistVar(*eVarToPushIter);
            removeExistVar(*eVarToPushIter);

            // erase the children that were connected to create a new child
            removeFromOrderedListOtherOrderedList(children, childrenToCombine);
            children.push_back(newChild);
            for (const Variable &eVarToBePushed : getExistVars()) {
                removeFromOrderedListOtherOrderedList(childrenContainingExistVar[eVarToBePushed], childrenToCombine);
                childrenContainingExistVar[eVarToBePushed].push_back(newChild);
            }
        }
        

        const VariableSet &univVars = getUnivVars();
        VariableSet dependentUnivVars; // universal variables on which some exist var depends
        for (const Variable &existVar : getExistVars()) {
            dependentUnivVars.insert(getExistVarDependencies(existVar).begin(), getExistVarDependencies(existVar).end());
        }
        for (const Variable &univVarToPush : univVars) {
            if (dependentUnivVars.contains(univVarToPush)) { // no leftover existential variable depends on univVarToPush --> can be pushed
                for (QuantifierTreeNode *child : children) {
                    // TODO for push universal variables the function pushVar should rename it
                    child->pushUnivVar(univVarToPush);
                }
                removeUnivVar(univVarToPush);
            }
        }
    } else {
        for (const Variable &existVarToPush : getExistVars()) {
            for (QuantifierTreeNode *child : children) {
                // there is no need to rename them, because universal vars on which this exist var depends will either stay as ancestor/same level
                // or will be pushed into some sibling and then we can say that this exist var is not there????
                child->pushExistVar(existVarToPush);
            }
            // TODO cannot remove because it will also be removed from dependencyMap???
            removeExistVar(existVarToPush);
        }
        // for each universal variable find set of children that contain it + all exist variables which depend on it and merge them into new QuantifierTree with same operator and remove it if it was pushed

    }
}