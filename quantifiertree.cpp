#include "quantifiertree.hpp"

void QuantifierTree::renameVar(Variable oldVar, Variable newVar) {
    for (QuantifierTreeNode *child : children) {
        child->renameVar(oldVar, newVar);
    }
}

void QuantifierTree::pushExistVar(Variable var) {
    addExistVar(var);
}

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

void QuantifierTree::localise() {
    if (isConj) {
        VariableSet existVars = getExistVars();
        // TODO for each exist var find set of children that contain it and merge them into new QuantifierTree with same operator
        VariableSet univVars = getUnivVars();
        VariableSet dependentUnivVars; // universal variables on which some exist var depends
        for (const Variable &existVar : getExistVars()) {
            // TODO move this bit into previous exist vars loop when exist var is not pushed
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
        // for each universal variable find set of children that contain it + all exist variables which depend on it and merge them into new QuantifierTree with same operator and remove it if it was pushed
    }
}