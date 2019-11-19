#ifndef QUANTIFIERTREE_HPP
#define QUANTIFIERTREE_HPP

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "variable.hpp"
#include "formula.hpp"

//typedef std::unordered_set<Variable> VariableSet;
//typedef std::unordered_map<Variable, VariableSet> DependencyMap;

class QuantifierTreeNode {
public:    
    virtual VariableSet getSupportSet() = 0;
    virtual void localise() = 0; // TODO ake parametre a return???
    virtual void pushExistVar(Variable var) = 0;
    virtual void pushUnivVar(Variable var) = 0;
};


// TODO zrobit bud ako samostatnu class kde si budem ukladat dependencies a nechat reference na to nejak v strome
class QuantifierTreeRoot : public QuantifierTreeNode {
private:
};

class QuantifierTree : public QuantifierTreeNode, public QuantifiedVariablesManipulator {
private:
    std::vector<QuantifierTreeNode*> children;
    bool isConj; // TODO change to operator, will have to learn how to do that tho
public:
    void localise();
    void pushExistVar(Variable var);
    void pushUnivVar(Variable var);
};

#endif