#ifndef DQBF_BDD_PARSER_HPP
#define DQBF_BDD_PARSER_HPP

#include "formula.hpp"
#include "quantifiertree.hpp"

/*
class ParserTree {
public:
    enum OperationType {
        AND,
        OR,
        XOR,
        EQIV,
        VAR,
        TRUE,
        FALSE
    };
private:
    OperationType opType;
    bool isNegated;

    // based on operation type one of these is used
    std::list<ParserTree*> parserTreeChildren;
    Variable *variableChild = nullptr;
public:
    ParserTree() = delete;*/
    /**
     * @brief Constructor for operators TRUE and FALSE
     * 
     * @param opType either TRUE or FALSE
     */
    //ParserTree(OperationType opType);
    /**
     * @brief Constructor for operators with operands
     * 
     * @param opType 
     * @param parserTreeChildren operands
     * @param isNegated is the result negated
     */
    //ParserTree(OperationType opType, const std::list<ParserTree*> &parserTreeChildren, bool isNegated = false);
    /**
     * @brief Constructor for variable type
     * 
     * @param var 
     * @param isNegated is the variable negated
     */
 /*   ParserTree(Variable &var, bool isNegated = false);
    ~ParserTree();

    void negate();
*/
    //void addParserTree(ParserTree *parserTreeToAdd);
    //void removeParserTree(ParserTree *parserTreeToRemove);

/*
    void addVariable(Variable varToAdd);
    void removeVariable(Variable varToRemove);
*/
    //Formula* transformToFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    //QuantifierTree* transformToQuantifierTree(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);

//};


class Parser {
protected:
    //ParserTree *parserTree = nullptr;
    // TODO initalize
public:
    Parser() = default;
    virtual ~Parser() = default;
    virtual void parse(std::string fileName) = 0;
    virtual Formula* getFormula() = 0;
    virtual QuantifierTreeNode* getQuantifierTree() = 0;
    //ParserTree* getParserTree();
};

#endif