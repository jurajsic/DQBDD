#ifndef DQBF_BDD_HQSPRE_INTERFACE_HPP
#define DQBF_BDD_HQSPRE_INTERFACE_HPP

#include <memory>
#include "parser.hpp"

class HQSPreInterface : public Parser {
private:
    // using pimpl idiom to hide implementation of HQSPre
    class HQSPreFormulaWrapper;
    std::unique_ptr<HQSPreFormulaWrapper> formulaPtr;

    Cudd &mgr;
    QuantifiedVariablesManipulator DQBFPrefix;
public:
    HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses a file in DQDIMACS format and runs HQSpre preprocessor
     * 
     * @param fileName name of the file to parse
     */
    void parse(std::string fileName);
    Formula* getFormula();
    QuantifierTreeNode* getQuantifierTree();
};

#endif