#ifndef DQBDD_PARSER_HPP
#define DQBDD_PARSER_HPP

#include "DQBDDformula.hpp"
#include "quantifiertree.hpp"

/**
 * @brief Base function for formula parsers to inherit from.
 */
class Parser {
public:
    Parser() = default;
    virtual ~Parser() = default;
    // returns true if resulting formula is trivial (equal to TRUE or FALSE) - can also mean that preprocessor solved
    virtual bool parse(std::string fileName) = 0;
    virtual Formula* getFormula() = 0;
    virtual QuantifierTreeNode* getQuantifierTree() = 0;
};

#endif