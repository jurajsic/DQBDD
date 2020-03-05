#ifndef DQBF_BDD_DQBDD_EXCEPTION_HPP
#define DQBF_BDD_DQBDD_EXCEPTION_HPP

#include <stdexcept>

class DQBDDexception : public std::logic_error
{
    using std::logic_error::logic_error;
};

#endif