#ifndef DQBDD_EXCEPTION_HPP
#define DQBDD_EXCEPTION_HPP

#include <stdexcept>

class DQBDDexception : public std::logic_error
{
    using std::logic_error::logic_error;
};

#endif