#include <functional>
#include "Variable.hpp"

Variable::Variable(int id) : id(id) {

}

int Variable::getId() const {
    return id;
}

bool Variable::operator==(const Variable &anotherVariable) const {
    return (id == anotherVariable.id);
}

namespace std
{
    template <>
    struct hash<Variable>
    {
        size_t operator()(const Variable& k) const
        {
            return hash<int>()(k.getId());
        }
    };
}
/*
ExistVariable::ExistVariable(int id) : Variable(id) {

}

UnivVariable::UnivVariable(int id) : Variable(id) {

}*/