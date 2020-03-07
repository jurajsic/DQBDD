/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DQBDD_VARIABLE_HPP
#define DQBDD_VARIABLE_HPP

#include <unordered_set>

#include "cuddObj.hh"

class Variable {
private:
    // each variable has unique ID based on Cudd manager
    unsigned int id;
    BDD representation;
    Cudd mgr;

public:
    operator BDD() const;
    Variable() = delete;
    Variable(int id, Cudd &mgr);
    //Variable(BDD repr);
    unsigned int getId() const;
    BDD getBDD() const;
    // get a new Variable that is at the same BDD level as this variable
    Variable newVarAtSameLevel();
    int getLevel() const;
    bool operator==(const Variable &anotherVariable) const;
    //BDD operator&(const Variable& other) const;
    //BDD operator|(const Variable& other) const;
    BDD operator&(const BDD& other) const;
    BDD operator|(const BDD& other) const;
    BDD operator!() const;
};

std::ostream& operator<<(std::ostream& os, const Variable& obj);

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

#endif