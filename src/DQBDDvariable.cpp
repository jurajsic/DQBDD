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

#include <functional>
#include <iostream>

#include "DQBDDvariable.hpp"

Variable::Variable(int id, Cudd &mgr) : id(id), mgr(mgr) {
    representation = mgr.bddVar(id);
}

/*
Variable::Variable(BDD repr) : representation(repr) {
    id = repr.NodeReadIndex();
}
*/

int Variable::getLevel() const {
    return mgr.ReadPerm(id);
}

Variable Variable::newVarAtSameLevel() {
    return Variable(mgr.bddNewVarAtLevel(getLevel()).NodeReadIndex(), mgr);
}

Variable::operator BDD() const {
    return getBDD();
}

unsigned int Variable::getId() const {
    return id;
}


BDD Variable::getBDD() const {
    return representation;
}


bool Variable::operator==(const Variable &anotherVariable) const {
    return (id == anotherVariable.id);
}

BDD Variable::operator&(const BDD& other) const {
    return getBDD() & other;
}

BDD Variable::operator|(const BDD& other) const {
    return getBDD() | other;
}

BDD Variable::operator!() const {
    return !getBDD();
}

std::ostream& operator<<(std::ostream& os, const Variable& obj)
{
    os << obj.getBDD();
    return os;
}