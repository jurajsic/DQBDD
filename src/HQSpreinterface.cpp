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

#include <easylogging++.hpp>
#include <formula.hpp>

#include "HQSpreinterface.hpp"
#include "DQBDDexceptions.hpp"

INITIALIZE_EASYLOGGINGPP

class HQSPreInterface::HQSPreFormulaWrapper {
public:
    hqspre::Formula formula;
    bool isSat = false;
    bool isUnSat = false;
    HQSPreFormulaWrapper() {}
    ~HQSPreFormulaWrapper() {}

    // following stuff is based on similar functions for AIGs from HQS

    /**
     * @brief gate table mapping HQSpre vars to BDDs (either variables or BDDs representing extracted gates)
     */
    std::vector<BDD> gate_table;

    BDD convertLiteral(const hqspre::Literal lit)
    {
        const hqspre::Variable var = hqspre::lit2var(lit);
        BDD res = gate_table[var];
        return (hqspre::isNegative(lit) ? !res : res);
    }

    BDD gateToBDD(const Cudd &mgr, const hqspre::Gate &g) {
        BDD result;

        switch (g._type) {
            case hqspre::GateType::AND_GATE:
            {
                result = mgr.bddOne();
                
                for (const hqspre::Literal lit : g._input_literals) {
                    result &= convertLiteral(lit);
                }
                break;
            }

            case hqspre::GateType::XOR_GATE:
            {
                result = convertLiteral(g._input_literals[0]).Xor(convertLiteral(g._input_literals[1]));
                break;
            }

            case hqspre::GateType::MUX_GATE:
            {
                // MUX(A,B,C) = (A AND B) OR (!A AND C)
                BDD A = convertLiteral(g._input_literals[0]);
                BDD B = convertLiteral(g._input_literals[1]);
                BDD C = convertLiteral(g._input_literals[2]);
                result = (A & B) | ((~A) & C);
                break;
            }

            default:
            {
                throw DQBDDexception("Invalid gate type encountered");
                break;
            }
        }

        return (hqspre::isNegative(g._output_literal) ? !result : result);
    }

    // stuff for transforming to 

    std::vector<const hqspre::Gate*> outputvarToGate;

    QuantifierTreeNode* literalToTree(Cudd &mgr, const hqspre::Literal lit, QuantifiedVariablesManager &qvm) {
        const hqspre::Variable var = hqspre::lit2var(lit);
        //std::cout << "Processing variable " << var << std::endl;

        QuantifierTreeNode *result;

        if (formula.isGateOutput(var)) {
            result = gateToTree(mgr, *outputvarToGate[var], qvm);
        } else {
            //std::cout << "not a gate" << std::endl;
            QuantifierTreeFormula *varFormula = new QuantifierTreeFormula(mgr, qvm);
            varFormula->setMatrix(Variable(var, mgr));
            result = varFormula;
        }

        if (hqspre::isNegative(lit)) {
            result->negate();
        }
        return result;
    }

    QuantifierTreeNode* gateToTree(Cudd &mgr, const hqspre::Gate &g, QuantifiedVariablesManager &qvm) {
        QuantifierTreeNode* result;
        switch (g._type) {
            case hqspre::GateType::AND_GATE:
            {
                //std::cout << "AND GATE" << std::endl;
                std::list<QuantifierTreeNode*> operands;
                for (const hqspre::Literal lit : g._input_literals) {
                    operands.push_back(literalToTree(mgr, lit, qvm));
                }

                if (operands.size() == 1) {
                    result = *operands.begin();
                } else {
                    result = new QuantifierTree(true, operands, qvm);
                }

                break;
            }

            case hqspre::GateType::XOR_GATE:
            {
                //std::cout << "XOR GATE" << std::endl;
                // A XOR B = (A AND !B) OR (!A AND B)
                QuantifierTreeNode *firstOp = literalToTree(mgr, g._input_literals[0], qvm);
                QuantifierTreeNode *firstOpNeg = literalToTree(mgr, g._input_literals[0], qvm);
                firstOpNeg->negate();
                QuantifierTreeNode *secondOp = literalToTree(mgr, g._input_literals[1], qvm);
                QuantifierTreeNode *secondOpNeg = literalToTree(mgr, g._input_literals[1], qvm);
                secondOpNeg->negate();
                QuantifierTreeNode *firstConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOp, secondOpNeg}, qvm);
                QuantifierTreeNode *secondConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOpNeg, secondOp}, qvm);
                result = new QuantifierTree(false, std::list<QuantifierTreeNode*>{firstConj, secondConj}, qvm);
                break;
            }

            case hqspre::GateType::MUX_GATE:
            {
                //std::cout << "MUX GATE" << std::endl;
                // MUX(A,B,C) = (A AND B) OR (!A AND C)
                QuantifierTreeNode *firstOp = literalToTree(mgr, g._input_literals[0], qvm);
                QuantifierTreeNode *firstOpNeg = literalToTree(mgr, g._input_literals[0], qvm);
                firstOpNeg->negate();
                QuantifierTreeNode *secondOp = literalToTree(mgr, g._input_literals[1], qvm);
                QuantifierTreeNode *thirdOp = literalToTree(mgr, g._input_literals[2], qvm);
                QuantifierTreeNode *firstConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOp, secondOp}, qvm);
                QuantifierTreeNode *secondConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOpNeg, thirdOp}, qvm);
                result = new QuantifierTree(false, std::list<QuantifierTreeNode*>{firstConj, secondConj}, qvm);
                break;
            }

            default:
            {
                throw DQBDDexception("Invalid gate type encountered");
                break;
            }
        }

        if (hqspre::isNegative(g._output_literal)) {
            result->negate();
        }
        return result;
    }
};

HQSPreInterface::HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : formulaPtr(nullptr), mgr(mgr), DQBFPrefix(qvmgr) {}

HQSPreInterface::~HQSPreInterface() = default;

// code based on hqsfork
bool HQSPreInterface::parse(std::string fileName) {
    formulaPtr.reset(new HQSPreFormulaWrapper());

    // Configure logging
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Enabled, "true");
    defaultConf.setGlobally(el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.set(el::Level::Verbose, el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
    defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
    el::Loggers::reconfigureAllLoggers(defaultConf);
    el::Loggers::setVerboseLevel(1);

    formulaPtr->formula.settings().consistency_check = false;

    // Parse the file
    std::string in_name(fileName);
    if (in_name == "") {
        throw DQBDDexception("No input file given.");
    }
    std::ifstream in(in_name);
    if (!in) {
        throw DQBDDexception("Could not open input file '");
    }
    in >> formulaPtr->formula;
    in.close();

    // do the preprocessing magic
    try {
        formulaPtr->formula.determineGates(true, true, true, false);
        if (formulaPtr->formula.getGates().size() > 5) {
            // First do full preprocessing on a copy of the formula
            hqspre::Formula formula2(formulaPtr->formula);
            //formula2.settings().bla              = false;
            //formula2.settings().ble              = false;
            //formula2.settings().pure_sat_timeout = 1000;
            formula2.preprocess();

            // Then do preprocessing, preserving gates
            formulaPtr->formula.settings().univ_expand      = 0; // maybe 2 is better???
            formulaPtr->formula.settings().bla              = false;
            formulaPtr->formula.settings().ble              = false;
            formulaPtr->formula.settings().preserve_gates   = true;
            formulaPtr->formula.settings().substitution     = false;
            formulaPtr->formula.settings().rewrite          = false;
            formulaPtr->formula.settings().resolution       = false;
            formulaPtr->formula.settings().max_loops        = 20;
            //formulaPtr->formula.settings().pure_sat_timeout = 1000;
        }
        formulaPtr->formula.preprocess();
        formulaPtr->formula.printStatistics();
    } catch (hqspre::SATException&) {
        formulaPtr->isSat = true;
        return true;
    } catch (hqspre::UNSATException&) {
        formulaPtr->isUnSat = true;
        return true;
    }

    // do gate extraction
    formulaPtr->formula.determineGates(true, true, true, false);
    formulaPtr->formula.enforceDQBF(true);
    formulaPtr->formula.unitPropagation();

    const auto gates = formulaPtr->formula.getGates();

    formulaPtr->gate_table = std::vector<BDD>(formulaPtr->formula.maxVarIndex() + 1);
    formulaPtr->outputvarToGate = std::vector<const hqspre::Gate*>(formulaPtr->formula.maxVarIndex() + 1);

    // Create the proper problem variables (without Tseitin variables)
    for (hqspre::Variable var = formulaPtr->formula.minVarIndex(); var <= formulaPtr->formula.maxVarIndex(); ++var) {
        if (formulaPtr->formula.isUniversal(var)) {
            Variable uVar = Variable(var, mgr);
            DQBFPrefix.addUnivVar(uVar);
            formulaPtr->gate_table[var] = uVar;
        } else if (formulaPtr->formula.isExistential(var) && !formulaPtr->formula.isGateOutput(var)) {
            Variable eVar = Variable(var, mgr);
            DQBFPrefix.addExistVar(eVar);
            for (auto dep : formulaPtr->formula.getDependencies(var)) {
                DQBFPrefix.addDependency(eVar, Variable(dep, mgr));
            }
            formulaPtr->gate_table[var] = eVar;
        } 
    }

    for (auto &gate : formulaPtr->formula.getGates()) {
        formulaPtr->outputvarToGate[hqspre::lit2var(gate._output_literal)] = &gate;
    }

    return false;
}

Formula* HQSPreInterface::getFormula() {
    if (formulaPtr == nullptr) {
        throw DQBDDexception("A file must be parsed before it is possible to get formula");
    }

    Formula *DQBFformula = new Formula(mgr, DQBFPrefix);
    DQBFPrefix.clear();

    if (formulaPtr->isSat) {
        DQBFformula->setMatrix(mgr.bddOne());
        return DQBFformula;
    }
    if (formulaPtr->isUnSat) {
        DQBFformula->setMatrix(mgr.bddZero());
        return DQBFformula;
    }

    // Create the gates
    for (const hqspre::Gate& g: formulaPtr->formula.getGates()) {
        const hqspre::Variable out_var = hqspre::lit2var(g._output_literal);
        formulaPtr->gate_table[out_var] = formulaPtr->gateToBDD(mgr, g);
        for (const auto c_nr: g._encoding_clauses) formulaPtr->formula.removeClause(c_nr);
    }

    // Convert the clauses
    BDD matrix = mgr.bddOne();
    for (hqspre::ClauseID c_nr = 0; c_nr <= formulaPtr->formula.maxClauseIndex(); c_nr++) {
        if (!formulaPtr->formula.clauseDeleted(c_nr)) {
            BDD clauseBDD = mgr.bddZero();
            const auto& clause = formulaPtr->formula.getClause(c_nr);
            for (const hqspre::Literal lit: clause) {
                clauseBDD |= formulaPtr->convertLiteral(lit);
            }
            matrix &= clauseBDD;
        }
    }

    DQBFformula->setMatrix(matrix);
    return DQBFformula;
}

QuantifierTreeNode* HQSPreInterface::getQuantifierTree() {
    if (formulaPtr == nullptr) {
        throw DQBDDexception("A file must be parsed first before it is possible to get quantifier tree");
    }

    if (formulaPtr->isSat || formulaPtr->isUnSat) {
        QuantifierTreeFormula *DQBFformula = new QuantifierTreeFormula(mgr, DQBFPrefix);
        DQBFPrefix.clear();
        if (formulaPtr->isSat) {
            DQBFformula->setMatrix(mgr.bddOne());
            return DQBFformula;
        }
        if (formulaPtr->isUnSat) {
            DQBFformula->setMatrix(mgr.bddZero());
            return DQBFformula;
        }
    }

    // delete clauses from which gates were generated
    for (const hqspre::Gate& g: formulaPtr->formula.getGates()) {
        for (const auto c_nr: g._encoding_clauses) formulaPtr->formula.removeClause(c_nr);
    }

    // find only clauses which were not deleted during preprocessing
    std::vector<hqspre::ClauseID> workingClausesNrs = {};
    for (hqspre::ClauseID c_nr = 0; c_nr <= formulaPtr->formula.maxClauseIndex(); c_nr++) {
        if (!formulaPtr->formula.clauseDeleted(c_nr)) {
            //std::cout << "Clause " << c_nr << " will be processed" << std::endl;
            workingClausesNrs.push_back(c_nr);
        }
    }

    // if we have only one clause, we make it a root...
    if (workingClausesNrs.size() == 1) {
        const auto& clause = formulaPtr->formula.getClause(workingClausesNrs[0]);
        std::list<QuantifierTreeNode*> literals;
        for (const hqspre::Literal lit: clause) {
            literals.push_back(formulaPtr->literalToTree(mgr, lit, *DQBFPrefix.getManager()));
        }
        QuantifierTree *clauseTree = new QuantifierTree(false, literals, DQBFPrefix);
        DQBFPrefix.clear();
        return clauseTree;
    }

    // ...otherwise we save all the clauses...
    std::list<QuantifierTreeNode*> clauses;
    for (hqspre::ClauseID c_nr : workingClausesNrs) {
        //std::cout << "Processing clause " << c_nr << " with literals:" << std::endl;
        const auto& clause = formulaPtr->formula.getClause(c_nr);
        std::list<QuantifierTreeNode*> literals;
        //for (const hqspre::Literal lit: clause) {
        //    std::cout << hqspre::lit2var(lit) << " ";
        //}
        //std::cout << std::endl;
        for (const hqspre::Literal lit: clause) {
            literals.push_back(formulaPtr->literalToTree(mgr, lit, *DQBFPrefix.getManager()));
        }
        QuantifierTree *clauseTree = new QuantifierTree(false, literals, *DQBFPrefix.getManager());
        clauses.push_back(clauseTree);
    }

    // ...and then make a root from their conjuction
    auto qt = new QuantifierTree(true, clauses, DQBFPrefix);
    DQBFPrefix.clear();
    return qt;
}