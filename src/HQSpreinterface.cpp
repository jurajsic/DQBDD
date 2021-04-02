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

#include<vector>

#include <easylogging++.hpp>
#include <formula.hpp>

#include "HQSpreinterface.hpp"
#include "DQBDDexceptions.hpp"

INITIALIZE_EASYLOGGINGPP

/**
 * Everything here is based on the way how preprocessing is done in HQS, especially
 * how the result of the preprocessing is tranformed into AIG.
 **/

class HQSPreInterface::HQSPreFormulaWrapper {
private:
    // this will map hqspre vars to dqbdd vars (we use it to create dqbdd vars starting from 1,2,...)
    std::unordered_map<hqspre::Variable,Variable> hqspreVarToDqbddVar;
public:
    hqspre::Formula formula;
    bool isSat = false;
    bool isUnSat = false;
    HQSPreFormulaWrapper() {}
    ~HQSPreFormulaWrapper() {}

    /* this function basically does hqspreVarToDqbddVar[hqspreVar] = dqbddVar but it assumes
     * that hqspreVarToDqbddVar[hqspreVar] does not exists, otherwise it throws exception
     */
    void mapHqspreVarToDqbddVar(const hqspre::Variable &hqspreVar, const Variable &dqbddVar) {
        if (!hqspreVarToDqbddVar.emplace(hqspreVar, dqbddVar).second) { // if hqspreVarToDqbddVar[hqspreVar] already exists 
            throw DQBDDexception("We cannot map two different DQBDD variables into one hqspre variable");
        }
    }

    /**
     * @brief Returns hqspreVarToDqbddVar[hqspreVar] if it exists, otherwise throws error
     */
    const Variable &getDqbddVarMappedIntoHqspreVar(const hqspre::Variable &hqspreVar) {
        auto foundIt = hqspreVarToDqbddVar.find(hqspreVar);
        if (foundIt != hqspreVarToDqbddVar.end()) {
            return foundIt->second;
        } else {
            throw DQBDDexception("HQSpre variable does not have any DQBDD variable mapped to it");
        }
    }

    /*****************************************************************************************/
    /*************** Stuff for transforming hqspre formula into BDD formula ******************/
    /*****************************************************************************************/

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

    
    /*****************************************************************************************/
    /************* Stuff for transforming hqspre formula into quantifier tree ****************/
    /*****************************************************************************************/

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
            varFormula->setMatrix(getDqbddVarMappedIntoHqspreVar(var));
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

// code based on HQS 
bool HQSPreInterface::parse(std::string fileName) {
    // if this interface parsed something before, we need to remove saved information
    formulaPtr.reset(new HQSPreFormulaWrapper());
    DQBFPrefix.clear();

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

    try {
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
        formulaPtr->formula.determineGates();
        auto numOfDetectedGates = formulaPtr->formula.getGates().size();
        VLOG(1) << "Detected " << numOfDetectedGates << " gates in the input formula";
        if (numOfDetectedGates > 5) { 
            /* If we have more than 5 gates, we would like to preserve them, but first we would like to
             * try HQSpre to solve the formula without preserving them
             */

            // First do full preprocessing (without preserving gates) on a copy of the formula
            hqspre::Formula formula2(formulaPtr->formula);
            // TODO decide what settings to use (here are except timeout same as in HQS 18-03-2021)
            // bla and ble are only useful for QBF (I think)
            formula2.settings().bla              = false;
            formula2.settings().ble              = false;
            formula2.settings().impl_chains      = 3;
            formula2.settings().max_substitution_cost = 250;
            formula2.settings().max_resolution_cost = 100;
            formula2.settings().vivify_fp        = true;
            // TODO timeout?
            //formula2.settings().pure_sat_timeout = 1000;
            // in HQS forksplitting is enabled, but it uses external MPhaseSAT64 solver, so here we better turn it off
            formula2.settings().enableFork       = false;
            VLOG(1) << "Start preprocessing without gate preservation";
            formula2.preprocess();

            // Then do preprocessing, preserving gates
            // TODO setting??? (here are except timeout same as in HQS 18-03-2021)
            formulaPtr->formula.settings().univ_expand      = 0; // maybe 2 is better???
            // bla and ble are only useful for QBF (I think)
            formulaPtr->formula.settings().bla              = false;
            formulaPtr->formula.settings().ble              = false;
            formulaPtr->formula.settings().preserve_gates   = true;
            formulaPtr->formula.settings().substitution     = false;
            formulaPtr->formula.settings().rewrite          = false;
            formulaPtr->formula.settings().resolution       = false;
            formulaPtr->formula.settings().max_loops        = 20;
            formulaPtr->formula.settings().enableFork       = false;
            // TODO timeout??
            //formulaPtr->formula.settings().pure_sat_timeout = 1000;

            VLOG(1) << "Start preprocessing with gate preservation";
        } else { // if we have no gates, there is no points in preserving them, run HQSpre normally
            // TODO decide what settings to use (here are except timeout same as in HQS 18-03-2021)
            // bla and ble are only useful for QBF (I think)
            formulaPtr->formula.settings().bla              = false;
            formulaPtr->formula.settings().ble              = false;
            formulaPtr->formula.settings().impl_chains      = 3;
            formulaPtr->formula.settings().max_substitution_cost = 250;
            formulaPtr->formula.settings().max_resolution_cost = 100;
            formulaPtr->formula.settings().vivify_fp        = true;
            // TODO timeout?
            //formulaPtr->formula.settings().pure_sat_timeout = 1000;
            // in HQS forking is enabled, but it uses external MPhaseSAT64 solver, so here we better turn it off
            formulaPtr->formula.settings().enableFork       = false;

            VLOG(1) << "Start preprocessing without gate preservation";
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

    // update gates (if they were not preserved they might have gotten deleted)
    formulaPtr->formula.determineGates();
        
    // enforce DQBF (to have DQBF prefix, it could have possibly changed to QBF)
    formulaPtr->formula.enforceDQBF(true);
    

    // Create the proper problem variables (without Tseitin variables)
    // we will use varNum as a way to have variables starting from one
    int varNum = 1;
    // first universal...
    for (hqspre::Variable var = formulaPtr->formula.minVarIndex(); var <= formulaPtr->formula.maxVarIndex(); ++var) {
        if (!formulaPtr->formula.varDeleted(var) && !formulaPtr->formula.isGateOutput(var) && formulaPtr->formula.isUniversal(var)) {
            Variable uVar = Variable(varNum, mgr);
            ++varNum;
            DQBFPrefix.addUnivVar(uVar);
            formulaPtr->mapHqspreVarToDqbddVar(var, uVar);
        }
    }
    // ...and then existential (so we can add dependencies without worrying whether univ var to be added to dependency was created or not)
    for (hqspre::Variable var = formulaPtr->formula.minVarIndex(); var <= formulaPtr->formula.maxVarIndex(); ++var) {
        if (!formulaPtr->formula.varDeleted(var) && !formulaPtr->formula.isGateOutput(var) && formulaPtr->formula.isExistential(var)) {
            Variable eVar = Variable(varNum, mgr);
            ++varNum;
            DQBFPrefix.addExistVar(eVar);
            formulaPtr->mapHqspreVarToDqbddVar(var, eVar);
            for (auto dep : formulaPtr->formula.getDependencies(var)) {
                if (!formulaPtr->formula.isUniversal(dep)) {
                    throw DQBDDexception("Existential variable is depending on non universal one in hqspre, this should not happen");
                }
                DQBFPrefix.addDependency(eVar, formulaPtr->getDqbddVarMappedIntoHqspreVar(dep));
            }
        }
    }
    
    VLOG(1) << "Preprocessing finished with formula with " << DQBFPrefix.getExistVars().size() << " existential vars, " << DQBFPrefix.getUnivVars().size() << " univ. vars, "
            << formulaPtr->formula.getGates().size() << " gates, and " << formulaPtr->formula.numClauses() << " clauses which potentionally encode the gates.";

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

    // gates need to be copied, otherwise it can get a bit funky for some reason
    const auto gates = formulaPtr->formula.getGates();

    formulaPtr->gate_table = std::vector<BDD>(formulaPtr->formula.maxVarIndex() + 1);

    // Load the gate_table with the primary inputs
    for (auto var = formulaPtr->formula.minVarIndex(); var <= formulaPtr->formula.maxVarIndex(); ++var) {
        if (!formulaPtr->formula.varDeleted(var) && !formulaPtr->formula.isGateOutput(var)) {
            formulaPtr->gate_table[var] = formulaPtr->getDqbddVarMappedIntoHqspreVar(var);
        }
    }

    // Create the gates
    for (const hqspre::Gate& g: gates) {
        const hqspre::Variable out_var = hqspre::lit2var(g._output_literal);
        formulaPtr->gate_table[out_var] = formulaPtr->gateToBDD(mgr, g);
        for (const auto c_nr: g._encoding_clauses) {
            formulaPtr->formula.removeClause(c_nr);
        }
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

    // gates need to be copied, otherwise it can get a bit funky for some reason
    const auto gates = formulaPtr->formula.getGates();

    formulaPtr->outputvarToGate = std::vector<const hqspre::Gate*>(formulaPtr->formula.maxVarIndex() + 1);

    for (auto &gate : gates) {
        formulaPtr->outputvarToGate[hqspre::lit2var(gate._output_literal)] = &gate;
        // delete clauses from which gates were generated
        for (const auto c_nr: gate._encoding_clauses) {
            formulaPtr->formula.removeClause(c_nr);
        }
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

void HQSPreInterface::turnIntoDQCIR(std::ostream &output) {
    if (formulaPtr == nullptr) {
        throw DQBDDexception("A file must be parsed before it is possible to turn it into DQCIR format");
    }

    if (formulaPtr->isSat) {
        output << "#QCIR-G14 0" << std::endl
               << "output(1)" << std::endl
               << "1 = and()" << std::endl;
        return;
    }
    if (formulaPtr->isUnSat) {
        output << "#QCIR-G14 0" << std::endl
               << "output(1)" << std::endl
               << "1 = or()" << std::endl;
        return;
    }


    output << "#QCIR-G14 " 
           // TODO what number here? number of normal variables, or number of all (normal + gate) variables?? the format is unclear about it
           //<< formulaPtr->hqspreVarToDqbddVar.size() 
           //<< formulaPtr->formula.maxVarIndex() 
           << std::endl;

    /*********************************/
    /**** Print quantifier prefix ****/
    /*********************************/
    // TODO if formula is QBF, we should print QBF style prefix (just QCIR)

    //first universal variables
    if (!DQBFPrefix.getUnivVars().empty()) {
        output << "forall(";
        for (auto uVarIter = DQBFPrefix.getUnivVars().begin(); uVarIter != DQBFPrefix.getUnivVars().end(); ++uVarIter) {
            if (uVarIter == DQBFPrefix.getUnivVars().begin()) {
                output << uVarIter->getId();
            } else {
                output << std::string(", ") << uVarIter->getId();
            }
        }
        output << std::string(")") << std::endl;
    }

    //then existential variables
    for (const Variable &eVar : DQBFPrefix.getExistVars()) {
        output << std::string("henkin(") << eVar.getId();
        for (const Variable &depVar : DQBFPrefix.getExistVarDependencies(eVar)) {
            output << std::string(", ") << depVar.getId();
        }
        output << std::string(")") << std::endl;
    }

    // gates need to be copied, otherwise it can get a bit funky for some reason
    const auto gates = formulaPtr->formula.getGates();

    std::vector<unsigned long> hqspreVarToOutputVar = std::vector<unsigned long>(formulaPtr->formula.maxVarIndex() + 1);
    auto hqspreLitToStringVar = [&hqspreVarToOutputVar](const hqspre::Literal lit) {
                                std::string var = std::to_string(hqspreVarToOutputVar[hqspre::lit2var(lit)]);
                                return (hqspre::isNegative(lit) ? (std::string("-") + var) : var);
                                };

    // Load the hqspreVarToOutputVar with the primary inputs
    unsigned int maxQuantifiedVarIndex = 0; 
    for (auto var = formulaPtr->formula.minVarIndex(); var <= formulaPtr->formula.maxVarIndex(); ++var) {
        if (!formulaPtr->formula.varDeleted(var) && !formulaPtr->formula.isGateOutput(var)) {
            auto varId = formulaPtr->getDqbddVarMappedIntoHqspreVar(var).getId();
            hqspreVarToOutputVar[var] = varId;
            if (varId > maxQuantifiedVarIndex) {
                maxQuantifiedVarIndex = varId;
            }
        }
    }

    unsigned long outputVar = maxQuantifiedVarIndex + 1; 
    output << "output(" << outputVar << ")" << std::endl;

    // Print the gates
    unsigned long gateVarId = outputVar + 1; 
    for (const hqspre::Gate& g: gates) {
        const hqspre::Variable out_var = hqspre::lit2var(g._output_literal);

        switch (g._type) {
            case hqspre::GateType::AND_GATE:
            {
                output << gateVarId << " = and(";
                for (auto inputLitIter = g._input_literals.begin(); inputLitIter != g._input_literals.end(); ++inputLitIter) {
                    if (inputLitIter == g._input_literals.begin()) {
                        output << hqspreLitToStringVar(*inputLitIter);
                    } else {
                        output << ", " << hqspreLitToStringVar(*inputLitIter);
                    }
                }
                output << ")" << std::endl;
                break;
            }

            case hqspre::GateType::XOR_GATE:
            {
                auto inputLit1 = g._input_literals[0];
                unsigned long inputVar1 = hqspreVarToOutputVar[hqspre::lit2var(inputLit1)];
                auto inputLit2 = g._input_literals[1];
                unsigned long inputVar2 = hqspreVarToOutputVar[hqspre::lit2var(inputLit2)];
                if ((hqspre::isNegative(inputLit1) && hqspre::isNegative(inputLit2))
                        || (!hqspre::isNegative(inputLit1) && !hqspre::isNegative(inputLit2))) {
                    // A XOR B = !A XOR !B = (A AND !B) OR (!A AND B)
                    output << gateVarId << " = and(" << inputVar1 << ", -" << inputVar2 << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = and(-" << inputVar1 << ", " << inputVar2 << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = or(" << (gateVarId - 1) << ", " << (gateVarId - 2) << ")" << std::endl;
                } else {
                    // A XOR !B = !A XOR B = (A AND B) OR (!A AND !B)
                    output << gateVarId << " = and(" << inputVar1 << ", " << inputVar2 << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = and(-" << inputVar1 << ", -" << inputVar2 << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = or(" << (gateVarId - 1) << ", " << (gateVarId - 2) << ")" << std::endl;
                }
                break;
            }

            case hqspre::GateType::MUX_GATE:
            {
                auto litA = g._input_literals[0];
                unsigned long varA = hqspreVarToOutputVar[hqspre::lit2var(litA)];;
                std::string varB = hqspreLitToStringVar(g._input_literals[1]);
                std::string varC = hqspreLitToStringVar(g._input_literals[2]);
                if (hqspre::isNegative(litA)) {
                    // MUX(!A,B,C) = (!A AND B) OR (A AND C)
                    output << gateVarId << " = and(-" << varA << ", " << varB << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = and(" << varA << ", " << varC << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = or(" << (gateVarId - 1) << ", " << (gateVarId - 2) << ")" << std::endl;
                } else {
                    // MUX(A,B,C) = (A AND B) OR (!A AND C)
                    output << gateVarId << " = and(" << varA << ", " << varB << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = and(-" << varA << ", " << varC << ")" << std::endl;
                    ++gateVarId;
                    output << gateVarId << " = or(" << (gateVarId - 1) << ", " << (gateVarId - 2) << ")" << std::endl;
                }
                break;
            }

            default:
            {
                throw DQBDDexception("Invalid gate type encountered");
                break;
            }
        }

        if (hqspre::isNegative(g._output_literal)) {
            ++gateVarId;
            output << gateVarId << " = and(-" << (gateVarId - 1)  << ")" << std::endl;
        }

        hqspreVarToOutputVar[out_var] = gateVarId;
        ++gateVarId;

        for (const auto c_nr: g._encoding_clauses) {
            formulaPtr->formula.removeClause(c_nr);
        }
    }

    // Convert the clauses
    std::list<unsigned long> clauseGateIDs;
    for (hqspre::ClauseID c_nr = 0; c_nr <= formulaPtr->formula.maxClauseIndex(); c_nr++) {
        if (!formulaPtr->formula.clauseDeleted(c_nr)) {
            output << gateVarId << " = or(";
            const auto& clause = formulaPtr->formula.getClause(c_nr);
            for (auto inputLitIter = clause.begin(); inputLitIter != clause.end(); ++inputLitIter) {
                if (inputLitIter == clause.begin()) {
                    output << hqspreLitToStringVar(*inputLitIter);
                } else {
                    output << ", " << hqspreLitToStringVar(*inputLitIter);
                }
            }
            output << ")" << std::endl;
            clauseGateIDs.push_back(gateVarId);
            ++gateVarId;
        }
    }

    // print out the output gate
    output << outputVar << " = and(";
    for (auto clauseIDiter = clauseGateIDs.begin(); clauseIDiter != clauseGateIDs.end(); ++ clauseIDiter) {
        if (clauseIDiter == clauseGateIDs.begin()) {
            output << *clauseIDiter;
        } else {
            output << ", " << *clauseIDiter;
        }
    }
    output << ")" << std::endl;
}