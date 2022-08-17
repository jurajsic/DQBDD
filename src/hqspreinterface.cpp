/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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

#include <vector>
#include <sstream>

#include <easylogging++.hpp>
#include <formula.hpp>

#include "hqspreinterface.hpp"
#include "dqbddexceptions.hpp"


/**
 * Everything here is based on the way how preprocessing is done in HQS, especially
 * how the result of the preprocessing is tranformed into AIG.
 **/

namespace dqbdd {

HQSPreInterface::HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : GateParser(mgr, qvmgr) {}

HQSPreResult HQSPreInterface::getPreprocessorResult() {
    return result;
}

// code based on HQS 
void HQSPreInterface::parse(std::string fileName) {
    result = HQSPreResult::UNKNOWN;

    // TODO create custom logger for HQSpre
    // Configure logging
    const el::Configurations *defConf = el::Loggers::defaultConfigurations(); // keep default logger so after preprocessing we can restore it
    el::Configurations hqspreConf; // configurations used for HQSpre
    // we want it to have same stuff as defConf, except for log messages
    for (auto conf : defConf->list()) { // we have to do this in a bit dirty way, we cannot use hqspreConf.setFromBase(defConf) as defConf is const
        hqspreConf.set(conf);
    }
    hqspreConf.setGlobally(el::ConfigurationType::Format, "%datetime %level [HQSpre] %msg");
    hqspreConf.set(el::Level::Verbose, el::ConfigurationType::Format, "%datetime %level-%vlevel [HQSpre] %msg");
    el::Loggers::reconfigureLogger("default", hqspreConf);

    hqspre::Formula formula;
    formula.settings().consistency_check = false;

    // Parse the file
    std::ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        std::string errorMes = "Could not open file '";
        errorMes += fileName + "'.";
        throw dqbddException(errorMes);
    }

    try {
        inputFile >> formula;
        inputFile.close();

        // do the preprocessing magic
        formula.determineGates();
        auto numOfDetectedGates = formula.getGates().size();
        VLOG(1) << "Detected " << numOfDetectedGates << " gates in the input formula";
        if (numOfDetectedGates > 5) { 
            /* If we have more than 5 gates, we would like to preserve them, but first we would like to
             * try HQSpre to solve the formula without preserving them
             */

            // First do full preprocessing (without preserving gates) on a copy of the formula
            hqspre::Formula formula2(formula);
            // TODO decide what settings to use (here are except timeout and forking same as in HQS 18-03-2021)
            // bla and ble are only useful for QBF (I think)
            formula2.settings().bla                     = false;
            formula2.settings().ble                     = false;
            formula2.settings().impl_chains             = 3;
            formula2.settings().max_substitution_cost   = 250;
            formula2.settings().max_resolution_cost     = 100;
            formula2.settings().vivify_fp               = true;
            // TODO timeout?
            //formula2.settings().pure_sat_timeout = 1000;
            // in HQS forksplitting is enabled, but it uses external MPhaseSAT64 solver, so here we better turn it off
            formula2.settings().enableFork       = false;
            VLOG(1) << "Start preprocessing without gate preservation";
            formula2.preprocess();

            // Then do preprocessing, preserving gates
            // TODO setting??? (here are except timeout and forking same as in HQS 18-03-2021)
            formula.settings().univ_expand      = 0; // maybe 2 is better???
            // bla and ble are only useful for QBF (I think)
            formula.settings().bla              = false;
            formula.settings().ble              = false;
            formula.settings().preserve_gates   = true;
            formula.settings().substitution     = false;
            formula.settings().rewrite          = false;
            formula.settings().resolution       = false;
            formula.settings().max_loops        = 20;
            formula.settings().enableFork       = false;
            // TODO timeout??
            //formulaPtr->formula.settings().pure_sat_timeout = 1000;

            VLOG(1) << "Start preprocessing with gate preservation";
        } else { // if we have no gates, there is no point in preserving them, run HQSpre normally
            // TODO decide what settings to use (here are except timeout and forking same as in HQS 18-03-2021)
            // bla and ble are only useful for QBF (I think)
            formula.settings().bla                      = false;
            formula.settings().ble                      = false;
            formula.settings().impl_chains              = 3;
            formula.settings().max_substitution_cost    = 250;
            formula.settings().max_resolution_cost      = 100;
            formula.settings().vivify_fp                = true;
            // TODO timeout?
            //formulaPtr->formula.settings().pure_sat_timeout = 1000;
            // in HQS forking is enabled, but it uses external MPhaseSAT64 solver, so here we better turn it off
            formula.settings().enableFork               = false;

            VLOG(1) << "Start preprocessing without gate preservation";
        }
        formula.preprocess();
        if (el::Loggers::verboseLevel() > 0) {
            formula.printStatistics();
        }
    } catch (hqspre::SATException&) {
        finishedParsing(true, addGate(GateType::AND));
        result = HQSPreResult::SAT;
        return;
    } catch (hqspre::UNSATException&) {
        finishedParsing(true, addGate(GateType::OR));
        result = HQSPreResult::UNSAT;
        return;
    }

    // update gates (if they were not preserved they might have gotten deleted)
    formula.determineGates();
        
    // enforce DQBF (to have DQBF prefix, it could have possibly changed to QBF)
    formula.enforceDQBF(true);
    
    // gates need to be copied, otherwise it can get a bit funky for some reason
    const auto gates = formula.getGates();

    // return back to previous configuration for logging
    el::Loggers::reconfigureLogger("default", *defConf);

    // mapping from hqspre (normal and gate) variables to vars in GateParser
    std::vector<unsigned long> hqspreVarToGateParserVar(formula.maxVarIndex() + 1);
    unsigned long gateParserMaxVar = 0; // TODO use maxGateID of parent class + right now variables are indexed from 1, so there is one not used variable 0 in CUDD
    auto getNewGateParserVar = [&gateParserMaxVar]() {
        ++gateParserMaxVar;
        return gateParserMaxVar;
    };

    // Create the proper problem variables (without Tseitin variables)
    // first universal...
    unsigned long numOfUnivVars = 0;
    for (hqspre::Variable var = formula.minVarIndex(); var <= formula.maxVarIndex(); ++var) {
        if (!formula.varDeleted(var) && !formula.isGateOutput(var) && formula.isUniversal(var)) {
            hqspreVarToGateParserVar[var] = getNewGateParserVar();
            addUnivVar(hqspreVarToGateParserVar[var]);
            ++numOfUnivVars;
        }
    }
    // ...and then existential (so we can add dependencies without worrying whether univ var to be added to dependency was created or not)
    // also the order of variables is important for creating BDDs, the default HQSpre order makes DQBDD really slow
    unsigned long numOfExistVars = 0;
    for (hqspre::Variable var = formula.minVarIndex(); var <= formula.maxVarIndex(); ++var) {
        if (!formula.varDeleted(var) && !formula.isGateOutput(var) && formula.isExistential(var)) {
            hqspreVarToGateParserVar[var] = getNewGateParserVar();
            std::vector<unsigned long> varDependencies;
            for (auto dep : formula.getDependencies(var)) {
                varDependencies.push_back(hqspreVarToGateParserVar[dep]);
            }
            addExistVar(hqspreVarToGateParserVar[var], varDependencies);
            ++numOfExistVars;
        }
    }

    VLOG(1) << "Preprocessing finished with formula with " << numOfExistVars << " existential vars, " << numOfUnivVars << " univ. vars, "
            << formula.getGates().size() << " gates, and " << formula.numClauses() << " clauses which potentionally encode the gates.";


    /* We will save in isHqspreVarOutputNegated the info about whether output of hqspre gates are negated.
     * The default value is false, as the output of normal variables are not negated.
     */
    std::vector<bool> isHqspreVarOutputNegated(formula.maxVarIndex() + 1, false);

    auto getGateParserLiteralFromHqspreInputLiteral = [&isHqspreVarOutputNegated, &hqspreVarToGateParserVar](const hqspre::Literal &inputLiteral)->GateLiteral {
        hqspre::Variable inputLitVar = hqspre::lit2var(inputLiteral);

        bool isInputNegated = hqspre::isNegative(inputLiteral);
        bool isOutputGateRepresentedByInputLitVarNegated = isHqspreVarOutputNegated[inputLitVar];
        
        if ( (isInputNegated && isOutputGateRepresentedByInputLitVarNegated) || // if both the input literal and the output of the gate represented by the inputLitVar are either both negated or...
            (!isInputNegated && !isOutputGateRepresentedByInputLitVarNegated))  // ...are both not negated then...
        {
            // ...the resulting input gate should not be negated
            return GateLiteral(true, hqspreVarToGateParserVar[inputLitVar]);
        } else {
            // otherwise it should be negated
            return GateLiteral(false, hqspreVarToGateParserVar[inputLitVar]);
        }
    };

    // we process the hqspre gates
    for (const hqspre::Gate &g : gates) {

        GateType newGateType;
        switch (g._type) {
            case hqspre::GateType::AND_GATE:
            {
                newGateType = GateType::AND;
                break;
            }

            case hqspre::GateType::XOR_GATE:
            {
                newGateType = GateType::XOR;
                break;
            }

            case hqspre::GateType::MUX_GATE:
            {
                newGateType = GateType::MUX;
                break;
            }

            default:
            {
                throw dqbddException("Invalid gate type encountered during processing output of HQSpre");
                break;
            }
        }

        std::vector<GateLiteral> newGateOperands;
        for (auto inputLitIter = g._input_literals.begin(); inputLitIter != g._input_literals.end(); ++inputLitIter) {
            newGateOperands.push_back(getGateParserLiteralFromHqspreInputLiteral(*inputLitIter));
        }


        const hqspre::Variable out_var = hqspre::lit2var(g._output_literal);
        unsigned long newGateID = getNewGateParserVar();
        hqspreVarToGateParserVar[out_var] = newGateID;
        isHqspreVarOutputNegated[out_var] = hqspre::isNegative(g._output_literal);

        addGate(newGateID, newGateType, newGateOperands);

        // remove clauses from which the hqspre gate is formed
        for (const auto c_nr: g._encoding_clauses) {
            formula.removeClause(c_nr);
        }
    }

    // now we convert the clauses
    std::vector<GateLiteral> clauseGateIDs;
    for (hqspre::ClauseID c_nr = 0; c_nr <= formula.maxClauseIndex(); c_nr++) {
        if (!formula.clauseDeleted(c_nr)) {
            
            const auto& clause = formula.getClause(c_nr);

            std::vector<GateLiteral> clauseLiterals;
            for (auto inputLitIter = clause.begin(); inputLitIter != clause.end(); ++inputLitIter) {
                clauseLiterals.push_back(getGateParserLiteralFromHqspreInputLiteral(*inputLitIter));
            }

            clauseGateIDs.push_back(GateLiteral(true, addGate(GateType::OR, clauseLiterals)));
        }
    }
    
    finishedParsing(true, addGate(GateType::AND, clauseGateIDs));

    return;
}


} // namespace dqbdd