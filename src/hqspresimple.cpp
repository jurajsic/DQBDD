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

// simpler HQSpre binary made for benchmarking

#include <easylogging++.hpp>
#include <formula.hpp>

INITIALIZE_EASYLOGGINGPP

int main(int, char **argv)
{
    /*
    Cudd mgr1;
    QuantifiedVariablesManager qvm1;
    std::vector<QuantifierTreeNode*> treeees = { };
    QuantifierTreeFormula *varFormula1 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula1->setMatrix(Variable(0, mgr1));
            treeees.push_back(varFormula1);
    QuantifierTreeFormula *varFormula12 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula12->setMatrix(!Variable(0, mgr1));
            treeees.push_back(varFormula12);
            QuantifierTreeFormula *varFormula2 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula2->setMatrix(Variable(1, mgr1));
            treeees.push_back(varFormula2);
            QuantifierTreeFormula *varFormula22 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula22->setMatrix(Variable(1, mgr1));
            treeees.push_back(varFormula22);
    QuantifierTreeNode *firstOp = varFormula1;
                QuantifierTreeNode *firstOpNeg = varFormula12;
                QuantifierTreeNode *secondOp = varFormula2;
                QuantifierTreeNode *secondOpNeg = varFormula22;
                QuantifierTreeNode *firstConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOp, secondOpNeg}, qvm1);
            treeees.push_back(firstConj);
                QuantifierTreeNode *secondConj = new QuantifierTree(true, std::list<QuantifierTreeNode*>{firstOpNeg, secondOp}, qvm1);
            treeees.push_back(secondConj);
    
            QuantifierTreeFormula *varFormula3 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula3->setMatrix(Variable(2, mgr1));
            treeees.push_back(varFormula3);
            QuantifierTreeFormula *varFormula32 = new QuantifierTreeFormula(mgr1, qvm1);
            varFormula32->setMatrix(!Variable(2, mgr1));
            treeees.push_back(varFormula32);
    QuantifierTreeNode *qtxor = new QuantifierTree(false, std::list<QuantifierTreeNode*>{firstConj, secondConj}, qvm1);
    QuantifierTreeNode *qtor = new QuantifierTree(false, std::list<QuantifierTreeNode*>{qtxor, varFormula3}, qvm1);
    //QuantifierTreeNode *qtor2 = new QuantifierTree(false, std::list<QuantifierTreeNode*>{qtor, varFormula32}, qvm1);
    //QuantifierTreeNode *qtor3 = new QuantifierTree(false, std::list<QuantifierTreeNode*>{qtxor, varFormula3, varFormula32}, qvm1);
    treeees.push_back(qtxor);
    treeees.push_back(qtor);
    //treeees.push_back(qtor2);
    //treeees.push_back(qtor3);
    //for (auto tree : treeees) {
    //    std::cout << *tree << " with support set " << tree->getSupportSet() <<std::endl;
    //}
    std::cout << *qtor << " with support set " << qtor->getSupportSet() <<std::endl;
    return 0;*/
                

    hqspre::Formula formula;

    // Configure logging
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Enabled, "true");
    defaultConf.setGlobally(el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.set(el::Level::Verbose, el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
    defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
    el::Loggers::reconfigureAllLoggers(defaultConf);
    el::Loggers::setVerboseLevel(0);

    formula.settings().consistency_check = false;

    // Parse the file
    std::string in_name(argv[1]);
    std::ifstream in(in_name);
    in >> formula;
    in.close();

    // do the preprocessing magic
    try {
        formula.settings().bla              = false;
        formula.settings().ble              = false;
        formula.settings().pure_sat_timeout = 1000;
        formula.preprocess();
        //formula.printStatistics();
    } catch (hqspre::SATException&) {
        std::cout << "p cnf 0 0\n" << std::endl;
        return 10;
    } catch (hqspre::UNSATException&) {
        std::cout << "p cnf 0 1\n0\n" << std::endl;
        return 20;
    }

    std::cout << formula;
    return 30;
}